from pathlib import Path
import numpy as np
import pandas as pd
import torch
from torch import nn
import matplotlib.pyplot as plt
import Entrenamiento_red as Mo

SEMILLAS = (1, 2, 3)
# Penalización de rugosidad fija para todas las ablaciones
LAMBDA = 100.0
EPOCAS, PACIENCIA, TASA = 3000, 200, 0.01
# Pentadas de cada temporada de lluvias para ubicar los dos picos (marzo-junio y septiembre-diciembre)
TEMPORADAS = ((13, 36), (49, 72))
# Resultados intermedios: cada etapa guarda un CSV y los modelos finales un .pt;
# si el archivo existe se carga en lugar de entrenar (borrar el archivo para repetir la etapa)
RESULTADOS = Path('Resultados_Ablacion_arquitectura')
CABEZAS = ('lineal', 'MLP', 'KAN')
COLORES = {'lineal': 'tab:red', 'MLP': 'tab:blue', 'KAN': 'tab:green', 'Satélite': 'darkorange'}
T = lambda x: torch.as_tensor(x, dtype=torch.float64)


# B-spline cúbica con nodos equiespaciados (recurrencia de Cox-de Boor)
def bspline(t, a, b, K, orden=4):
    cortes = np.linspace(a, b, K - orden + 2)
    nudos = np.r_[[a] * (orden - 1), cortes, [b] * (orden - 1)]
    t = np.asarray(t, float)
    B = np.array([(t >= nudos[j]) & (t < nudos[j + 1]) for j in range(len(nudos) - 1)], float).T
    B[t == b, np.max(np.nonzero(nudos[:-1] < nudos[1:])[0])] = 1
    for k in range(1, orden):
        N = np.zeros((len(t), len(nudos) - 1 - k))
        for j in range(len(nudos) - 1 - k):
            d1, d2 = nudos[j + k] - nudos[j], nudos[j + k + 1] - nudos[j + 1]
            if d1 > 0: N[:, j] += (t - nudos[j]) / d1 * B[:, j]
            if d2 > 0: N[:, j] += (nudos[j + k + 1] - t) / d2 * B[:, j + 1]
        B = N
    return B


# Datos de la red (raíz después de suavizar), lluvia observada, formas del ONI y referencias
def preparar_datos(carpeta_curvas='Curvas', carpeta_datos='data'):
    carpeta_curvas, carpeta_datos = Path(carpeta_curvas), Path(carpeta_datos)
    datos = Mo.cargar_datos_red(carpeta_curvas, 'raiz')
    meta = datos['meta']
    tr, va = datos['particion'] == 'entrenamiento', datos['particion'] == 'validacion'
    t_fina, filas = datos['t'], datos['filas_pentadas']
    D = dict(datos=datos, meta=meta, tr=tr, va=va, t_fina=t_fina, filas=filas,
             phi_s_f=datos['phi_sat'], phi_e_f=datos['phi_est'],
             phi_s=datos['phi_sat'][filas], phi_e=datos['phi_est'][filas])
    D['W_s'], D['R_s'] = Mo.matrices_integracion(t_fina, D['phi_s_f'])
    D['W_e'], _ = Mo.matrices_integracion(t_fina, D['phi_e_f'])

    leer = lambda f: pd.read_csv(carpeta_curvas / f, dtype={'curva': str}).set_index('curva').loc[meta.curva].to_numpy(float)
    D['observado'] = leer('observado_estaciones.csv')[va]
    D['satelite_obs'] = leer('observado_satelite.csv')[va]
    D['estacion_mm'] = leer('coeficientes_estaciones.csv')[va] @ D['phi_e_f'].T
    D['satelite_mm'] = leer('coeficientes_satelite.csv')[va] @ D['phi_s_f'].T
    D['climatologia'] = np.nanmean(leer('observado_estaciones.csv')[tr], axis=0)
    D['curvas_va'] = meta.curva[va].to_numpy()

    # ONI: promedio anual, máximo del año, ventana DEF y curva de las 12 ventanas en [0, 72]
    enso = pd.read_csv(carpeta_datos / 'enso_evaluacion.csv')
    tabla = enso.pivot(index='YEAR', columns='mes_central', values='ONI')
    D['oni_meses'] = tabla.loc[meta.YEAR].to_numpy(float)
    D['t_meses'] = 6 * np.arange(1, 13) - 3
    D['phi_o_f'] = bspline(t_fina, 0, 72, 6)
    D['c_oni'] = np.linalg.lstsq(bspline(D['t_meses'], 0, 72, 6), D['oni_meses'].T, rcond=None)[0].T
    D['W_o'], D['R_o'] = Mo.matrices_integracion(t_fina, D['phi_o_f'])
    D['fase_def'] = enso[enso.mes_central == 1].set_index('YEAR').fase_enso.loc[meta.YEAR[va]].to_numpy()
    estandarizar = lambda x: (x - x[tr].mean(0)) / x[tr].std(0)
    D['altitud'] = meta.altitud_completa.to_numpy(float)
    D['covariables'] = {
        'promedio anual': estandarizar(np.c_[D['altitud'], meta.oni_promedio.to_numpy(float)]),
        'máximo': estandarizar(np.c_[D['altitud'], meta.oni_max_intensidad.to_numpy(float)]),
        'DEF': estandarizar(np.c_[D['altitud'], tabla.loc[meta.YEAR, 1].to_numpy(float)]),
        'curva': estandarizar(D['altitud'][:, None]),
    }
    # Predicciones de las regresiones funcionales del script de R
    reg = pd.read_csv(carpeta_curvas / 'regresion_funcional_predicciones.csv', dtype={'curva': str})
    D['regresion'] = {m: reg[reg.modelo == m].set_index('curva').loc[D['curvas_va']].drop(columns='modelo').to_numpy(float)
                      for m in reg.modelo.unique()}
    RESULTADOS.mkdir(exist_ok=True)
    return D


# Neuronas funcionales z = A W c, con penalización de rugosidad de las w_k(t)
class Codificador(nn.Module):
    def __init__(self, gram, rugosidad, q):
        super().__init__()
        self.register_buffer('gram', T(gram)); self.register_buffer('R', T(rugosidad))
        self.A = nn.Parameter(torch.empty(q, gram.shape[0], dtype=torch.float64)); nn.init.xavier_uniform_(self.A)
    def forward(self, c): return (c @ self.gram) @ self.A.T
    def rugosidad(self): return ((self.A @ self.R) * self.A).sum()


# Red: codificador del satélite (+ codificador del ONI si es curva) + covariables escalares -> cabeza -> M puntajes
# La KAN admite varias capas ocultas (anchos) y cualquier número de intervalos
class Red(nn.Module):
    def __init__(self, D, cabeza, q, m, n_escalares, anchos=(2,), intervalos=12, q_oni=0, ocultas_mlp=None):
        super().__init__()
        self.sat = Codificador(D['W_s'], D['R_s'], q)
        self.oni = Codificador(D['W_o'], D['R_o'], q_oni) if q_oni else None
        n = q + q_oni + n_escalares
        self.es_kan = cabeza == 'KAN'
        if cabeza == 'lineal':
            self.cabeza = nn.Linear(n, m, dtype=torch.float64)
        elif cabeza == 'MLP':
            self.cabeza = nn.Sequential(nn.Linear(n, ocultas_mlp), nn.Tanh(), nn.Linear(ocultas_mlp, m)).double()
        else:
            tam = [n, *anchos, m]
            self.cabeza = nn.ModuleList([Mo.CapaKAN(a, b, intervalos) for a, b in zip(tam[:-1], tam[1:])])
    def entradas(self, c, co, x):
        return torch.cat([self.sat(c)] + ([self.oni(co)] if self.oni is not None else []) + [x], dim=1)
    def forward(self, c, co, x):
        h = self.entradas(c, co, x)
        if not self.es_kan: return self.cabeza(h)
        for capa in self.cabeza: h = capa(h)
        return h
    @torch.no_grad()
    def calibrar(self, c, co, x):
        h = self.entradas(c, co, x)
        for capa in self.cabeza: capa.calibrar(h); h = capa(h)
    def penalizacion_rugosidad(self):
        return self.sat.rugosidad() + (self.oni.rugosidad() if self.oni is not None else 0)


# Neuronas ocultas del MLP para igualar los parámetros de la KAN con la misma entrada y salida
def ocultas_igualadas(n, m, anchos=(2,), intervalos=12):
    n_kan = sum(p.numel() for a, b in zip([n, *anchos], [*anchos, m]) for p in Mo.CapaKAN(a, b, intervalos).parameters())
    return max(1, round((n_kan - m) / (n + m + 1)))

def construir_red(D, cabeza, q, m, oni='promedio anual', anchos=(2,), intervalos=12, q_oni=2):
    n_esc = D['covariables'][oni].shape[1]; qo = q_oni if oni == 'curva' else 0
    ocultas = ocultas_igualadas(q + qo + n_esc, m, anchos, intervalos) if cabeza == 'MLP' else None
    return Red(D, cabeza, q, m, n_esc, anchos, intervalos, qo, ocultas)

# FPCA de las curvas de estaciones de entrenamiento, con puntajes escalados por su desviación
def fpca_estaciones(D, m):
    curvas = D['datos']['c_est'] @ D['phi_e'].T
    fp = Mo.fit_fpca(curvas[D['tr']], m)
    sd = ((curvas[D['tr']] - fp['media']) @ fp['base'].T).std(0)
    return fp, dict(fp, base=fp['base'] * sd[:, None]), sd

def argumentos(D, idx, oni):
    return T(D['datos']['c_sat'][idx]), T(D['c_oni'][idx]), T(D['covariables'][oni][idx])

# Entrenamiento con lote completo, AdamW y parada temprana con la pérdida L² de validación
def entrenar(D, cabeza, q, m, semilla, oni='promedio anual', anchos=(2,), intervalos=12, q_oni=2):
    tr, va = D['tr'], D['va']
    c_obj = D['datos']['c_est']
    fp, fp_n, _ = fpca_estaciones(D, m)
    mu_c, psi_c, _ = Mo.coeficientes_fpca(fp_n, D['phi_e'])
    Mo.establecer_semilla(semilla)
    red = construir_red(D, cabeza, q, m, oni, anchos, intervalos, q_oni)
    if cabeza == 'KAN': red.calibrar(*argumentos(D, tr, oni))
    perdida = Mo.PerdidaFuncional(D['W_e'], mu_c, psi_c, LAMBDA)
    opt = torch.optim.AdamW(red.parameters(), lr=TASA, weight_decay=0)
    y_tr, y_va = T(c_obj[tr]), T(c_obj[va])
    mejor, ep_m, estado, sin = np.inf, 0, None, 0
    for ep in range(1, EPOCAS + 1):
        opt.zero_grad(set_to_none=True)
        perdida(red(*argumentos(D, tr, oni)), y_tr, red).backward(); opt.step()
        with torch.no_grad(): v = perdida.error_l2(red(*argumentos(D, va, oni)), y_va).item()
        if v < mejor - 1e-10:
            mejor, ep_m, sin = v, ep, 0; estado = {k: w.clone() for k, w in red.state_dict().items()}
        else:
            sin += 1
            if sin >= PACIENCIA: break
    red.load_state_dict(estado)
    with torch.no_grad(): c_hat = perdida.coeficientes(red(*argumentos(D, va, oni))).numpy()
    return red, c_hat, dict(L2_validacion=mejor, epoca=ep_m, parametros_cabeza=sum(p.numel() for p in red.cabeza.parameters()),
                            varianza_explicada=fp['explicada'][:m].sum())

# Coeficientes en la escala del modelo -> curva en mm (recortada en 0) en las 72 pentadas y en la malla fina
def a_mm(D, c_hat):
    return np.maximum(c_hat @ D['phi_e'].T, 0)**2, np.maximum(c_hat @ D['phi_e_f'].T, 0)**2


# Métricas en validación: RECM y sesgo frente a lo observado; L² relativo, KGE y fecha de los picos frente a la curva de la estación
def metricas(D, pred_pent, pred_fina=None):
    obs = D['observado']; m = ~np.isnan(obs)
    out = {'RECM_mm': np.sqrt(np.mean((pred_pent[m] - obs[m])**2)), 'sesgo_%': 100 * (pred_pent[m].sum() / obs[m].sum() - 1)}
    if pred_fina is not None:
        t, ref_f = D['t_fina'], D['estacion_mm']
        out['L2_relativo'] = np.mean(np.sqrt(np.trapezoid((pred_fina - ref_f)**2, t, axis=1) / np.trapezoid(ref_f**2, t, axis=1)))
        ref = ref_f[:, D['filas']]
        r = np.array([np.corrcoef(a, b)[0, 1] for a, b in zip(pred_pent, ref)])
        kge = 1 - np.sqrt((r - 1)**2 + (pred_pent.std(1) / ref.std(1) - 1)**2 + (pred_pent.mean(1) / ref.mean(1) - 1)**2)
        out['KGE'] = np.median(kge)
        for j, (a, b) in enumerate(TEMPORADAS, 1):
            v = (t >= a) & (t <= b)
            out[f'pico{j}_pentadas'] = np.mean(np.abs(t[v][pred_fina[:, v].argmax(1)] - t[v][ref_f[:, v].argmax(1)]))
    return out

def referencias(D):
    filas = {'Satélite sin corregir (péntadas)': metricas(D, D['satelite_obs']),
             'Satélite sin corregir (curva)': metricas(D, D['satelite_mm'][:, D['filas']], D['satelite_mm']),
             'Climatología de entrenamiento': metricas(D, np.tile(D['climatologia'], (D['va'].sum(), 1)))}
    for m, p in D['regresion'].items(): filas[f'Regresión: {m}'] = metricas(D, p)
    return pd.DataFrame(filas).T


# Cambio relativo de los coeficientes del codificador respecto a su inicialización y fracción de entradas
# de entrenamiento fuera del soporte de las B-splines (promedio de las capas de la KAN)
@torch.no_grad()
def diagnostico_codificador(D, red, cabeza, q, m, semilla, oni='promedio anual', anchos=(2,), intervalos=12):
    Mo.establecer_semilla(semilla)
    A0 = construir_red(D, cabeza, q, m, oni, anchos, intervalos).sat.A
    cambio = float((red.sat.A - A0).norm() / A0.norm())
    if not red.es_kan: return cambio, np.nan
    h, fuera = red.entradas(*argumentos(D, D['tr'], oni)), []
    for capa in red.cabeza:
        u = (h - capa.centro) / capa.escala
        fuera.append(float((u.abs() > 1 + 6 / capa.intervalos).double().mean())); h = capa(h)
    return cambio, float(np.mean(fuera))

# Arquitecturas candidatas de la KAN con 3 semillas: métricas, cambio del codificador y entradas fuera del soporte
def diagnostico_arquitecturas(D, registros, candidatos, q, m):
    for anchos, intervalos in candidatos:
        for s in SEMILLAS:
            red, c_hat, info = entrenar(D, 'KAN', q, m, s, anchos=anchos, intervalos=intervalos)
            cambio, fuera = diagnostico_codificador(D, red, 'KAN', q, m, s, anchos=anchos, intervalos=intervalos)
            registros.append(dict(anchos=str(anchos), intervalos=intervalos, semilla=s, **info, cambio_codificador=cambio,
                                  fuera_soporte=fuera, **metricas(D, *a_mm(D, c_hat))))

# Persistencia: una configuración entrenada con las 3 semillas se agrega a 'registros'
def correr(D, registros, cabeza, q, m, **kw):
    for s in SEMILLAS:
        _, c_hat, info = entrenar(D, cabeza, q, m, s, **kw)
        registros.append(dict(cabeza=cabeza, q=q, M=m, semilla=s, **{k: (str(v) if isinstance(v, tuple) else v) for k, v in kw.items()},
                              **info, **metricas(D, *a_mm(D, c_hat))))

# Etapa persistente: si su CSV existe se carga; si no, se entrena con 'funcion(registros)' y se guarda
def etapa(nombre, funcion):
    archivo = RESULTADOS / f'{nombre}.csv'
    if archivo.exists(): return pd.read_csv(archivo)
    registros = []; funcion(registros)
    df = pd.DataFrame(registros); df.to_csv(archivo, index=False); return df

# Modelo persistente: pesos y predicción de validación en un .pt
def red_guardada(D, nombre, cabeza, q, m, semilla, **kw):
    archivo = RESULTADOS / f'{nombre}.pt'
    if archivo.exists():
        d = torch.load(archivo, weights_only=False)
        red = construir_red(D, cabeza, q, m, **kw); red.load_state_dict(d['estado']); return red, d['c_hat']
    red, c_hat, _ = entrenar(D, cabeza, q, m, semilla, **kw)
    torch.save({'estado': red.state_dict(), 'c_hat': c_hat}, archivo); return red, c_hat

# Media de 3 semillas de la configuración: modelo de la primera semilla y predicción promedio en mm
def modelo_final(D, nombre, cabeza, cfg, oni):
    cfg = dict(cfg); q, m = cfg.pop('q'), cfg.pop('M')
    preds, redes = [], []
    for s in SEMILLAS:
        red, c_hat = red_guardada(D, f'{nombre}_semilla{s}', cabeza, q, m, s, oni=oni, **cfg)
        preds.append(a_mm(D, c_hat)); redes.append(red)
    return redes[0], (np.mean([p[0] for p in preds], 0), np.mean([p[1] for p in preds], 0))

def resumen(df, por, cols=('RECM_mm', 'sesgo_%', 'L2_validacion')):
    g = df.groupby(list(por))
    r = g[list(cols)].mean()
    r['RECM_sd'] = g['RECM_mm'].std()
    return r.round(2)

def texto_anchos(texto):
    return tuple(int(a) for a in str(texto).strip('()').split(',') if a.strip())


# Gráficas
def grafica_oni_curva(D, anios=(1997, 1998, 1999)):
    fig, ax = plt.subplots(figsize=(7, 3.2))
    for y in anios:
        i = np.flatnonzero(D['meta'].YEAR.to_numpy() == y)[0]
        linea, = ax.plot(D['t_fina'], D['phi_o_f'] @ D['c_oni'][i], label=str(y))
        ax.plot(D['t_meses'], D['oni_meses'][i], 'o', ms=3, color=linea.get_color())
    ax.axhline(0, color='k', lw=.5); ax.set(xlabel='Péntada', ylabel='ONI', title='ONI como curva'); ax.legend(frameon=False)
    return fig

def grafica_q_M(e1, referencia, limites=(27.5, 30.5)):
    fig, axes = plt.subplots(1, 3, figsize=(14, 3.8), sharey=True)
    for ax, cabeza in zip(axes, CABEZAS):
        tabla = e1[e1.cabeza == cabeza].groupby(['q', 'M']).RECM_mm.mean().unstack()
        im = ax.imshow(tabla, cmap='viridis_r', vmin=limites[0], vmax=limites[1], aspect='auto')
        ax.set(xticks=range(tabla.shape[1]), xticklabels=tabla.columns, yticks=range(tabla.shape[0]),
               yticklabels=tabla.index, xlabel='M', title=cabeza)
        for i in range(tabla.shape[0]):
            for j in range(tabla.shape[1]): ax.text(j, i, f'{tabla.iat[i, j]:.2f}', ha='center', va='center', fontsize=8, color='w')
    axes[0].set_ylabel('q')
    fig.colorbar(im, ax=axes, label=f'RECM en validación (mm); satélite = {referencia:.2f}')
    return fig

def grafica_arquitectura_kan(e2, referencia, titulo):
    fig, ax = plt.subplots(figsize=(9, 3.6))
    for iv, g in e2.groupby('intervalos'):
        s = g.groupby('anchos', sort=False).RECM_mm.agg(['mean', 'std'])
        ax.errorbar(range(len(s)), s['mean'], s['std'], marker='o', capsize=3, label=f'{int(iv)} intervalos')
    ax.set_xticks(range(len(s))); ax.set_xticklabels(s.index)
    ax.axhline(referencia, ls='--', c='grey', label='Satélite sin corregir')
    ax.set(xlabel='Capas ocultas de la KAN', ylabel='RECM validación (mm)', title=titulo); ax.legend(frameon=False)
    return fig

def grafica_pesos_oni(D, red):
    w = D['phi_o_f'] @ red.oni.A.detach().numpy().T
    fig, ax = plt.subplots(figsize=(7, 3.2))
    for k in range(w.shape[1]): ax.plot(D['t_fina'], w[:, k], lw=2, label=f'$w^{{ONI}}_{k+1}(t)$')
    ax.axhline(0, c='k', lw=.5)
    ax.set(xlabel='Péntada (mes central de la ventana ONI)', ylabel='peso', title='Funciones de pesos del ONI (KAN)')
    ax.legend(frameon=False)
    return fig

def grafica_media_sesgo(D, predicciones):
    t, obs = np.arange(1, 73), D['observado']; m_obs = np.nanmean(obs, 0)
    fig, axes = plt.subplots(1, 2, figsize=(14, 4))
    axes[0].plot(t, m_obs, 'k', lw=2, label='Observado')
    axes[0].plot(t, np.nanmean(D['satelite_obs'], 0), c=COLORES['Satélite'], label='Satélite')
    for nombre, p in predicciones.items():
        axes[0].plot(t, p[0].mean(0), c=COLORES.get(nombre), label=nombre)
        axes[1].plot(t, 100 * (p[0].mean(0) / m_obs - 1), c=COLORES.get(nombre), label=nombre)
    axes[1].plot(t, 100 * (np.nanmean(D['satelite_obs'], 0) / m_obs - 1), c=COLORES['Satélite'], label='Satélite')
    axes[0].set(xlabel='Péntada', ylabel='mm', title='Media en validación'); axes[0].legend(frameon=False)
    axes[1].axhline(0, c='k', lw=.5); axes[1].set(xlabel='Péntada', ylabel='%', title='Sesgo por péntada'); axes[1].legend(frameon=False)
    return fig

def grafica_casos(D, predicciones, referencia='KAN'):
    t = np.arange(1, 73)
    err = np.sqrt(np.trapezoid((predicciones[referencia][1] - D['estacion_mm'])**2, D['t_fina'], axis=1))
    casos = {'Mejor': np.argmin(err), 'Típica': np.argsort(err)[len(err) // 2], 'Peor': np.argmax(err)}
    fig, axes = plt.subplots(1, 3, figsize=(15, 3.8))
    for ax, (nombre, i) in zip(axes, casos.items()):
        ax.plot(t, D['observado'][i], 'o', ms=3, c='grey', label='Observado')
        ax.plot(D['t_fina'], D['estacion_mm'][i], 'k', lw=2, label='Estación (MCP)')
        ax.plot(D['t_fina'], D['satelite_mm'][i], c=COLORES['Satélite'], label='Satélite')
        for c, p in predicciones.items(): ax.plot(D['t_fina'], p[1][i], c=COLORES.get(c), label=c)
        ax.set(title=f'{nombre}: {D["curvas_va"][i]}', xlabel='Péntada', ylabel='mm')
    axes[0].legend(frameon=False, fontsize=7)
    return fig

def grafica_pesos_satelite(D, red):
    w = D['phi_s_f'] @ red.sat.A.detach().numpy().T
    media = (D['datos']['c_sat'][D['tr']] @ D['phi_s_f'].T).mean(0)
    q = w.shape[1]; ncol = 4; nfil = int(np.ceil(q / ncol))
    fig, axes = plt.subplots(nfil, ncol, figsize=(14, 2.6 * nfil), squeeze=False)
    for k, ax in enumerate(axes.flat):
        if k >= q: ax.axis('off'); continue
        ax.fill_between(D['t_fina'], 0, media / media.max() * np.abs(w[:, k]).max(), color='0.85')
        ax.plot(D['t_fina'], w[:, k], c='tab:red', lw=2); ax.axhline(0, c='k', lw=.5)
        ax.set_title(f'$w_{{{k+1}}}(t)$', fontsize=10); ax.set_xticks([0, 18, 36, 54, 72])
    fig.supxlabel('Péntada'); fig.tight_layout()
    return fig

def grafica_fpca(D, m):
    fp, _, _ = fpca_estaciones(D, m)
    fig, ax = plt.subplots(figsize=(7, 3.4))
    for k in range(m): ax.plot(np.arange(1, 73), fp['base'][k], label=f'$\\psi_{k+1}$')
    ax.axhline(0, c='k', lw=.5); ax.set(title='Funciones propias de la FPCA (escala √P)', xlabel='Péntada'); ax.legend(frameon=False, fontsize=8)
    return fig

def nombres_entradas(D, red, oni):
    nombres = [f'z{k+1}' for k in range(red.sat.A.shape[0])]
    if red.oni is not None: nombres += [f'oni{k+1}' for k in range(red.oni.A.shape[0])]
    return nombres + (['altitud', 'ONI'] if D['covariables'][oni].shape[1] == 2 else ['altitud'])

@torch.no_grad()
def graficas_activaciones_kan(D, red, oni):
    nombres = nombres_entradas(D, red, oni)
    h = red.entradas(*argumentos(D, D['tr'], oni))
    figuras = []
    for n_capa, capa in enumerate(red.cabeza):
        fig, axes = plt.subplots(capa.salidas, capa.entradas, figsize=(1.9 * capa.entradas, 1.6 * capa.salidas), squeeze=False)
        for o in range(capa.salidas):
            for e in range(capa.entradas):
                x = torch.linspace(h[:, e].min(), h[:, e].max(), 200, dtype=torch.float64)
                ax = axes[o, e]; ax.plot(x, capa.arista(o, e, x), lw=1.5); ax.tick_params(labelsize=6)
                if o == 0: ax.set_title(nombres[e] if n_capa == 0 else f'oculto {e+1}', fontsize=8)
                if e == 0: ax.set_ylabel(f'salida {o+1}', fontsize=8)
        fig.suptitle(f'KAN, capa {n_capa + 1}', fontsize=10); fig.tight_layout()
        figuras.append(fig)
        h = capa(h)
    return figuras

# Amplitud de cada arista de la primera capa: rango de phi(x) sobre los datos de entrenamiento
@torch.no_grad()
def amplitud_aristas(D, red, oni):
    h = red.entradas(*argumentos(D, D['tr'], oni)); capa = red.cabeza[0]
    amp = [[float(np.ptp(capa.arista(o, e, h[:, e]).numpy())) for e in range(capa.entradas)] for o in range(capa.salidas)]
    return pd.DataFrame(amp, index=[f'oculto {o+1}' for o in range(capa.salidas)], columns=nombres_entradas(D, red, oni)).round(2)

@torch.no_grad()
def grafica_puntajes(D, red, m, oni):
    fp, _, sd = fpca_estaciones(D, m)
    curvas = D['datos']['c_est'] @ D['phi_e'].T
    real = ((curvas[D['va']] - fp['media']) @ fp['base'].T) / sd
    pred = red(*argumentos(D, D['va'], oni)).numpy()
    fig, axes = plt.subplots(1, m, figsize=(3 * m, 3), squeeze=False)
    for j, ax in enumerate(axes[0]):
        ax.plot(real[:, j], pred[:, j], '.', ms=2, alpha=.5); lim = np.abs(real[:, j]).max()
        ax.plot([-lim, lim], [-lim, lim], 'k--', lw=.8)
        ax.set(title=f'$s_{j+1}$: r = {np.corrcoef(real[:, j], pred[:, j])[0, 1]:.2f}', xlabel='real', ylabel='predicho' if j == 0 else '')
    fig.tight_layout()
    return fig

def grafica_error_pentada(D, predicciones, referencia='KAN'):
    t, obs = np.arange(1, 73), D['observado']
    fig, axes = plt.subplots(1, 2, figsize=(14, 4))
    for c, p in predicciones.items():
        axes[0].plot(t, np.sqrt(np.nanmean((p[0] - obs)**2, 0)), c=COLORES.get(c), label=c)
    axes[0].plot(t, np.sqrt(np.nanmean((D['satelite_obs'] - obs)**2, 0)), c=COLORES['Satélite'], label='Satélite')
    axes[0].set(title='RECM por péntada', xlabel='Péntada', ylabel='mm'); axes[0].legend(frameon=False)
    res = predicciones[referencia][0] - obs; m_res, s_res = np.nanmean(res, 0), np.nanstd(res, 0)
    axes[1].fill_between(t, m_res - s_res, m_res + s_res, color=COLORES[referencia], alpha=.2, label='± 1 desv.')
    axes[1].plot(t, m_res, c=COLORES[referencia], lw=2, label=f'Residuo medio {referencia}')
    axes[1].plot(t, np.nanmean(D['satelite_obs'] - obs, 0), c=COLORES['Satélite'], label='Residuo medio satélite')
    axes[1].axhline(0, c='k', lw=.5); axes[1].set(title='Residuo (predicho − observado)', xlabel='Péntada', ylabel='mm')
    axes[1].legend(frameon=False)
    return fig

def grafica_mapa_residuos(D, prediccion, titulo):
    res = prediccion - D['observado']; orden = np.argsort(D['altitud'][D['va']])
    fig, ax = plt.subplots(figsize=(12, 5))
    im = ax.imshow(res[orden], aspect='auto', cmap='RdBu_r', vmin=-60, vmax=60, interpolation='nearest')
    ax.set(xlabel='Péntada', ylabel='Curvas de validación (por altitud, de menor a mayor)', title=titulo)
    fig.colorbar(im, ax=ax, label='predicho − observado (mm)')
    return fig

# Error y sesgo por curva según altitud y fase ENSO (ventana DEF, criterio NOAA)
def diagnostico_altitud_enso(D, prediccion, nombre='KAN'):
    obs = D['observado']
    rmse = lambda p: np.sqrt(np.nanmean((p - obs)**2, 1))
    sesgo = lambda p: 100 * (np.nansum(np.where(np.isnan(obs), np.nan, p), 1) / np.nansum(obs, 1) - 1)
    diag = pd.DataFrame({'altitud': D['altitud'][D['va']], 'fase': D['fase_def'],
                         f'RECM_{nombre}': rmse(prediccion), 'RECM_satelite': rmse(D['satelite_obs']),
                         f'sesgo_{nombre}': sesgo(prediccion), 'sesgo_satelite': sesgo(D['satelite_obs'])})
    fig, axes = plt.subplots(1, 3, figsize=(15, 3.8))
    axes[0].plot(diag.altitud, diag.RECM_satelite - diag[f'RECM_{nombre}'], '.', ms=3, alpha=.6); axes[0].axhline(0, c='k', lw=.5)
    axes[0].set(xlabel='Altitud (m)', ylabel=f'RECM satélite − RECM {nombre} (mm)', title=f'Ganancia de la {nombre} según altitud')
    axes[1].plot(diag.altitud, diag[f'sesgo_{nombre}'], '.', ms=3, alpha=.6, c=COLORES.get(nombre), label=nombre)
    axes[1].plot(diag.altitud, diag.sesgo_satelite, '.', ms=3, alpha=.4, c=COLORES['Satélite'], label='Satélite')
    axes[1].axhline(0, c='k', lw=.5)
    axes[1].set(xlabel='Altitud (m)', ylabel='Sesgo del acumulado (%)', title='Sesgo según altitud', ylim=(-100, 150)); axes[1].legend(frameon=False)
    grupos = ['La Niña', 'Neutral', 'El Niño']
    axes[2].boxplot([diag.loc[diag.fase == g, f'sesgo_{nombre}'] for g in grupos], tick_labels=grupos, showfliers=False)
    axes[2].axhline(0, c='k', lw=.5); axes[2].set(ylabel=f'Sesgo del acumulado {nombre} (%)', title='Sesgo según fase ENSO (ventana DEF)')
    tabla = diag.groupby('fase')[[f'RECM_{nombre}', 'RECM_satelite', f'sesgo_{nombre}', 'sesgo_satelite']].median().round(1)
    return fig, tabla
