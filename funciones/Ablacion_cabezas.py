from pathlib import Path
import numpy as np
import pandas as pd
import torch
import matplotlib.pyplot as plt

import Entrenamiento_red as Mo
import Ablacion_arquitectura as Aa

# Ablación de las cabezas lineal, MLP y KAN con los modelos entrenados en Entrenamiento_red.py (objetivo directo)
# Mismo codificador, mismo decodificador FPCA, parámetros de KAN y MLP igualados y 3 semillas por cabeza
CABEZAS = ('lineal', 'MLP', 'KAN')
SATELITE = 'Satélite sin corregir'
SEMILLAS = Mo.SEMILLAS
CARPETA_MODELOS = Path('Resultados_Modelos1')
REMUESTREOS = 2000
SEMILLA_BOOTSTRAP = 2026
# Métricas del enunciado (más la RECM) y si un valor menor es mejor
METRICAS = {'L2_relativo': True, 'KGE': False, 'sesgo_abs_%': True, 'pico1_pentadas': True, 'pico2_pentadas': True, 'RECM_mm': True}
COLORES = {'lineal': 'tab:red', 'MLP': 'tab:blue', 'KAN': 'tab:green', SATELITE: 'darkorange'}


# Curvas predichas en mm (malla fina) de cada cabeza y semilla para las curvas de validación
def predicciones_validacion(D):
    datos, va = D['datos'], D['va']
    W_s, R_s = Mo.matrices_integracion(datos['t'], datos['phi_sat'])
    c_sat = torch.as_tensor(datos['c_sat'][va])
    pred = {}
    for cabeza in CABEZAS:
        for s in SEMILLAS:
            g = torch.load(CARPETA_MODELOS / f'red_{cabeza}_semilla{s:04d}.pt', weights_only=False)
            modelo = Mo.crear_modelo(cabeza, W_s, R_s, g['q'], g['M'])
            modelo.load_state_dict(g['estado']); modelo.eval()
            cov = torch.as_tensor((datos['cov'][va] - g['media_covariables']) / g['desviacion_covariables'])
            with torch.no_grad(): puntajes = modelo(c_sat, cov).numpy()
            c_hat = g['media_coeficientes'] + puntajes @ g['base_coeficientes']
            pred[(cabeza, s)] = Mo.reconstruir_mm(c_hat, datos['phi_est'], g['escala'])[0]
    pred[(SATELITE, 1)] = D['satelite_mm']
    return pred


# Métricas de cada curva: L² relativo, KGE, acumulados, error en la fecha de los picos y error cuadrático en las péntadas
def metricas_por_curva(D, curva_fina):
    t, ref_f, obs = D['t_fina'], D['estacion_mm'], D['observado']
    pent, ref = curva_fina[:, D['filas']], ref_f[:, D['filas']]
    m = ~np.isnan(obs)
    tabla = pd.DataFrame({'L2_relativo': np.sqrt(np.trapezoid((curva_fina - ref_f)**2, t, axis=1) / np.trapezoid(ref_f**2, t, axis=1))})
    r = np.array([np.corrcoef(a, b)[0, 1] for a, b in zip(pent, ref)])
    tabla['KGE'] = 1 - np.sqrt((r - 1)**2 + (pent.std(1) / ref.std(1) - 1)**2 + (pent.mean(1) / ref.mean(1) - 1)**2)
    tabla['acumulado_predicho'] = np.where(m, pent, 0).sum(1)
    tabla['acumulado_observado'] = np.nansum(obs, 1)
    for j, (a, b) in enumerate(Aa.TEMPORADAS, 1):
        v = (t >= a) & (t <= b)
        tabla[f'pico{j}_pentadas'] = np.abs(t[v][curva_fina[:, v].argmax(1)] - t[v][ref_f[:, v].argmax(1)])
    tabla['error_cuadratico'] = np.where(m, (pent - np.nan_to_num(obs))**2, 0).sum(1)
    tabla['n_pentadas'] = m.sum(1)
    return tabla


# Métricas de un conjunto de curvas (filas de 'tabla', con repeticiones si viene de un remuestreo)
def agregar(tabla):
    sesgo = 100 * (tabla.acumulado_predicho.sum() / tabla.acumulado_observado.sum() - 1)
    return {'L2_relativo': tabla.L2_relativo.mean(), 'KGE': tabla.KGE.median(), 'sesgo_%': sesgo, 'sesgo_abs_%': abs(sesgo),
            'pico1_pentadas': tabla.pico1_pentadas.mean(), 'pico2_pentadas': tabla.pico2_pentadas.mean(),
            'RECM_mm': np.sqrt(tabla.error_cuadratico.sum() / tabla.n_pentadas.sum())}


def tablas_por_curva(D, pred):
    return {clave: metricas_por_curva(D, curva) for clave, curva in pred.items()}


# Métricas de cada modelo en validación: media y desviación estándar de las semillas
def tabla_metricas(por_curva):
    filas = []
    for modelo in (*CABEZAS, SATELITE):
        valores = pd.DataFrame([agregar(t) for (m, s), t in por_curva.items() if m == modelo])
        fila = {'modelo': modelo}
        for c in ('L2_relativo', 'KGE', 'sesgo_%', 'pico1_pentadas', 'pico2_pentadas', 'RECM_mm'):
            fila[c] = valores[c].mean(); fila[f'{c}_desv'] = valores[c].std() if len(valores) > 1 else np.nan
        filas.append(fila)
    return pd.DataFrame(filas).set_index('modelo')


# Bootstrap por estaciones: se remuestrean estaciones con reemplazo y cada una aporta todas sus curvas
# En cada remuestreo, la métrica de una cabeza es el promedio de sus 3 semillas
def bootstrap(D, por_curva, remuestreos=REMUESTREOS, semilla=SEMILLA_BOOTSTRAP):
    estaciones = D['meta'].STATION_ID.to_numpy()[D['va']]
    unicas = np.unique(estaciones)
    filas_estacion = {e: np.flatnonzero(estaciones == e) for e in unicas}
    rng = np.random.default_rng(semilla)
    registros = []
    for b in range(remuestreos):
        elegidas = rng.choice(unicas, size=len(unicas), replace=True)
        idx = np.concatenate([filas_estacion[e] for e in elegidas])
        for modelo in (*CABEZAS, SATELITE):
            valores = pd.DataFrame([agregar(t.iloc[idx]) for (m, s), t in por_curva.items() if m == modelo]).mean()
            registros.append({'remuestreo': b, 'modelo': modelo, **valores.to_dict()})
    return pd.DataFrame(registros)


# Diferencia A - B de cada métrica con su IC del 95 % (percentiles 2,5 y 97,5 del bootstrap)
# A supera a B solo si el intervalo excluye el cero en la dirección favorable a A
def comparar(por_curva, remuestreos_bootstrap):
    observado = {m: pd.DataFrame([agregar(t) for (mm, s), t in por_curva.items() if mm == m]).mean() for m in (*CABEZAS, SATELITE)}
    ancho = remuestreos_bootstrap.set_index(['remuestreo', 'modelo'])
    pares = [('MLP', 'lineal'), ('KAN', 'lineal'), ('KAN', 'MLP')] + [(c, SATELITE) for c in CABEZAS]
    filas = []
    for a, b in pares:
        for metrica, menor_mejor in METRICAS.items():
            dif = ancho.xs(a, level='modelo')[metrica].to_numpy() - ancho.xs(b, level='modelo')[metrica].to_numpy()
            inf, sup = np.percentile(dif, [2.5, 97.5])
            if inf > 0 or sup < 0:
                a_mejor = (sup < 0) if menor_mejor else (inf > 0)
                conclusion = f'{a} supera a {b}' if a_mejor else f'{b} supera a {a}'
            else:
                conclusion = 'Sin diferencia'
            filas.append({'A': a, 'B': b, 'metrica': metrica, 'diferencia_A_menos_B': observado[a][metrica] - observado[b][metrica],
                          'IC95_inferior': inf, 'IC95_superior': sup, 'conclusion': conclusion})
    return pd.DataFrame(filas)


def grafica_intervalos(comparacion):
    metricas = list(METRICAS)
    fig, axes = plt.subplots(1, len(metricas), figsize=(4 * len(metricas), 3.6), sharey=True)
    pares = comparacion[['A', 'B']].drop_duplicates().apply(lambda f: f'{f.A} − {f.B}', axis=1).tolist()
    for ax, metrica in zip(axes, metricas):
        sub = comparacion[comparacion.metrica == metrica].reset_index(drop=True)
        for i, f in sub.iterrows():
            color = 'k' if f.conclusion == 'Sin diferencia' else 'tab:purple'
            ax.plot([f.IC95_inferior, f.IC95_superior], [i, i], c=color, lw=2)
            ax.plot(f.diferencia_A_menos_B, i, 'o', c=color)
        ax.axvline(0, c='0.5', lw=0.8, ls='--')
        ax.set_title(metrica + (' (menor es mejor)' if METRICAS[metrica] else ' (mayor es mejor)'), fontsize=9)
    axes[0].set_yticks(range(len(pares))); axes[0].set_yticklabels(pares); axes[0].invert_yaxis()
    fig.supxlabel('Diferencia A − B con IC bootstrap del 95 % (morado: el intervalo excluye el cero)')
    fig.tight_layout()
    return fig


def grafica_media(D, pred):
    t = np.arange(1, 73); m_obs = np.nanmean(D['observado'], 0)
    fig, axes = plt.subplots(1, 2, figsize=(14, 4))
    axes[0].plot(t, m_obs, 'k', lw=2, label='Observado')
    for modelo in (*CABEZAS, SATELITE):
        curvas = [c[:, D['filas']] for (m, s), c in pred.items() if m == modelo]
        media = np.mean(curvas, axis=0).mean(0)
        axes[0].plot(t, media, c=COLORES[modelo], label=modelo)
        axes[1].plot(t, 100 * (media / m_obs - 1), c=COLORES[modelo], label=modelo)
    axes[0].set(xlabel='Péntada', ylabel='mm', title='Media en validación'); axes[0].legend(frameon=False)
    axes[1].axhline(0, c='k', lw=.5); axes[1].set(xlabel='Péntada', ylabel='%', title='Sesgo por péntada'); axes[1].legend(frameon=False)
    return fig
