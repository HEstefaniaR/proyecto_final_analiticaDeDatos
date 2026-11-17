import argparse
import random
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from torch import nn

NEURONAS = (6, 7, 8)
COMPONENTES = (4, 5, 6)
SEMILLAS = (1, 2, 3)  # 0001, 0002, 0003
# Parámetros ajustables: cambiar OCULTAS tras conocer los parámetros de la cabeza KAN
OCULTAS = 8
ACTIVACION = nn.Tanh
OCULTAS_KAN = 2
INTERVALOS_KAN = 12


# Datos exportados desde R
def cargar_datos(carpeta='Curvas'):
    carpeta = Path(carpeta)
    archivos = ('metadatos_curvas.csv', 'coeficientes_satelite.csv', 'base_satelite_malla_fina.csv', 'parametros_construccion.csv')
    faltantes = [f for f in archivos if not (carpeta / f).is_file()]
    if faltantes: raise FileNotFoundError(f'Faltan archivos en {carpeta}: {faltantes}')
    meta = pd.read_csv(carpeta / archivos[0], dtype={'curva': str, 'STATION_ID': str})
    coefs = pd.read_csv(carpeta / archivos[1], dtype={'curva': str})
    base = pd.read_csv(carpeta / archivos[2]); parametros = pd.read_csv(carpeta / archivos[3])
    if meta.curva.duplicated().any() or coefs.curva.duplicated().any() or set(meta.curva) != set(coefs.curva):
        raise ValueError('Identificadores de curvas duplicados o no coincidentes.')
    columnas = [c for c in coefs if c != 'curva']
    if columnas != [c for c in base if c != 't']: raise ValueError('Base y coeficientes satelitales incompatibles.')
    k_sat = parametros.loc[parametros.fuente.eq('satelite'), ['K', 'orden']]
    k_est = parametros.loc[parametros.fuente.eq('estaciones'), 'K']
    if len(k_sat) != 1 or len(k_est) != 1 or int(k_sat.iloc[0].K) != len(columnas) or int(k_sat.iloc[0].orden) != 4:
        raise ValueError('Parámetros de construcción incompatibles con la exportación R.')
    datos = meta.merge(coefs, on='curva', validate='one_to_one', sort=False)
    requeridas = ('altitud_completa', 'oni_promedio', 'particion', 'YEAR')
    if any(c not in datos for c in requeridas): raise ValueError(f'Faltan columnas: {requeridas}')
    anios = datos.YEAR.to_numpy(); particion = datos.particion.to_numpy()
    esperada = np.where(anios <= 1995, 'entrenamiento', np.where(anios <= 2000, 'validacion', 'prueba'))
    if np.any((anios < 1981) | (anios > 2003)) or np.any(particion != esperada):
        raise ValueError('Particiones incorrectas: 1981-95 / 1996-2000 / 2001-03.')
    c = datos[columnas].to_numpy(dtype=np.float64)
    cov = datos[['altitud_completa', 'oni_promedio']].to_numpy(dtype=np.float64)
    t = base.t.to_numpy(dtype=np.float64); phi = base[columnas].to_numpy(dtype=np.float64)
    if not all(np.isfinite(x).all() for x in (c, cov, t, phi)): raise ValueError('Hay entradas faltantes o no finitas.')
    if not np.all(np.diff(t) > 0) or abs(t[0]) > 1e-9 or abs(t[-1]-72) > 1e-9:
        raise ValueError('La malla de integración debe cubrir [0,72] en orden creciente.')
    return {'curva': datos.curva.to_numpy(), 'particion': particion, 'c': c, 'cov': cov, 't': t, 'phi': phi, 'k_est': int(k_est.iloc[0])}

# Gram W y rugosidad R por cuadratura
def matrices_integracion(t, phi):
    dt = np.diff(t); pesos = np.r_[dt[0]/2, (dt[:-1]+dt[1:])/2, dt[-1]/2]
    W = phi.T @ (pesos[:, None]*phi)
    segunda = np.gradient(np.gradient(phi, t, axis=0, edge_order=2), t, axis=0, edge_order=2)
    R = segunda.T @ (pesos[:, None]*segunda)
    return (W+W.T)/2, (R+R.T)/2

# Estandarización ajustada solo con entrenamiento
class NormalizadorCovariables:
    def ajustar(self, cov, particion):
        train = cov[np.asarray(particion) == 'entrenamiento']
        if len(train) == 0: raise ValueError('No hay curvas de entrenamiento.')
        self.media = train.mean(axis=0); self.desviacion = np.maximum(train.std(axis=0), 1e-12)
        return self
    def transformar(self, cov): return (cov-self.media)/self.desviacion

# Neuronas funcionales z = A W c
class CodificadorFuncional(nn.Module):
    def __init__(self, gram, neuronas=8, rugosidad=None):
        super().__init__()
        if neuronas not in NEURONAS: raise ValueError(f'Neuronas permitidas: {NEURONAS}')
        W = torch.as_tensor(gram, dtype=torch.float64)
        if W.ndim != 2 or W.shape[0] != W.shape[1]: raise ValueError('Gram debe ser cuadrada.')
        self.register_buffer('gram', W.clone()); self.n_bases = W.shape[0]; self.n_neuronas = neuronas
        R = torch.as_tensor(rugosidad, dtype=torch.float64) if rugosidad is not None else torch.zeros_like(W)
        if R.shape != W.shape: raise ValueError('Rugosidad incompatible con Gram.')
        self.register_buffer('rugosidad', R.clone())
        self.A = nn.Parameter(torch.empty(neuronas, self.n_bases, dtype=torch.float64))
        nn.init.xavier_uniform_(self.A)
    def forward(self, c, cov):
        if c.ndim != 2 or c.shape[1] != self.n_bases or cov.shape != (len(c), 2):
            raise ValueError('Se esperan coeficientes (batch,K) y covariables (batch,2).')
        return torch.cat(((c @ self.gram) @ self.A.T, cov), dim=1)
    def penalizacion_rugosidad(self): return ((self.A @ self.rugosidad)*self.A).sum()
    @torch.no_grad()
    def funciones_pesos(self, phi):
        p = torch.as_tensor(phi, dtype=self.A.dtype, device=self.A.device)
        return (p @ self.A.T).cpu().numpy()

# Cabeza lineal de q+2 entradas a M coordenadas FPCA
class CabezaLineal(nn.Module):
    def __init__(self, neuronas=8, m_fpca=5):
        super().__init__()
        if neuronas not in NEURONAS or m_fpca not in COMPONENTES: raise ValueError('Se requiere q=6..8 y M=4..6.')
        self.capa = nn.Linear(neuronas+2, m_fpca, dtype=torch.float64)
    def forward(self, h): return self.capa(h)

# Codificador funcional + cabeza lineal
class ModeloCodificadorCabezaLineal(nn.Module):
    def __init__(self, gram, rugosidad, neuronas=8, m_fpca=5):
        super().__init__()
        self.codificador = CodificadorFuncional(gram, neuronas, rugosidad)
        self.cabeza = CabezaLineal(neuronas, m_fpca)
    def forward(self, c, cov): return self.cabeza(self.codificador(c, cov))

# Semillas e identidad algebraica de la cabeza lineal
def establecer_semilla(semilla):
    random.seed(semilla); np.random.seed(semilla); torch.manual_seed(semilla)

def verificar_equivalencia_lineal(modelo, c, cov, media_fpca, psi_fpca, tolerancia=1e-4):
    with torch.no_grad():
        A, W = modelo.codificador.A, modelo.codificador.gram
        B, b = modelo.cabeza.capa.weight, modelo.cabeza.capa.bias
        mu = torch.as_tensor(media_fpca, dtype=A.dtype, device=A.device)
        psi = torch.as_tensor(psi_fpca, dtype=A.dtype, device=A.device)
        c = torch.as_tensor(c, dtype=A.dtype, device=A.device)
        cov = torch.as_tensor(cov, dtype=A.dtype, device=A.device)
        q = modelo.codificador.n_neuronas
        if psi.ndim != 2 or psi.shape[1] != B.shape[0] or mu.shape != (psi.shape[0],):
            raise ValueError('FPCA incompatible con la cantidad de salidas M.')
        pred_red = mu + modelo(c, cov) @ psi.T
        beta = W @ A.T @ B[:, :q].T @ psi.T
        pred_equiv = (mu+b @ psi.T) + c @ beta + cov @ B[:, q:].T @ psi.T
        error = (pred_red-pred_equiv).abs().max().item()
        if not np.isfinite(error) or error > tolerancia: raise AssertionError(f'Equivalencia incumplida: {error:.3e}')
        return error

# Revisión de configuraciones de la cabeza lineal sin entrenamiento
def revisar_configuraciones_cabeza_lineal(datos, W, R, cov_norm, normalizador, salida, neuronas, componentes, semillas):
    idx = np.flatnonzero(datos['particion'] == 'entrenamiento')[:16]
    c = torch.as_tensor(datos['c'][idx], dtype=torch.float64)
    cov = torch.as_tensor(cov_norm[idx], dtype=torch.float64)
    errores = []
    for q in neuronas:
        for m in componentes:
            for seed in semillas:
                establecer_semilla(seed)
                modelo = ModeloCodificadorCabezaLineal(W, R, q, m)
                if modelo(c, cov).shape != (len(idx), m): raise AssertionError('Dimensiones de salida incorrectas.')
                parametros = [p.detach().clone() for p in modelo.parameters()]
                establecer_semilla(seed); replica = ModeloCodificadorCabezaLineal(W, R, q, m)
                if not all(torch.equal(x, y) for x, y in zip(parametros, replica.parameters())):
                    raise AssertionError('Inicialización no reproducible.')
                generador = torch.Generator().manual_seed(seed)
                mu = torch.randn(72, generator=generador, dtype=torch.float64)
                psi = torch.randn(72, m, generator=generador, dtype=torch.float64)
                errores.append(verificar_equivalencia_lineal(modelo, c, cov, mu, psi))
    salida.mkdir(parents=True, exist_ok=True)
    np.savez(salida/'normalizador_covariables.npz', media=normalizador.media, desviacion=normalizador.desviacion)
    print(f'Configuraciones comprobadas: {len(errores)} ({len(neuronas)} q × {len(componentes)} M × {len(semillas)} semillas).')
    print(f'Máximo error de equivalencia sintética: {max(errores):.3e} (no es FPCA real).')
    print('Sin entrenamiento ni selección de q o M.')
    print('Único archivo generado:', salida/'normalizador_covariables.npz')

# Cabeza MLP: recibe q características funcionales, altitud y ONI; produce M coordenadas FPCA
class CabezaMLP(nn.Module):
    def __init__(self, neuronas=8, m_fpca=5, ocultas=OCULTAS):
        super().__init__()
        if neuronas not in NEURONAS or m_fpca not in COMPONENTES or ocultas < 1:
            raise ValueError('Se requiere q=6..8, M=4..6 y ocultas >= 1.')
        self.red = nn.Sequential(nn.Linear(neuronas+2, ocultas), ACTIVACION(), nn.Linear(ocultas, m_fpca)).double()
    def forward(self, h): return self.red(h)

# Codificador funcional + cabeza MLP; el decodificador FPCA se conecta después
class ModeloCodificadorCabezaMLP(nn.Module):
    def __init__(self, gram, rugosidad, neuronas=8, m_fpca=5, ocultas=OCULTAS):
        super().__init__()
        self.codificador = CodificadorFuncional(gram, neuronas, rugosidad)
        self.cabeza = CabezaMLP(neuronas, m_fpca, ocultas)
    def forward(self, c, cov): return self.cabeza(self.codificador(c, cov))

def contar_parametros(modulo): return sum(p.numel() for p in modulo.parameters() if p.requires_grad)


# Aristas KAN con SiLU y B-splines cúbicos

class CapaKAN(nn.Module):
    def __init__(self, entradas, salidas, intervalos=INTERVALOS_KAN):
        super().__init__()
        if entradas < 1 or salidas < 1 or intervalos < 1:
            raise ValueError('Entradas, salidas e intervalos deben ser positivos.')
        self.entradas = entradas
        self.salidas = salidas
        self.intervalos = intervalos
        self.grado = 3
        nudos = -1 + 2 * torch.arange(-3, intervalos + 4, dtype=torch.float64) / intervalos
        self.register_buffer('nudos', nudos)
        self.register_buffer('centro', torch.zeros(entradas, dtype=torch.float64))
        self.register_buffer('escala', torch.ones(entradas, dtype=torch.float64))
        self.peso_base = nn.Parameter(torch.empty(salidas, entradas, dtype=torch.float64))
        self.peso_spline = nn.Parameter(torch.ones(salidas, entradas, dtype=torch.float64))
        self.coeficientes = nn.Parameter(torch.empty(salidas, entradas, intervalos + 3, dtype=torch.float64))
        nn.init.xavier_uniform_(self.peso_base)
        nn.init.normal_(self.coeficientes, std=0.1 / np.sqrt(entradas))

    def bases(self, x):
        # Recurrencia de Cox-de Boor sobre nudos uniformes extendidos
        valores = x.unsqueeze(-1)
        nudos = self.nudos
        bases = ((valores >= nudos[:-1]) & (valores < nudos[1:])).to(x.dtype)
        for grado in range(1, self.grado + 1):
            izquierda = (valores - nudos[:-(grado + 1)]) / (nudos[grado:-1] - nudos[:-(grado + 1)])
            derecha = (nudos[grado + 1:] - valores) / (nudos[grado + 1:] - nudos[1:-grado])
            bases = izquierda * bases[..., :-1] + derecha * bases[..., 1:]
        return bases

    def forward(self, x):
        normalizada = (x - self.centro) / self.escala
        base = nn.functional.silu(normalizada) @ self.peso_base.T
        coeficientes = self.peso_spline.unsqueeze(-1) * self.coeficientes
        spline = torch.einsum('big,oig->bo', self.bases(normalizada), coeficientes)
        return base + spline

    @torch.no_grad()
    def calibrar(self, x):
        # Escala fija de cada entrada, calculada antes del ajuste con entrenamiento
        if x.ndim != 2 or x.shape[1] != self.entradas or len(x) == 0 or not torch.isfinite(x).all():
            raise ValueError('La calibración requiere entradas finitas (batch, entradas).')
        minimo = x.amin(dim=0)
        maximo = x.amax(dim=0)
        self.centro.copy_((minimo + maximo) / 2)
        amplitud = (maximo - minimo) / 2
        self.escala.copy_(torch.where(amplitud > 0, amplitud, torch.ones_like(amplitud)))

    @torch.no_grad()
    def arista(self, salida, entrada, x):
        normalizada = (x - self.centro[entrada]) / self.escala[entrada]
        base = self.peso_base[salida, entrada] * nn.functional.silu(normalizada)
        spline = self.bases(normalizada) @ self.coeficientes[salida, entrada]
        return base + self.peso_spline[salida, entrada] * spline


# Composición de dos capas KAN, de q+2 entradas a M coordenadas FPCA
class CabezaKAN(nn.Module):
    def __init__(self, neuronas=8, m_fpca=5, ocultas=OCULTAS_KAN, intervalos=INTERVALOS_KAN):
        super().__init__()
        if neuronas not in NEURONAS or m_fpca not in COMPONENTES:
            raise ValueError('Se requiere q=6..8 y M=4..6.')
        self.capas = nn.ModuleList([
            CapaKAN(neuronas + 2, ocultas, intervalos),
            CapaKAN(ocultas, m_fpca, intervalos),
        ])

    def forward(self, h):
        for capa in self.capas:
            h = capa(h)
        return h

    @torch.no_grad()
    def calibrar(self, h):
        for capa in self.capas:
            capa.calibrar(h)
            h = capa(h)


class ModeloCodificadorCabezaKAN(nn.Module):
    def __init__(self, gram, rugosidad, neuronas=8, m_fpca=5, ocultas=OCULTAS_KAN, intervalos=INTERVALOS_KAN):
        super().__init__()
        self.codificador = CodificadorFuncional(gram, neuronas, rugosidad)
        self.cabeza = CabezaKAN(neuronas, m_fpca, ocultas, intervalos)

    def forward(self, c, cov):
        return self.cabeza(self.codificador(c, cov))


@torch.no_grad()
def ajustar_simbolica(capa, salida, entrada, x):
    # Polinomio cúbico en la entrada normalizada de una arista
    if x.ndim != 1 or len(x) < 8 or not torch.isfinite(x).all():
        raise ValueError('El ajuste simbólico requiere una rejilla de al menos ocho puntos finitos.')
    u = (x - capa.centro[entrada]) / capa.escala[entrada]
    valores = capa.arista(salida, entrada, x)
    matriz = torch.stack((torch.ones_like(u), u, u.square(), u.pow(3)), dim=1)
    if torch.linalg.matrix_rank(matriz[::2]) < 4:
        raise ValueError('El ajuste cúbico requiere cuatro columnas linealmente independientes.')
    coeficientes = torch.linalg.lstsq(matriz[::2], valores[::2]).solution
    ajuste = matriz @ coeficientes
    variacion = (valores[1::2] - valores[1::2].mean()).square().sum()
    # R² queda indefinido para una activación constante
    r2 = float('nan')
    if variacion > 0:
        r2 = (1 - (valores[1::2] - ajuste[1::2]).square().sum() / variacion).item()
    return coeficientes, ajuste, r2


# Decodificador FPCA
# Curvas de estaciones evaluadas en las 72 pentadas y coeficientes en su base
def cargar_curvas_estaciones(carpeta='Curvas'):
    carpeta = Path(carpeta)
    meta = pd.read_csv(carpeta / 'metadatos_curvas.csv', dtype={'curva': str})
    curvas = pd.read_csv(carpeta / 'curvas_estaciones_72_pentadas.csv', dtype={'curva': str})
    coefs = pd.read_csv(carpeta / 'coeficientes_estaciones.csv', dtype={'curva': str})
    base = pd.read_csv(carpeta / 'base_estaciones_malla_fina.csv')
    datos = meta[['curva', 'particion']].merge(curvas, on='curva', validate='one_to_one')
    datos = datos.merge(coefs, on='curva', validate='one_to_one')
    pentadas = [c for c in curvas if c != 'curva']
    columnas = [c for c in coefs if c != 'curva']
    # Base de estaciones en las pentadas 1, ..., 72 (la malla fina va de 0 a 72 con paso 0,1)
    filas = [int(np.argmin(np.abs(base.t.to_numpy() - t))) for t in range(1, 73)]
    phi = base.loc[filas, columnas].to_numpy(dtype=np.float64)
    return {'curva': datos.curva.to_numpy(), 'particion': datos.particion.to_numpy(),
            'y': datos[pentadas].to_numpy(dtype=np.float64), 'c': datos[columnas].to_numpy(dtype=np.float64),
            'phi': phi}

# FPCA discretizada: media mu(t) -> SVD -> autofunciones phi_k(t) y varianza explicada
# Las curvas ya vienen suavizadas, por eso no se suavizan de nuevo
def fit_fpca(curvas, n_componentes):
    Y = np.asarray(curvas, dtype=np.float64)
    media = Y.mean(0)
    _, sing, Vt = np.linalg.svd(Y - media, full_matrices=False)
    Vt = Vt * np.sign(Vt.sum(1, keepdims=True) + 1e-12)   # fija el signo (la SVD es ambigua en +-)
    explicada = sing**2 / np.sum(sing**2)
    return {'media': media, 'base': Vt[:n_componentes], 'explicada': explicada}

# Decodificador funcional: puntajes s (B, M) -> curva y(t) = mu(t) + sum_m s_m psi_m(t)  (B, 72)
# mu y psi son buffers: viajan en el state_dict pero NO se entrenan
class DecodificadorFPCA(nn.Module):
    def __init__(self, n_puntos=72, m_fpca=5):
        super().__init__()
        if m_fpca not in COMPONENTES: raise ValueError(f'Componentes permitidas: {COMPONENTES}')
        self.register_buffer('media', torch.zeros(n_puntos, dtype=torch.float64))
        self.register_buffer('base', torch.zeros(m_fpca, n_puntos, dtype=torch.float64))
    def fijar_fpca(self, fpca):
        self.media.copy_(torch.as_tensor(fpca['media'])); self.base.copy_(torch.as_tensor(fpca['base']))
        return self
    def forward(self, puntajes): return self.media + puntajes @ self.base
    def proyectar(self, curvas): return (curvas - self.media) @ self.base.T

# Coeficientes de mu y psi en la base de estaciones (para la perdida (c_hat - c)' W (c_hat - c))
# Las curvas viven en el espacio de las K funciones base, asi que hay a lo sumo K componentes
def coeficientes_fpca(fpca, phi):
    mu_c = np.linalg.lstsq(phi, fpca['media'], rcond=None)[0]
    psi_c = np.linalg.lstsq(phi, fpca['base'].T, rcond=None)[0].T
    error = max(np.abs(phi @ mu_c - fpca['media']).max(), np.abs(psi_c @ phi.T - fpca['base']).max())
    return mu_c, psi_c, error

# Varianza explicada y piso de reconstruccion en validacion para cada M
def tabla_componentes(curvas, fpca):
    entrenamiento = curvas['y'][curvas['particion'] == 'entrenamiento']
    validacion = curvas['y'][curvas['particion'] == 'validacion']
    filas = []
    for m in range(1, fpca['base'].shape[0] + 1):
        psi = fpca['base'][:m]
        reconstruida = fpca['media'] + (validacion - fpca['media']) @ psi.T @ psi
        filas.append({'M': m, 'varianza_explicada_acumulada': fpca['explicada'][:m].sum(),
                      'piso_RMSE_validacion_mm': np.sqrt(np.mean((reconstruida - validacion)**2))})
    return pd.DataFrame(filas), len(entrenamiento)


# Ejecución: codificador funcional + cabeza lineal
def main_cabeza_lineal():
    parser = argparse.ArgumentParser(description='Codificador funcional y cabeza lineal')
    parser.add_argument('--datos', default='Curvas')
    parser.add_argument('--neuronas', type=int, choices=NEURONAS, help='Si se omite, prueba q=6,7,8.')
    parser.add_argument('--m-fpca', type=int, choices=COMPONENTES, help='Si se omite, prueba M=4,5,6.')
    parser.add_argument('--salida', default='Resultados_Modelos1')
    args = parser.parse_args()
    datos = cargar_datos(args.datos); W, R = matrices_integracion(datos['t'], datos['phi'])
    if np.linalg.eigvalsh(W).min() <= 0: raise AssertionError('Gram no es definida positiva.')
    normalizador = NormalizadorCovariables().ajustar(datos['cov'], datos['particion'])
    cov_norm = normalizador.transformar(datos['cov'])
    for nombre in ('entrenamiento', 'validacion', 'prueba'):
        print(f'{nombre}: {(datos["particion"] == nombre).sum()} curvas')
    print(f'Bases: satélite={W.shape[0]}, estaciones={datos["k_est"]}; Gram simétrica={np.allclose(W, W.T)}')
    print(f'Semillas: {", ".join(f"{s:04d}" for s in SEMILLAS)}')
    revisar_configuraciones_cabeza_lineal(datos, W, R, cov_norm, normalizador, Path(args.salida),
                                          (args.neuronas,) if args.neuronas else NEURONAS,
                                          (args.m_fpca,) if args.m_fpca else COMPONENTES, SEMILLAS)


# Ejecución: codificador funcional + cabeza MLP
def main_cabeza_mlp():
    parser = argparse.ArgumentParser(description='Cabeza MLP: pruebas de integración sin entrenamiento')
    parser.add_argument('--datos', default='Curvas')
    parser.add_argument('--neuronas', type=int, choices=NEURONAS)
    parser.add_argument('--m-fpca', type=int, choices=COMPONENTES)
    parser.add_argument('--ocultas', type=int, default=OCULTAS, help='Ajustar según los parámetros de la cabeza KAN')
    parser.add_argument('--parametros-kan', type=int, help='Parámetros de la CABEZA KAN para la misma q y M')
    args = parser.parse_args()
    if args.ocultas < 1: parser.error('--ocultas debe ser >= 1')
    if args.parametros_kan is not None and (args.neuronas is None or args.m_fpca is None or args.parametros_kan < 1):
        parser.error('Para comparar con KAN indica --neuronas, --m-fpca y --parametros-kan positivo.')
    datos = cargar_datos(args.datos); W, R = matrices_integracion(datos['t'], datos['phi'])
    normalizador = NormalizadorCovariables().ajustar(datos['cov'], datos['particion'])
    cov_norm = normalizador.transformar(datos['cov'])
    idx = np.flatnonzero(datos['particion'] == 'entrenamiento')[:16]
    c = torch.as_tensor(datos['c'][idx], dtype=torch.float64)
    cov = torch.as_tensor(cov_norm[idx], dtype=torch.float64)
    opciones_q = (args.neuronas,) if args.neuronas else NEURONAS
    opciones_m = (args.m_fpca,) if args.m_fpca else COMPONENTES
    for q in opciones_q:
        for m in opciones_m:
            for semilla in SEMILLAS:
                establecer_semilla(semilla); modelo = ModeloCodificadorCabezaMLP(W, R, q, m, args.ocultas)
                s = modelo(c, cov)
                if s.shape != (len(idx), m) or not torch.isfinite(s).all(): raise AssertionError('Salida MLP incorrecta.')
                s.sum().backward()
                if modelo.codificador.A.grad is None or modelo.cabeza.red[0].weight.grad is None:
                    raise AssertionError('No fluye el gradiente hacia codificador y cabeza.')
                parametros = [p.detach().clone() for p in modelo.parameters()]
                establecer_semilla(semilla); copia = ModeloCodificadorCabezaMLP(W, R, q, m, args.ocultas)
                if not all(torch.equal(a, b) for a, b in zip(parametros, copia.parameters())):
                    raise AssertionError('Inicialización no reproducible.')
            n = contar_parametros(modelo.cabeza)
            print(f'q={q}, M={m}, ocultas={args.ocultas}: cabeza MLP={n} parámetros; 3 semillas correctas')
            if args.parametros_kan is not None:
                diferencia = abs(n-args.parametros_kan)/args.parametros_kan*100
                print(f'KAN={args.parametros_kan}; diferencia={diferencia:.2f}% -> {"CUMPLE" if diferencia <= 10 else "NO CUMPLE"} (máximo 10%)')
    print('Comprobaciones finalizadas. No se entrenó ni se generaron archivos.')
    return datos, W, R, cov_norm, opciones_q, opciones_m


# Ejecución: decodificador FPCA ajustado solo con entrenamiento
def main_decodificador():
    parser = argparse.ArgumentParser(description='Decodificador FPCA de las curvas de estaciones')
    parser.add_argument('--datos', default='Curvas')
    parser.add_argument('--salida', default='Resultados_Modelos1')
    args = parser.parse_args()
    curvas = cargar_curvas_estaciones(args.datos)
    entrenamiento = curvas['y'][curvas['particion'] == 'entrenamiento']
    fpca = fit_fpca(entrenamiento, curvas['phi'].shape[1])
    tabla, n_entrenamiento = tabla_componentes(curvas, fpca)
    mu_c, psi_c, error_base = coeficientes_fpca(fpca, curvas['phi'])
    print(f'FPCA con {n_entrenamiento} curvas de entrenamiento (1981-1995)')
    print(tabla.round(4).to_string(index=False))
    print(f'Diferencia máxima al pasar mu y psi a la base de estaciones: {error_base:.2e}')

    # Verificacion: con M = todas las componentes, decodificar los puntajes reproduce la curva
    decodificador = DecodificadorFPCA(72, max(COMPONENTES)).fijar_fpca(fpca)
    y = torch.as_tensor(entrenamiento[:16])
    print(f'Error de proyectar y decodificar (M = {max(COMPONENTES)}): '
          f'{(decodificador(decodificador.proyectar(y)) - y).abs().max().item():.2e}')

    salida = Path(args.salida); salida.mkdir(parents=True, exist_ok=True)
    np.savez(salida / 'decodificador_fpca.npz', media=fpca['media'], base=fpca['base'],
             explicada=fpca['explicada'], media_coeficientes=mu_c, base_coeficientes=psi_c)
    Path('Tablas').mkdir(exist_ok=True)
    tabla.to_csv('Tablas/07_decodificador_fpca_componentes.csv', index=False)
    print('Archivos:', salida / 'decodificador_fpca.npz', 'y Tablas/07_decodificador_fpca_componentes.csv')


if __name__ == '__main__':
    datos, W, R, cov_norm, opciones_q, opciones_m = main_cabeza_mlp()

    # Resolución de la malla para el seno de la función sintética
    torch.set_num_threads(1)
    x_malla = torch.linspace(-1, 1, 1001, dtype=torch.float64)
    x_intermedia = (x_malla[:-1] + x_malla[1:]) / 2
    resolucion = []
    for intervalos in (4, 8, 12, 16):
        capa = CapaKAN(1, 1, intervalos)
        coeficientes = torch.linalg.lstsq(capa.bases(x_malla), torch.sin(2 * torch.pi * x_malla)).solution
        error = capa.bases(x_intermedia) @ coeficientes - torch.sin(2 * torch.pi * x_intermedia)
        resolucion.append({'intervalos': intervalos, 'RMSE_aproximacion': error.square().mean().sqrt().item()})
    tabla_resolucion = pd.DataFrame(resolucion)
    print(tabla_resolucion.to_string(index=False))
    # Doce intervalos reducen el error de aproximación del seno a cerca de 0,00115

    # Comprobación sintética con ajuste, validación y verificación independientes
    generador = torch.Generator().manual_seed(2026)
    x_sintetica = 4 * torch.rand(2200, 8, generator=generador, dtype=torch.float64) - 2
    funcion = torch.sin(torch.pi * x_sintetica[:, 0]) + x_sintetica[:, 1].square()
    escalas_salida = torch.tensor([[0.5, 0.75, 1.0, 1.25]], dtype=torch.float64)
    y_sintetica = funcion[:, None] * escalas_salida
    x_ajuste, x_validacion, x_verificacion = x_sintetica[:1400], x_sintetica[1400:1800], x_sintetica[1800:]
    y_ajuste, y_validacion, y_verificacion = y_sintetica[:1400], y_sintetica[1400:1800], y_sintetica[1800:]
    resultados_sinteticos = []
    for semilla in SEMILLAS:
        establecer_semilla(semilla)
        cabeza = CabezaKAN(neuronas=6, m_fpca=4)
        cabeza.calibrar(x_ajuste)
        optimizador = torch.optim.AdamW(cabeza.parameters(), lr=0.01, weight_decay=0)
        mejor_validacion = float('inf')
        sin_mejora = 0
        for epoca in range(1, 4001):
            cabeza.train()
            optimizador.zero_grad(set_to_none=True)
            perdida = (cabeza(x_ajuste) - y_ajuste).square().mean()
            perdida.backward()
            optimizador.step()
            cabeza.eval()
            with torch.no_grad():
                perdida_validacion = (cabeza(x_validacion) - y_validacion).square().mean().item()
            if perdida_validacion < mejor_validacion - 1e-10:
                mejor_validacion = perdida_validacion
                mejor_estado = {nombre: valor.detach().clone() for nombre, valor in cabeza.state_dict().items()}
                mejor_epoca = epoca
                sin_mejora = 0
            else:
                sin_mejora += 1
            if sin_mejora >= 300:
                break
        cabeza.load_state_dict(mejor_estado)
        with torch.no_grad():
            error_verificacion = cabeza(x_verificacion) - y_verificacion
            rmse_salida = error_verificacion.square().mean(dim=0).sqrt()
            rmse = error_verificacion.square().mean().sqrt().item()
        resultados_sinteticos.append({
            'semilla': semilla, 'mejor_epoca': mejor_epoca,
            'RMSE_validacion': np.sqrt(mejor_validacion), 'RMSE_verificacion': rmse,
            'RMSE_max_salida': rmse_salida.max().item(),
        })
        print(f'KAN sintética, semilla {semilla:04d}, RMSE de verificación = {rmse:.6f}')
        if semilla == SEMILLAS[0]:
            cabeza_sintetica = cabeza
    tabla_sintetica = pd.DataFrame(resultados_sinteticos)
    print(tabla_sintetica.to_string(index=False))
    print(tabla_sintetica[['RMSE_verificacion', 'RMSE_max_salida']].agg(['mean', 'std']).to_string())
    assert np.isfinite(tabla_sintetica['RMSE_max_salida']).all()
    assert tabla_sintetica['RMSE_max_salida'].le(0.01).all(), 'La KAN supera el RMSE sintético permitido.'
    tabla_sintetica.to_csv('Tablas/06_KAN_sintetica.csv', index=False)
    # Las tres semillas quedan por debajo de 0,01 en cada una de las cuatro salidas sintéticas

    # Activaciones de la primera capa en el rango del ajuste sintético
    import matplotlib.pyplot as plt
    primera_capa = cabeza_sintetica.capas[0]
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.2))
    for entrada, ax in enumerate(axes):
        rejilla = torch.linspace(x_ajuste[:, entrada].min(), x_ajuste[:, entrada].max(), 200, dtype=torch.float64)
        for salida in range(OCULTAS_KAN):
            valores = primera_capa.arista(salida, entrada, rejilla)
            ax.plot(rejilla.numpy(), (valores - valores.mean()).numpy(), label=f'Arista hacia h{salida + 1}')
        ax.set_title(f'x{entrada + 1}')
        ax.legend(fontsize=8)
    plt.suptitle('Primera capa KAN · función sintética · semilla 0001')
    plt.tight_layout()
    fig.savefig('Figuras/06_KAN_activaciones_sinteticas.png', dpi=150)
    plt.close(fig)
    # En la semilla 0001, las aristas hacia h2 recuperan el seno de x1 y la curvatura de x2

    # Ajuste simbólico de la arista de x2 hacia la primera neurona oculta
    rejilla = torch.linspace(x_ajuste[:, 1].min(), x_ajuste[:, 1].max(), 200, dtype=torch.float64)
    coef_simbolicos, _, r2_simbolico = ajustar_simbolica(primera_capa, 0, 1, rejilla)
    tabla_simbolica = pd.DataFrame([{
        'semilla': SEMILLAS[0], 'entrada': 1, 'salida': 0,
        'centro': primera_capa.centro[1].item(), 'escala': primera_capa.escala[1].item(),
        'a0': coef_simbolicos[0].item(), 'a1': coef_simbolicos[1].item(),
        'a2': coef_simbolicos[2].item(), 'a3': coef_simbolicos[3].item(),
        'R2_verificacion': r2_simbolico,
    }])
    print('Expresión de la arista sintética: a0 + a1*u + a2*u² + a3*u³, u = (x-centro)/escala')
    print(tabla_simbolica.to_string(index=False))
    # El polinomio aproxima la arista sintética con R² cercano a 0,908 en los puntos intercalados
    tabla_simbolica.to_csv('Tablas/06_KAN_simbolica_sintetica.csv', index=False)

    # Integración del codificador y la KAN con entradas de entrenamiento
    idx_entrenamiento = datos['particion'] == 'entrenamiento'
    c_entrenamiento = torch.as_tensor(datos['c'][idx_entrenamiento], dtype=torch.float64)
    cov_entrenamiento = torch.as_tensor(cov_norm[idx_entrenamiento], dtype=torch.float64)
    comprobaciones = []
    parametros_cabezas = []
    for q in opciones_q:
        for m in opciones_m:
            for semilla in SEMILLAS:
                establecer_semilla(semilla)
                modelo = ModeloCodificadorCabezaKAN(W, R, q, m)
                with torch.no_grad():
                    h_entrenamiento = modelo.codificador(c_entrenamiento, cov_entrenamiento)
                    modelo.cabeza.calibrar(h_entrenamiento)
                c_lote, cov_lote = c_entrenamiento[:16], cov_entrenamiento[:16]
                coordenadas = modelo(c_lote, cov_lote)
                coordenadas.square().mean().backward()
                gradientes = all(
                    p.grad is not None and torch.isfinite(p.grad).all().item() and p.grad.abs().sum().item() > 0
                    for p in modelo.parameters()
                )
                establecer_semilla(semilla)
                replica = ModeloCodificadorCabezaKAN(W, R, q, m)
                with torch.no_grad():
                    replica.cabeza.calibrar(replica.codificador(c_entrenamiento, cov_entrenamiento))
                    error_replica = (coordenadas - replica(c_lote, cov_lote)).abs().max().item()
                comprobaciones.append({
                    'q': q, 'M': m, 'semilla': semilla,
                    'dimensiones_correctas': coordenadas.shape == (len(c_lote), m),
                    'salida_finita': torch.isfinite(coordenadas).all().item(),
                    'gradientes_finitos_no_nulos': gradientes, 'error_replica': error_replica,
                })
            n_kan = contar_parametros(modelo.cabeza)
            ocultas_mlp = max(1, round((n_kan - m) / (q + m + 3)))
            n_mlp = contar_parametros(CabezaMLP(q, m, ocultas_mlp))
            parametros_cabezas.append({
                'q': q, 'M': m, 'ocultas_KAN': OCULTAS_KAN, 'intervalos': INTERVALOS_KAN,
                'grado': 3, 'parametros_KAN': n_kan, 'ocultas_MLP': ocultas_mlp,
                'parametros_MLP': n_mlp, 'diferencia_pct': 100 * abs(n_mlp - n_kan) / n_kan,
                'parametros_codificador': contar_parametros(modelo.codificador),
            })
    tabla_comprobaciones = pd.DataFrame(comprobaciones)
    print(tabla_comprobaciones.to_string(index=False))
    assert tabla_comprobaciones[['dimensiones_correctas', 'salida_finita', 'gradientes_finitos_no_nulos']].all().all()
    assert tabla_comprobaciones['error_replica'].eq(0).all()
    # El gradiente alcanza el codificador y las dos capas KAN con inicializaciones reproducibles
    tabla_comprobaciones.to_csv('Tablas/06_KAN_comprobaciones.csv', index=False)
    tabla_parametros = pd.DataFrame(parametros_cabezas)
    print(tabla_parametros.to_string(index=False))
    assert tabla_parametros['diferencia_pct'].le(10).all()
    tabla_parametros.to_csv('Tablas/06_KAN_parametros.csv', index=False)
    print('La KAN devuelve M coordenadas latentes. Entrenamiento con lluvia pendiente del FPCA y la pérdida funcional.')
