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
OCULTAS_KAN = 4
INTERVALOS_KAN = 12

# Entrenamiento de la red funcional
# Escala del modelo: 'raiz' (raiz de P) o 'log1p' (log(1 + P)), aplicada a la curva ya suavizada en mm;
# con None no se entrena
ESCALA = 'raiz'
# Componentes del decodificador FPCA (M = 5 explica el 96 % de la varianza)
COMPONENTES_RED = (5,)
# Malla de lambda para la penalización de rugosidad; se elige por validación cruzada
LAMBDAS_RUGOSIDAD = (1e-2, 1e-1, 1.0, 1e1, 1e2, 1e3, 1e4)
EPOCAS_MAXIMAS = 3000
PACIENCIA = 200
TASA_APRENDIZAJE = 0.01
# Grupos de la validación cruzada por estaciones (columna grupo_vc exportada desde R)
GRUPOS_VC = 5


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
    requeridas = ('altitud_completa', 'oni_promedio', 'particion', 'YEAR', 'grupo_vc')
    if any(c not in datos for c in requeridas): raise ValueError(f'Faltan columnas: {requeridas}')
    anios = datos.YEAR.to_numpy(); particion = datos.particion.to_numpy()
    esperada = np.where(anios <= 1997, 'entrenamiento', 'prueba')
    if np.any((anios < 1981) | (anios > 2003)) or np.any(particion != esperada):
        raise ValueError('Particiones incorrectas: 1981-1997 / 1998-2003.')
    if datos.groupby('STATION_ID').particion.nunique().max() > 1:
        raise ValueError('Hay estaciones en más de un conjunto: la partición debe ser por estaciones.')
    grupos = datos.loc[datos.particion == 'entrenamiento', ['STATION_ID', 'grupo_vc']]
    if set(grupos.grupo_vc) != set(range(1, GRUPOS_VC + 1)) or grupos.groupby('STATION_ID').grupo_vc.nunique().max() > 1:
        raise ValueError(f'Cada estación de entrenamiento debe estar en un solo grupo de 1 a {GRUPOS_VC}.')
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
    # particion puede ser la columna de partición o una máscara de las curvas de ajuste
    def ajustar(self, cov, particion):
        filas = np.asarray(particion)
        train = cov[filas if filas.dtype == bool else filas == 'entrenamiento']
        if len(train) == 0: raise ValueError('No hay curvas de entrenamiento.')
        self.media = train.mean(axis=0); self.desviacion = np.maximum(train.std(axis=0), 1e-12)
        return self
    def transformar(self, cov): return (cov-self.media)/self.desviacion

# Neuronas funcionales
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
    print(f'Configuraciones comprobadas: {len(errores)} ({len(neuronas)} q x {len(componentes)} M x {len(semillas)} semillas).')
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
    # R2 queda indefinido para una activación constante
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

# Decodificador funcional: de los puntajes a la curva
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

# Varianza explicada y piso de reconstruccion en entrenamiento para cada M
def tabla_componentes(curvas, fpca):
    entrenamiento = curvas['y'][curvas['particion'] == 'entrenamiento']
    filas = []
    for m in range(1, fpca['base'].shape[0] + 1):
        psi = fpca['base'][:m]
        reconstruida = fpca['media'] + (entrenamiento - fpca['media']) @ psi.T @ psi
        filas.append({'M': m, 'varianza_explicada_acumulada': fpca['explicada'][:m].sum(),
                      'piso_RMSE_entrenamiento_mm': np.sqrt(np.mean((reconstruida - entrenamiento)**2))})
    return pd.DataFrame(filas), len(entrenamiento)


# Datos para entrenar: coeficientes en la escala del modelo
# La curva MCP en mm se recorta en 0, se transforma y se escribe de nuevo en la misma base por mínimos cuadrados
TRANSFORMACIONES = {'raiz': (np.sqrt, np.square), 'log1p': (np.log1p, np.expm1)}

def coeficientes_en_escala(c_mm, phi_fina, escala):
    curvas = TRANSFORMACIONES[escala][0](np.maximum(c_mm @ phi_fina.T, 0))
    return np.linalg.lstsq(phi_fina, curvas.T, rcond=None)[0].T

def cargar_datos_red(carpeta='Curvas', escala='raiz'):
    carpeta = Path(carpeta)
    leer = lambda archivo: pd.read_csv(carpeta / archivo, dtype={'curva': str, 'STATION_ID': str})
    meta = leer('metadatos_curvas.csv')
    orden = meta.curva.to_numpy()
    coeficientes = lambda archivo: leer(archivo).set_index('curva').loc[orden].to_numpy(dtype=np.float64)
    base_s = pd.read_csv(carpeta / 'base_satelite_malla_fina.csv')
    base_e = pd.read_csv(carpeta / 'base_estaciones_malla_fina.csv')
    t = base_s.t.to_numpy(dtype=np.float64)
    filas = [int(np.argmin(np.abs(t - k))) for k in range(1, 73)]
    phi_sat = base_s.drop(columns='t').to_numpy(dtype=np.float64)
    phi_est = base_e.drop(columns='t').to_numpy(dtype=np.float64)
    return {'meta': meta, 'particion': meta.particion.to_numpy(), 'grupo_vc': meta.grupo_vc.to_numpy(),
            'cov': meta[['altitud_completa', 'oni_promedio']].to_numpy(dtype=np.float64),
            'c_sat': coeficientes_en_escala(coeficientes('coeficientes_satelite.csv'), phi_sat, escala),
            'c_est': coeficientes_en_escala(coeficientes('coeficientes_estaciones.csv'), phi_est, escala),
            'est_mm': np.maximum(coeficientes('coeficientes_estaciones.csv') @ phi_est.T, 0),
            't': t, 'filas_pentadas': filas, 'phi_sat': phi_sat, 'phi_est': phi_est}

# MLP con el número de neuronas ocultas que deja sus parámetros más cerca de los de la KAN
# Neuronas ocultas del MLP para igualar los parámetros de la KAN
def neuronas_ocultas_mlp(q, m):
    n_kan = contar_parametros(CabezaKAN(q, m))
    return max(1, round((n_kan - m) / (q + m + 3)))

def crear_modelo(cabeza, gram, rugosidad, q, m):
    if cabeza == 'lineal': return ModeloCodificadorCabezaLineal(gram, rugosidad, q, m)
    if cabeza == 'MLP': return ModeloCodificadorCabezaMLP(gram, rugosidad, q, m, neuronas_ocultas_mlp(q, m))
    if cabeza == 'KAN': return ModeloCodificadorCabezaKAN(gram, rugosidad, q, m)
    raise ValueError('Cabezas permitidas: lineal, MLP y KAN.')

# FPCA de las curvas de estaciones de ajuste en las 72 pentadas (por defecto, todo el entrenamiento)
# Los puntajes se escalan por su desviación estándar para que las M salidas de la cabeza tengan magnitud similar;
# el decodificador multiplica de nuevo por esa desviación, así que la curva reconstruida es la misma
def preparar_fpca(datos, m, filas=None):
    if filas is None: filas = datos['particion'] == 'entrenamiento'
    phi = datos['phi_est'][datos['filas_pentadas']]
    curvas = datos['c_est'] @ phi.T
    entrenamiento = curvas[filas]
    fpca = fit_fpca(entrenamiento, m)
    desviacion = ((entrenamiento - fpca['media']) @ fpca['base'].T).std(axis=0)
    fpca_escalada = dict(fpca, base=fpca['base'] * desviacion[:, None])
    mu_c, psi_c, error = coeficientes_fpca(fpca_escalada, phi)
    return fpca_escalada, mu_c, psi_c, error

# Pérdida L2 funcional más la penalización de rugosidad de las funciones de pesos
class PerdidaFuncional(nn.Module):
    def __init__(self, gram_estaciones, mu_c, psi_c, lam):
        super().__init__()
        self.register_buffer('gram', torch.as_tensor(gram_estaciones, dtype=torch.float64))
        self.register_buffer('mu_c', torch.as_tensor(mu_c, dtype=torch.float64))
        self.register_buffer('psi_c', torch.as_tensor(psi_c, dtype=torch.float64))
        self.lam = lam
    def coeficientes(self, puntajes): return self.mu_c + puntajes @ self.psi_c
    def error_l2(self, puntajes, c):
        diferencia = self.coeficientes(puntajes) - c
        return ((diferencia @ self.gram) * diferencia).sum(dim=1).mean()
    def forward(self, puntajes, c, codificador):
        return self.error_l2(puntajes, c) + self.lam * codificador.penalizacion_rugosidad()

# RECM (en las péntadas) y L2 relativo en mm frente a la curva suavizada de la estación
def metricas_mm(c_hat, filas, datos, escala):
    pred, _ = reconstruir_mm(c_hat, datos['phi_est'], escala)
    ref, t, p = datos['est_mm'][filas], datos['t'], datos['filas_pentadas']
    recm = np.sqrt(np.mean((pred[:, p] - ref[:, p])**2))
    l2 = np.mean(np.sqrt(np.trapezoid((pred - ref)**2, t, axis=1) / np.trapezoid(ref**2, t, axis=1)))
    return recm, l2

# Entrenamiento con lote completo y AdamW
# tr y va son máscaras de las curvas de ajuste y de validación (en la validación cruzada, el grupo que queda fuera)
# Con va: parada temprana con la pérdida L2 de validación y se guarda el estado con menor pérdida
# Sin va (ajuste final): se entrena el número fijo de épocas y se guarda el último estado
# Con con_metricas=True también se guardan por época la RECM y el L2 relativo de ajuste y validación
def entrenar_red(cabeza, q, m, lam, semilla, datos, gram_s, rugosidad_s, gram_e, fpca, mu_c, psi_c, cov_norm, tr, va=None,
                 epocas=EPOCAS_MAXIMAS, paciencia=PACIENCIA, con_metricas=False, escala=ESCALA):
    T = lambda x: torch.as_tensor(x, dtype=torch.float64)
    c_tr, h_tr, y_tr = T(datos['c_sat'][tr]), T(cov_norm[tr]), T(datos['c_est'][tr])
    if va is not None: c_va, h_va, y_va = T(datos['c_sat'][va]), T(cov_norm[va]), T(datos['c_est'][va])
    establecer_semilla(semilla)
    modelo = crear_modelo(cabeza, gram_s, rugosidad_s, q, m)
    if cabeza == 'KAN':
        with torch.no_grad(): modelo.cabeza.calibrar(modelo.codificador(c_tr, h_tr))
    perdida = PerdidaFuncional(gram_e, mu_c, psi_c, lam)
    optimizador = torch.optim.AdamW(modelo.parameters(), lr=TASA_APRENDIZAJE, weight_decay=0)
    mejor, mejor_epoca, mejor_estado, sin_mejora = np.inf, 0, None, 0
    historia = []
    for epoca in range(1, epocas + 1):
        modelo.train(); optimizador.zero_grad(set_to_none=True)
        valor = perdida(modelo(c_tr, h_tr), y_tr, modelo.codificador)
        valor.backward(); optimizador.step()
        modelo.eval()
        validacion = np.nan
        if va is not None:
            with torch.no_grad(): validacion = perdida.error_l2(modelo(c_va, h_va), y_va).item()
        fila = [epoca, valor.item(), validacion]
        if con_metricas:
            with torch.no_grad():
                salida_tr = modelo(c_tr, h_tr)
                fila.append(perdida.error_l2(salida_tr, y_tr).item())
                fila += metricas_mm(perdida.coeficientes(salida_tr).numpy(), tr, datos, escala)
                if va is not None: fila += metricas_mm(perdida.coeficientes(modelo(c_va, h_va)).numpy(), va, datos, escala)
                else: fila += [np.nan, np.nan]
        historia.append(tuple(fila))
        if va is None:
            mejor_epoca = epoca
            mejor_estado = {k: v.detach().clone() for k, v in modelo.state_dict().items()}
        elif validacion < mejor - 1e-10:
            mejor, mejor_epoca, sin_mejora = validacion, epoca, 0
            mejor_estado = {k: v.detach().clone() for k, v in modelo.state_dict().items()}
        else:
            sin_mejora += 1
            if sin_mejora >= paciencia: break
    modelo.load_state_dict(mejor_estado)
    return modelo, perdida, {'perdida_validacion': mejor if va is not None else np.nan, 'mejor_epoca': mejor_epoca,
                             'epocas': epoca,
                             'parametros_cabeza': contar_parametros(modelo.cabeza),
                             'parametros_codificador': contar_parametros(modelo.codificador)}, historia

# Curva predicha en la malla fina: en la escala del modelo se recorta en 0 y luego se vuelve a mm
def reconstruir_mm(c_escala, phi_fina, escala):
    curva = c_escala @ phi_fina.T
    recortados = 100 * np.mean(curva < 0)
    return TRANSFORMACIONES[escala][1](np.maximum(curva, 0)), recortados

# Equivalencia de la cabeza lineal con la regresión función a función (integral por cuadratura)
@torch.no_grad()
def verificar_regresion_funcional(modelo, datos, fpca, cov_norm, idx):
    q = modelo.codificador.n_neuronas
    A = modelo.codificador.A.numpy()
    B, b = modelo.cabeza.capa.weight.numpy(), modelo.cabeza.capa.bias.numpy()
    psi, mu = fpca['base'], fpca['media']
    beta = (datos['phi_sat'] @ A.T) @ B[:, :q].T @ psi                       # beta(s, t): malla fina x 72 pentadas
    X = datos['c_sat'][idx] @ datos['phi_sat'].T                              # X(s) en la malla fina
    integral = np.trapezoid(X[:, :, None] * beta[None, :, :], datos['t'], axis=1)
    regresion = mu + b @ psi + cov_norm[idx] @ B[:, q:].T @ psi + integral
    red = (mu + modelo(torch.as_tensor(datos['c_sat'][idx]), torch.as_tensor(cov_norm[idx])).numpy() @ psi)
    return np.abs(regresion - red).max(), beta


# Ejecución: codificador funcional + cabeza lineal
def main_cabeza_lineal(opciones=()):
    parser = argparse.ArgumentParser(description='Codificador funcional y cabeza lineal')
    parser.add_argument('--datos', default='Curvas')
    parser.add_argument('--neuronas', type=int, choices=NEURONAS, help='Si se omite, prueba q=6,7,8.')
    parser.add_argument('--m-fpca', type=int, choices=COMPONENTES, help='Si se omite, prueba M=4,5,6.')
    parser.add_argument('--salida', default='Resultados_Modelos1')
    args = parser.parse_args(list(opciones))
    datos = cargar_datos(args.datos); W, R = matrices_integracion(datos['t'], datos['phi'])
    if np.linalg.eigvalsh(W).min() <= 0: raise AssertionError('Gram no es definida positiva.')
    normalizador = NormalizadorCovariables().ajustar(datos['cov'], datos['particion'])
    cov_norm = normalizador.transformar(datos['cov'])
    for nombre in ('entrenamiento', 'prueba'):
        print(f'{nombre}: {(datos["particion"] == nombre).sum()} curvas')
    print(f'Bases: satélite={W.shape[0]}, estaciones={datos["k_est"]}; Gram simétrica={np.allclose(W, W.T)}')
    print(f'Semillas: {", ".join(f"{s:04d}" for s in SEMILLAS)}')
    revisar_configuraciones_cabeza_lineal(datos, W, R, cov_norm, normalizador, Path(args.salida),
                                          (args.neuronas,) if args.neuronas else NEURONAS,
                                          (args.m_fpca,) if args.m_fpca else COMPONENTES, SEMILLAS)


# Ejecución: codificador funcional + cabeza MLP
def main_cabeza_mlp(opciones=()):
    parser = argparse.ArgumentParser(description='Cabeza MLP: pruebas de integración sin entrenamiento')
    parser.add_argument('--datos', default='Curvas')
    parser.add_argument('--neuronas', type=int, choices=NEURONAS)
    parser.add_argument('--m-fpca', type=int, choices=COMPONENTES)
    parser.add_argument('--ocultas', type=int, default=OCULTAS, help='Ajustar según los parámetros de la cabeza KAN')
    parser.add_argument('--parametros-kan', type=int, help='Parámetros de la CABEZA KAN para la misma q y M')
    args = parser.parse_args(list(opciones))
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
def main_decodificador(opciones=()):
    parser = argparse.ArgumentParser(description='Decodificador FPCA de las curvas de estaciones')
    parser.add_argument('--datos', default='Curvas')
    parser.add_argument('--salida', default='Resultados_Modelos1')
    args = parser.parse_args(list(opciones))
    curvas = cargar_curvas_estaciones(args.datos)
    entrenamiento = curvas['y'][curvas['particion'] == 'entrenamiento']
    fpca = fit_fpca(entrenamiento, curvas['phi'].shape[1])
    tabla, n_entrenamiento = tabla_componentes(curvas, fpca)
    mu_c, psi_c, error_base = coeficientes_fpca(fpca, curvas['phi'])
    print(f'FPCA con {n_entrenamiento} curvas de entrenamiento (1981-1997)')
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


# Gráfico de sedimentación de la FPCA: proporción de varianza de cada componente y acumulada,
# con la línea del 95 % y el M usado por la red
def grafica_sedimentacion(explicada, m_elegido, archivo):
    import matplotlib.pyplot as plt
    k = len(explicada)
    componentes = np.arange(1, k + 1)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.bar(componentes, 100 * explicada, color='#9ab8d8', label='Varianza de cada componente')
    ax.plot(componentes, 100 * np.cumsum(explicada), 'o-', color='#1f4e79', label='Varianza acumulada')
    ax.axhline(95, color='#c0392b', ls='--', lw=1.2, label='95 %')
    ax.axvline(m_elegido, color='#555555', ls=':', lw=1.2, label=f'M = {m_elegido}')
    for x, y in zip(componentes, 100 * np.cumsum(explicada)):
        ax.annotate(f'{y:.1f} %', (x, y), textcoords='offset points', xytext=(0, -15), ha='center', fontsize=8)
    ax.set(xlabel='Componente', ylabel='Varianza explicada (%)', ylim=(0, 110), xticks=componentes,
           title='FPCA de las curvas de estaciones (entrenamiento)')
    ax.legend(loc='center right', fontsize=8)
    fig.tight_layout(); fig.savefig(archivo, dpi=150); plt.close(fig)


# Historial por época del ajuste final de cada cabeza (mismas semillas y épocas que los modelos guardados)
COLUMNAS_HISTORIAL = ['epoca', 'perdida_entrenamiento', 'L2_validacion', 'L2_entrenamiento',
                      'RECM_entrenamiento', 'L2_relativo_entrenamiento', 'RECM_validacion', 'L2_relativo_validacion']

def historial_entrenamiento(escala=ESCALA, carpeta='Curvas', semillas=SEMILLAS, prefijo='08_red'):
    datos = cargar_datos_red(carpeta, escala)
    gram_s, rugosidad_s = matrices_integracion(datos['t'], datos['phi_sat'])
    gram_e, _ = matrices_integracion(datos['t'], datos['phi_est'])
    tr = datos['particion'] == 'entrenamiento'
    cov_norm = NormalizadorCovariables().ajustar(datos['cov'], tr).transformar(datos['cov'])
    elegidas = pd.read_csv(f'Tablas/{prefijo}_configuracion_elegida.csv').set_index('cabeza')
    filas = []
    for cabeza in ('lineal', 'MLP', 'KAN'):
        q, m, lam = int(elegidas.loc[cabeza, 'q']), int(elegidas.loc[cabeza, 'M']), float(elegidas.loc[cabeza, 'lambda'])
        epocas = int(elegidas.loc[cabeza, 'epocas_finales'])
        fpca, mu_c, psi_c, _ = preparar_fpca(datos, m)
        for semilla in semillas:
            _, _, resumen, historia = entrenar_red(cabeza, q, m, lam, semilla, datos, gram_s, rugosidad_s, gram_e,
                                                   fpca, mu_c, psi_c, cov_norm, tr, None, epocas,
                                                   con_metricas=True, escala=escala)
            for h in historia:
                filas.append({'cabeza': cabeza, 'semilla': semilla, **dict(zip(COLUMNAS_HISTORIAL, h)),
                              'mejor_epoca': resumen['mejor_epoca']})
            print(f'{cabeza} semilla={semilla:04d}: {resumen["mejor_epoca"]} épocas', flush=True)
    tabla = pd.DataFrame(filas)
    tabla.to_csv(f'Tablas/{prefijo}_historial_entrenamiento.csv', index=False)
    return tabla

# Pérdida y métricas por época: una fila por cabeza, entrenamiento en azul y validación en naranja si la hay
# Eje y en escala log porque las primeras épocas tienen errores muy grandes
def grafica_historial(tabla, archivo):
    import matplotlib.pyplot as plt
    paneles = [('Pérdida L2 (escala raíz)', 'L2_entrenamiento', 'L2_validacion'),
               ('RECM (mm)', 'RECM_entrenamiento', 'RECM_validacion'),
               ('L2 relativo', 'L2_relativo_entrenamiento', 'L2_relativo_validacion')]
    cabezas = ('lineal', 'MLP', 'KAN')
    fig, ejes = plt.subplots(len(cabezas), len(paneles), figsize=(14, 10), sharex='row')
    for i, cabeza in enumerate(cabezas):
        sub = tabla[tabla.cabeza == cabeza]
        for j, (titulo, col_tr, col_va) in enumerate(paneles):
            ax = ejes[i, j]
            for k, (semilla, g) in enumerate(sub.groupby('semilla')):
                g = g[g.epoca <= g.mejor_epoca]
                ax.plot(g.epoca, g[col_tr], color='#1f77b4', lw=1, alpha=0.8, label='Entrenamiento' if k == 0 else None)
                if g[col_va].notna().any():
                    ax.plot(g.epoca, g[col_va], color='#ff7f0e', lw=1, alpha=0.8, label='Validación' if k == 0 else None)
            ax.set_yscale('log')
            ax.set_title(f'{cabeza}: {titulo}', fontsize=10)
            if i == len(cabezas) - 1: ax.set_xlabel('Época')
            if i == 0 and j == 0: ax.legend(fontsize=8)
    fig.suptitle('Ajuste final de la red con todo el entrenamiento (configuración elegida, 3 semillas)')
    fig.tight_layout(); fig.savefig(archivo, dpi=150); plt.close(fig)
    return archivo

# Ejecución: búsqueda de q y lambda por validación cruzada de 5 grupos de estaciones y ajuste final
# En cada grupo, el normalizador, la FPCA y la calibración de la KAN se ajustan solo con los otros cuatro
# El grupo que queda fuera sirve para la parada temprana y para medir la pérdida
# Las curvas de prueba no se usan
def main_entrenamiento(escala=ESCALA, carpeta='Curvas', neuronas=NEURONAS, componentes=COMPONENTES_RED,
                       lambdas=LAMBDAS_RUGOSIDAD, semillas=SEMILLAS, epocas=EPOCAS_MAXIMAS, n_folds=GRUPOS_VC,
                       prefijo='08_red'):
    if escala not in TRANSFORMACIONES: raise ValueError(f'Escalas permitidas: {tuple(TRANSFORMACIONES)}')
    datos = cargar_datos_red(carpeta, escala)
    gram_s, rugosidad_s = matrices_integracion(datos['t'], datos['phi_sat'])
    gram_e, _ = matrices_integracion(datos['t'], datos['phi_est'])
    entrenamiento = datos['particion'] == 'entrenamiento'
    grupos = datos['grupo_vc']
    if set(grupos[entrenamiento]) != set(range(1, n_folds + 1)):
        raise ValueError(f'Los grupos de validación cruzada deben ir de 1 a {n_folds}.')
    cabezas = ('lineal', 'MLP', 'KAN')
    Path('Tablas').mkdir(exist_ok=True)

    # Parámetros de KAN y MLP para cada q y M (sum(p.numel()) de los parámetros entrenables)
    parametros = []
    for q in neuronas:
        for m in componentes:
            n_kan, n_mlp = contar_parametros(CabezaKAN(q, m)), contar_parametros(CabezaMLP(q, m, neuronas_ocultas_mlp(q, m)))
            parametros.append({'q': q, 'M': m, 'parametros_lineal': contar_parametros(CabezaLineal(q, m)),
                               'parametros_KAN': n_kan, 'ocultas_MLP': neuronas_ocultas_mlp(q, m), 'parametros_MLP': n_mlp,
                               'diferencia_pct': 100 * abs(n_mlp - n_kan) / n_kan})
    tabla_parametros = pd.DataFrame(parametros)
    assert tabla_parametros['diferencia_pct'].le(10).all(), 'KAN y MLP difieren en más de 10 % de parámetros.'
    tabla_parametros.to_csv(f'Tablas/{prefijo}_parametros_cabezas.csv', index=False)

    # Varianza explicada acumulada de la FPCA para M = 1, ..., K con todo el entrenamiento (justificación de M)
    k_est = datos['phi_est'].shape[1]
    explicada = preparar_fpca(datos, k_est)[0]['explicada']
    pd.DataFrame({'M': range(1, k_est + 1), 'varianza_explicada_acumulada': np.cumsum(explicada[:k_est])}) \
        .to_csv(f'Tablas/{prefijo}_fpca_varianza.csv', index=False)
    Path('Figuras').mkdir(exist_ok=True)
    grafica_sedimentacion(explicada[:k_est], max(componentes), f'Figuras/{prefijo}_fpca_sedimentacion.png')

    # Búsqueda de hiperparámetros: las mismas q, M, lambda, semillas y grupos para las tres cabezas
    # Se guarda la predicción de cada grupo cuando queda fuera (coeficientes en la escala del modelo)
    busqueda, fuera_de_grupo = [], {}
    for grupo in range(1, n_folds + 1):
        tr = entrenamiento & (grupos != grupo)
        va = entrenamiento & (grupos == grupo)
        cov_norm = NormalizadorCovariables().ajustar(datos['cov'], tr).transformar(datos['cov'])
        fpcas = {m: preparar_fpca(datos, m, tr) for m in componentes}
        for cabeza in cabezas:
            for q in neuronas:
                for m in componentes:
                    fpca, mu_c, psi_c, _ = fpcas[m]
                    for lam in lambdas:
                        for semilla in semillas:
                            modelo, perdida, resumen, _ = entrenar_red(cabeza, q, m, lam, semilla, datos, gram_s,
                                                                       rugosidad_s, gram_e, fpca, mu_c, psi_c,
                                                                       cov_norm, tr, va, epocas)
                            clave = (cabeza, q, m, lam, semilla)
                            if clave not in fuera_de_grupo:
                                fuera_de_grupo[clave] = np.full(datos['c_est'].shape, np.nan)
                            with torch.no_grad():
                                fuera_de_grupo[clave][va] = perdida.coeficientes(
                                    modelo(torch.as_tensor(datos['c_sat'][va]), torch.as_tensor(cov_norm[va]))).numpy()
                            busqueda.append({'fold': grupo, 'cabeza': cabeza, 'q': q, 'M': m, 'lambda': lam,
                                             'semilla': semilla, 'curvas_ajuste': int(tr.sum()),
                                             'curvas_validacion': int(va.sum()), **resumen})
                            print(f'fold {grupo} {cabeza} q={q} M={m} lambda={lam} semilla={semilla:04d}: '
                                  f'validación={resumen["perdida_validacion"]:.4f} (época {resumen["mejor_epoca"]})',
                                  flush=True)
    tabla_busqueda = pd.DataFrame(busqueda)
    tabla_busqueda.to_csv(f'Tablas/{prefijo}_busqueda_hiperparametros.csv', index=False)
    resumen = (tabla_busqueda.groupby(['cabeza', 'q', 'M', 'lambda'])
               .agg(mean=('perdida_validacion', 'mean'), std=('perdida_validacion', 'std'),
                    epocas_finales=('mejor_epoca', 'median')).reset_index())
    resumen['epocas_finales'] = resumen['epocas_finales'].round().astype(int)
    # q y M comunes a las tres cabezas (mismo codificador y mismo decodificador): para cada (q, M) se toma
    # el mejor lambda de cada cabeza y se promedia la pérdida de validación cruzada de las tres cabezas
    # lambda se elige por cabeza porque es parte de la pérdida, no de la arquitectura
    mejor_lambda = resumen.loc[resumen.groupby(['cabeza', 'q', 'M'])['mean'].idxmin()]
    q_comun, m_comun = mejor_lambda.groupby(['q', 'M'])['mean'].mean().idxmin()
    elegidas = mejor_lambda[(mejor_lambda.q == q_comun) & (mejor_lambda.M == m_comun)].set_index('cabeza')
    elegidas.to_csv(f'Tablas/{prefijo}_configuracion_elegida.csv')
    print(elegidas.to_string())

    # Predicciones de validación cruzada de la configuración elegida (cada curva de entrenamiento predicha
    # por el modelo que no la vio); las usan la ablación de cabezas y la elección de la semilla en la evaluación
    carpeta_modelos = Path('Resultados_Modelos1'); carpeta_modelos.mkdir(exist_ok=True)
    predicciones_vc = {}
    for cabeza in cabezas:
        q, m, lam = int(elegidas.loc[cabeza, 'q']), int(elegidas.loc[cabeza, 'M']), float(elegidas.loc[cabeza, 'lambda'])
        for semilla in semillas:
            predicciones_vc[f'{cabeza}_semilla{semilla:04d}'] = fuera_de_grupo[(cabeza, q, m, lam, semilla)]
    np.savez(carpeta_modelos / 'predicciones_vc.npz', curva=datos['meta'].curva.to_numpy(), escala=escala,
             **predicciones_vc)

    # Ajuste final de cada cabeza con todo el entrenamiento: la mediana de las mejores épocas de los grupos
    normalizador = NormalizadorCovariables().ajustar(datos['cov'], entrenamiento)
    cov_norm = normalizador.transformar(datos['cov'])
    fpcas = {m: preparar_fpca(datos, m) for m in componentes}
    modelos = {}
    for cabeza in cabezas:
        q, m, lam = int(elegidas.loc[cabeza, 'q']), int(elegidas.loc[cabeza, 'M']), float(elegidas.loc[cabeza, 'lambda'])
        fpca, mu_c, psi_c, _ = fpcas[m]
        for semilla in semillas:
            modelos[(cabeza, semilla)] = entrenar_red(cabeza, q, m, lam, semilla, datos, gram_s, rugosidad_s, gram_e,
                                                      fpca, mu_c, psi_c, cov_norm, entrenamiento, None,
                                                      int(elegidas.loc[cabeza, 'epocas_finales']))[:2]

    # Equivalencia de la cabeza lineal elegida con la regresión función a función (cada semilla)
    # Se comprueba en 100 curvas de entrenamiento
    m = int(elegidas.loc['lineal', 'M'])
    idx_eq = np.flatnonzero(entrenamiento)[:100]
    equivalencia = []
    for semilla in semillas:
        modelo = modelos[('lineal', semilla)][0]
        error_integral, _ = verificar_regresion_funcional(modelo, datos, fpcas[m][0], cov_norm, idx_eq)
        error_coeficientes = verificar_equivalencia_lineal(modelo, datos['c_sat'][idx_eq], cov_norm[idx_eq],
                                                           fpcas[m][0]['media'], fpcas[m][0]['base'].T)
        equivalencia.append({'semilla': semilla, 'diferencia_max_integral': error_integral,
                             'diferencia_max_coeficientes': error_coeficientes})
    tabla_equivalencia = pd.DataFrame(equivalencia)
    assert tabla_equivalencia[['diferencia_max_integral', 'diferencia_max_coeficientes']].le(1e-4).all().all()
    tabla_equivalencia.to_csv(f'Tablas/{prefijo}_equivalencia_lineal.csv', index=False)
    print(tabla_equivalencia.to_string(index=False))

    # Valores recortados en 0 al volver a mm (entrenamiento, ajuste final de cada cabeza)
    recortes = []
    for cabeza in cabezas:
        for semilla in semillas:
            modelo, perdida = modelos[(cabeza, semilla)]
            with torch.no_grad():
                c_hat = perdida.coeficientes(modelo(torch.as_tensor(datos['c_sat'][entrenamiento]),
                                                    torch.as_tensor(cov_norm[entrenamiento]))).numpy()
            _, recortados = reconstruir_mm(c_hat, datos['phi_est'], escala)
            recortes.append({'cabeza': cabeza, 'semilla': semilla, 'valores_recortados_pct': recortados})
    pd.DataFrame(recortes).to_csv(f'Tablas/{prefijo}_valores_recortados.csv', index=False)

    # Pesos de los modelos finales (3 cabezas x semillas) con la FPCA y el normalizador de covariables
    for cabeza in cabezas:
        q, m, lam = int(elegidas.loc[cabeza, 'q']), int(elegidas.loc[cabeza, 'M']), float(elegidas.loc[cabeza, 'lambda'])
        fpca, mu_c, psi_c, _ = fpcas[m]
        for semilla in semillas:
            modelo = modelos[(cabeza, semilla)][0]
            torch.save({'estado': modelo.state_dict(), 'cabeza': cabeza, 'q': q, 'M': m, 'lambda': lam,
                        'semilla': semilla, 'escala': escala, 'epocas': int(elegidas.loc[cabeza, 'epocas_finales']),
                        'media_fpca': fpca['media'], 'base_fpca': fpca['base'],
                        'media_coeficientes': mu_c, 'base_coeficientes': psi_c,
                        'media_covariables': normalizador.media, 'desviacion_covariables': normalizador.desviacion},
                       carpeta_modelos / f'red_{cabeza}_semilla{semilla:04d}.pt')
    return {'datos': datos, 'modelos': modelos, 'elegidas': elegidas}
