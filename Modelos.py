# Mariana: codificador funcional y cabeza lineal; FPCA, pérdida y entrenamiento se integran después.
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

class ModeloMariana(nn.Module):
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

# Revisión de configuraciones sin entrenamiento

def revisar_configuraciones(datos, W, R, cov_norm, normalizador, salida, neuronas, componentes, semillas):
    idx = np.flatnonzero(datos['particion'] == 'entrenamiento')[:16]
    c = torch.as_tensor(datos['c'][idx], dtype=torch.float64)
    cov = torch.as_tensor(cov_norm[idx], dtype=torch.float64)
    errores = []
    for q in neuronas:
        for m in componentes:
            for seed in semillas:
                establecer_semilla(seed)
                modelo = ModeloMariana(W, R, q, m)
                if modelo(c, cov).shape != (len(idx), m): raise AssertionError('Dimensiones de salida incorrectas.')
                parametros = [p.detach().clone() for p in modelo.parameters()]
                establecer_semilla(seed); replica = ModeloMariana(W, R, q, m)
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

# Ejecución
def main():
    parser = argparse.ArgumentParser(description='Codificador y cabeza lineal, parte de Mariana')
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
    revisar_configuraciones(datos, W, R, cov_norm, normalizador, Path(args.salida),
                            (args.neuronas,) if args.neuronas else NEURONAS,
                            (args.m_fpca,) if args.m_fpca else COMPONENTES, SEMILLAS)

if __name__ == '__main__': main()