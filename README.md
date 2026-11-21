# proyecto_final_analiticaDeDatos — Corrección funcional de la lluvia satelital en el Valle del Cauca

<img width="100" height="30" alt="uao-color" src="https://github.com/user-attachments/assets/b3112555-4eb4-4302-9a71-b88e845d8065" />

**Universidad Autónoma de Occidente**<br>
**Materia:** Analítica de Datos<br>
**Estudiantes:** Mariana Cifuentes, John David Franco, Estefanía Hernández, Natalia Paredes<br>
**Profesores:** Carlos Andrés Ferro Sánchez y Johann Alexis Ospina Galindez

---

## Descripción del repositorio

Cada año de lluvia en una localización del Valle del Cauca (1981–2003) se trata como un dato funcional sobre 72 péntadas. Hay dos fuentes:
- **Estaciones en tierra:** lluvia diaria, que se suma por péntada.
- **Satélite CHIRPS:** total por péntada en el píxel de cada estación.

El proyecto sigue cuatro fases:
1. Construcción de las curvas.
2. Análisis exploratorio funcional.
3. Comparación funcional entre fuentes.
4. Modelo que reconstruye la curva de la estación a partir de la curva satelital.

**El modelo es una red neuronal funcional** con tres bloques:
- **Codificador:** neuronas funcionales con pesos w_k(t) en B-splines, z = A W c.
- **Cabeza:** el bloque que se compara en la ablación. Hay tres:
  - **lineal:** equivale a una regresión funcional función a función;
  - **MLP:** red densa;
  - **KAN:** B-splines en las aristas.
- **Decodificador:** FPCA de las curvas de estaciones de entrenamiento, congelado.

Las tres cabezas usan el mismo codificador y el mismo decodificador. KAN y MLP tienen el mismo número de parámetros (±10 %), y cada configuración se entrena con 3 semillas. Todo se compara contra la curva satelital sin corregir y contra la regresión funcional función a función.

Los pesos finales están en el repositorio, así que no hace falta reentrenar para evaluar:

| Modelo | Pesos |
| --- | --- |
| Cabeza lineal (3 semillas) | `Resultados_Modelos1/red_lineal_semilla000{1,2,3}.pt` |
| Cabeza MLP (3 semillas) | `Resultados_Modelos1/red_MLP_semilla000{1,2,3}.pt` |
| Cabeza KAN (3 semillas) | `Resultados_Modelos1/red_KAN_semilla000{1,2,3}.pt` |

### Ruta rápida: solo evaluación

Las curvas están en `Curvas/` y los pesos en `Resultados_Modelos1/`.

1. Instalar el entorno (requiere `uv` y Python 3.10 a 3.13):

```bash
git clone https://github.com/HEstefaniaR/proyecto_final_analiticaDeDatos.git
cd proyecto_final_analiticaDeDatos
uv sync
```

2. Abrir `Evaluacion_interpretacion.ipynb` 

El notebook hace lo siguiente:
- elige la semilla de cada cabeza con la validación cruzada del notebook 01 (`Resultados_Modelos1/predicciones_vc.npz`);
- calcula en prueba las cuatro métricas con IC bootstrap por estaciones;
- separa los resultados por fase ENSO;
- interpreta la KAN;
- escribe `predicciones_estaciones_reservadas.csv`.

## Estructura del repositorio

```
proyecto_final_analiticaDeDatos/
├── Proyecto_Final_AnaliticaDeDatos.R   # preprocesamiento, partición, construcción, EDA, 
│                                       # comparación funcional, regresión funcional y curvas 
├── funciones/                          # Funciones de Python 
│   ├── Entrenamiento_red.py            
│   ├── Ablacion_arquitectura.py        
│   └── Ablacion_cabezas.py             
├── 01_Entrenamiento_red.ipynb          # Entrenamiento y selección con validación cruzada
├── 02_Ablacion_arquitectura.ipynb      
├── 03_Ablacion_cabezas.ipynb           
├── 04_Evaluacion_interpretacion.ipynb   
├── 05_Diagnostico_atenuacion_KAN.ipynb # Atenuación de la KAN y variantes con L1, entropía, poda y armónicos
├── data/                               # Datos crudos (estaciones, catálogos, CHIRPS, ONI, cartografía)
├── Curvas/                             # Curvas exportadas (coeficientes, bases, metadatos)
├── Resultados_Modelos1/                # Pesos de las tres cabezas × 3 semillas
├── Tablas/ 
├── Figuras/
├── predicciones_estaciones_reservadas.csv  # Predicciones de las estaciones de prueba (KAN)
├── visualizador.py                     # Visualizador en Streamlit de las estaciones de prueba
├── documentos/AI Usage.md              # Declaración de uso de IA e informe
├── pyproject.toml                      
└── README.md
```


Los nombres de las tablas y figuras empiezan por el número de la fase:

| Prefijo | Fase |
| --- | --- |
| `01` | Preprocesamiento y partición |
| `02` | ONI |
| `03` | Emparejamiento estación-píxel |
| `04` | Construcción |
| `05` | Exploratorio |
| `06` | Comparación y regresión, y las comprobaciones de la KAN |
| `08` | Red |
| `09` | Ablación de cabezas |
| `10` | Diagnóstico de atenuación de la KAN |


---

## Parte 1 - Entorno

**Python.** Requiere `uv` y Python entre 3.10 y 3.13. `requirements.txt` lista los mismos paquetes.

```bash
uv sync
```

```r
install.packages(c("fda", "readxl", "terra", "geodata", "tidyverse"))
```


## Parte 2 - Datos

Todo está en `data/`:

| Archivo | Contenido | Fuente |
| --- | --- | --- |
| `P_1980_1990.csv`, `P_1991_2003.csv` | Lluvia diaria de las estaciones | CVC |
| `CVC_Wheather_Stations_Long_Lat.csv` | Coordenadas de las estaciones | CVC |
| `CatalogoEstacionesElevacion.xlsx` | Altitud de las estaciones | CVC |
| `ImagenesSatelitePentadal/pptAAAAPP.tif` | Lluvia por péntada (año AAAA, péntada PP) | CHIRPS |
| `ONI_1981_2003.xlsx` | Índice ONI por ventana trimestral | NOAA CPC |
| `Paisajes_IGAC_hidrografia_IDEAM.gpkg` | Paisajes y relieve, y zonificación hidrográfica | IGAC e IDEAM |
| `oni_anual.csv`, `enso_evaluacion.csv` | ONI anual y fase ENSO por mes | Los genera el script de R |

## Parte 3 - Script de R
1. **Carga de datos** y días faltantes por año y estación.
2. **Partición por estaciones.**
   - **Sorteo por píxel CHIRPS:** las estaciones que comparten píxel ven la misma curva del satélite, así que se sortean juntas. Dentro de cada tercil de altitud, 15 % de los píxeles a prueba, elegidos entre los que tienen datos desde 1998; el resto queda en entrenamiento. Ningún píxel queda repartido entre entrenamiento y prueba.
   - **Años de cada grupo:** entrenamiento 1981–1997 y prueba 1998–2003, que incluye La Niña de 1998–2001 y El Niño de 2002–2003.
   - **Validación cruzada:** los píxeles de entrenamiento con datos hasta 1997 se reparten en 5 grupos dentro de cada tercil de altitud; cada píxel, con todas sus estaciones, queda en un solo grupo (`01_particion_estaciones.csv`, columna `grupo_vc`). La regresión funcional y la red usan los mismos grupos.
   - **Quedan fuera del modelo:** las estaciones sin coordenadas o sin altitud. Sí se usan en el EDA y en las comparaciones.
3. **Preprocesamiento:**
   - errores aislados de más de 400 mm;
   - totales mensuales registrados en un solo día;
   - imputación de días faltantes por IDW con las K estaciones de entrenamiento más cercanas que tienen dato ese día (peso 1/distancia, distancia euclidiana con las coordenadas), solo en estación-año con menos del 20 % de días faltantes; las vecinas son de la misma zona, valle (hasta 1239 m) o montaña, como en Russi y Larrahondo (2017). K se elige en cada zona por validación cruzada dejando una estación fuera, solo con las estaciones de entrenamiento y los días de 1981–1997 (`01_preprocesamiento_imputacion_seleccion_K.csv`), y ese K se aplica a todas las estaciones. La estación de la costa (planicie del Pacífico) no se imputa ni sirve de vecina (`01_preprocesamiento_imputacion.csv`);
   - cobertura mínima por estación y año (95 % de los días, huecos y bordes de máximo 7 días), después de imputar.

   Las reglas se revisan solo con estaciones de entrenamiento y se aplican a todos los datos. La trazabilidad queda en `01_preprocesamiento_trazabilidad.csv`.
4. **ONI y fases ENSO** (criterio NOAA).
5. **Emparejamiento estación-píxel** con las imágenes CHIRPS.
6. **Descriptivo de los datos emparejados**, antes de construir las curvas, con todos los pares limpios e imputados de 1981 a 2003.
   - Estadísticas por péntada de estación y satélite, para todos los años y por año: media, DE, coeficiente de variación, cuartiles, mínimo, máximo y porcentaje de péntadas y de días faltantes (`03_descriptivo_estadisticas.csv`, `03_descriptivo_CV_por_anio.png`).
   - Dispersión estación frente a satélite por año (`03_descriptivo_dispersion_por_anio.png`).
   - Correlación estación-satélite por péntada y de la lluvia media de cada curva con la altitud (`03_descriptivo_correlacion.csv`, `03_descriptivo_altitud.png`).
7. **Construcción de datos funcionales.**
   - Matriz Φ y mínimos cuadrados programados paso a paso.
   - Comparación de bases B-spline, Fourier y monomial, con MCO, MCP y penalización, usando ECM, VC y VCG.
   - B-spline cúbica con MCP y K elegido por VCG (K = 6 en estaciones y K = 11 en satélite).
   - Verificación de S, covarianza residual y bandas de confianza punto a punto.
8. **Análisis exploratorio funcional.**
   - Media verificada contra `fda`, varianza, covarianza, correlación y derivadas.
   - Curvas de cada año con su media y su desviación estándar, un panel por año con la fase ENSO según el ONI promedio (`05_exploratorio_curvas_por_anio_terrestre.png` y `05_exploratorio_curvas_por_anio_satelital.png`).
   - Profundidades de Fraiman–Muniz y MBD, y boxplot funcional con las curvas atípicas.
9. **Comparación funcional por fuente.**
   - Diferencia d(t) con su banda ± 2 EE y los tramos donde la banda no contiene el cero.
   - Correlación cruzada y prueba t.
   - ANOVA funcional por fase ENSO.
10. **Regresión funcional.** Modelos función a función (`linmod`) y concurrente con altitud y ONI (`fRegress`); λ se elige con validación cruzada de 5 grupos de estaciones y el modelo final se ajusta con todo el entrenamiento. La red se compara solo con el función a función.
11. **Curvas para Python.** Exporta `Curvas/` en la escala √P, solo con las curvas de entrenamiento y prueba, y el grupo de validación cruzada de cada curva de entrenamiento.

## Parte 4 - La red funcional
- **Codificador:** z_k = a_kᵀ W c, con W, la matriz de Gram de la base del satélite, calculada por cuadratura.
  - Neuronas funcionales: q ∈ {6, 7, 8}.
  - Al vector z se agregan la altitud y el ONI promedio del año.
  - Las funciones de pesos tienen penalización de rugosidad.
- **Cabezas:**
  - **lineal:** Σ_k β_k z_k;
  - **MLP:** una capa oculta con Tanh (53 neuronas con q = 8, para igualar los parámetros de la KAN);
  - **KAN:** una sola capa, sin capa oculta, de las q + 2 entradas a los M puntajes FPCA, con 12 intervalos y φ(x) = w_b·SiLU(x) + w_s·Σ c_i B_i(x) en cada arista. Es la arquitectura de menor pérdida de validación cruzada en la ablación del notebook 02.
- **Decodificador:**
  - **FPCA con M = 5** (96,5 % de la varianza): media y funciones propias de las curvas de estaciones de entrenamiento, congeladas.
- **Pérdida:** norma L² funcional (ĉ − c)ᵀ W (ĉ − c) más la penalización de rugosidad, en la escala √P. La predicción se recorta en 0 antes de volver a mm.
- **Entrenamiento:** lote completo y AdamW con tasa 0,01.
  - **Validación cruzada de 5 grupos de estaciones:** en cada grupo, el modelo se ajusta con los otros cuatro, hasta 3000 épocas, con parada temprana (paciencia de 200 épocas) sobre la pérdida del grupo que queda fuera. El normalizador de covariables, la FPCA y la calibración de la KAN se calculan solo con las curvas de ajuste.
  - **Selección:** q y λ con la pérdida L² media de los 5 grupos y las 3 semillas.
  - **Ajuste final:** cada cabeza se reajusta con todo el entrenamiento durante la mediana de las mejores épocas de los grupos. La prueba no se usa.

## Parte 5 - Notebooks

| Notebook | Qué hace |
| --- | --- |
| `01_Entrenamiento_red.ipynb` | Comprobaciones del codificador y de la MLP; KAN verificada con una función sintética (RMSE ≤ 0,01), activaciones y ajuste simbólico; búsqueda de q y λ por validación cruzada de 5 grupos de estaciones con 3 semillas y ajuste final con todo el entrenamiento; equivalencia numérica de la cabeza lineal con la regresión función a función; pérdida y métricas por época. Guarda los pesos de los modelos finales y las predicciones de validación cruzada de la configuración elegida en `Resultados_Modelos1/`. |
| `02_Ablacion_arquitectura.ipynb` | Ablación con validación cruzada de 5 grupos de estaciones de entrenamiento: q y M, capas e intervalos de la KAN y forma del ONI; diagnóstico del mejor modelo. No usa la prueba. |
| `03_Ablacion_cabezas.ipynb` | Lineal, MLP y KAN con las predicciones de validación cruzada del notebook 01: métricas, climatología, desfase de los picos y bootstrap por estaciones (2000 remuestreos).  |
| `Evaluacion_interpretacion.ipynb` | Evaluación final en prueba contra la curva suavizada de la estación: L² relativo, KGE, sesgo anual y error de la fecha de los dos picos, con IC bootstrap por estaciones; resultados por fase ENSO; w_k(t), activaciones de la KAN, contraste con la cabeza lineal, ajuste simbólico y peores casos; CSV de las estaciones de prueba. |
| `05_Diagnostico_atenuacion_KAN.ipynb` | Atenuación de la KAN con las predicciones de validación cruzada: techo de la FPCA, puntajes FPCA y picos frente al oráculo. Variantes de la KAN (L1, L1 con entropía, poda con reajuste y su control, armónicos del satélite) comparadas con la validación cruzada de 5 grupos. No usa la prueba. |


## Parte 6 - Predicciones de las estaciones de prueba

`predicciones_estaciones_reservadas.csv` tiene las predicciones de la KAN seleccionada para las estaciones de prueba (curvas estación-año de 1998–2003).

| Columna | Contenido |
| --- | --- |
| `station_id` | Código de la estación |
| `year` | Año |
| `day` | Día del año en que empieza la péntada (1, 6, …, 26, 32, …) |
| `value_mm` | Lluvia total predicha en la péntada (mm) |


**Visualizador.** `visualizador.py` muestra en un mapa del Valle del Cauca las 26 estaciones de prueba (1998–2003), coloreadas por la RECM de la KAN, con las estaciones de entrenamiento en gris como referencia. La barra de años, o un clic sobre la serie del ONI promedio, elige el año; al pasar el cursor por una estación se ve su curva (KAN, estación y satélite), y al pasar por la curva, la péntada con sus días. Usa los archivos de `Curvas/`, `predicciones_estaciones_reservadas.csv` y `data/oni_anual.csv`, y necesita internet para el mapa base (CARTO) y el relieve (OpenTopoMap):

```bash
uv run streamlit run visualizador.py
```

## Parte 7 - Resultados de referencia

**Partición** (`Tablas/01_particion_estaciones.csv`)

| Grupo | Estaciones sorteadas | Estaciones con curvas | Píxeles | Curvas estación-año |
| --- | --- | --- | --- | --- |
| Entrenamiento (1981–1997) | 160 | 143 | 120 | 1958 |
| Prueba (1998–2003) | 29 | 26 | 22 | 127 |

Las curvas de entrenamiento se reparten en 5 grupos de validación cruzada de 316 a 433 curvas.

**Configuración elegida con validación cruzada:** q = 8 y M = 5. λ es 10 en la lineal y en la KAN y 100 en la MLP; el ajuste final usa 2756, 700 y 91 épocas, y las cabezas KAN y MLP tienen 850 y 853 parámetros (`Tablas/08_red_configuracion_elegida.csv`, `Tablas/08_red_parametros_cabezas.csv`).

**Validación cruzada** (`Tablas/09_ablacion_cabezas_metricas.csv`, media de 3 semillas)

| Modelo | L² relativo | KGE | Sesgo del acumulado | RECM (mm) |
| --- | --- | --- | --- | --- |
| Satélite sin corregir | 0,575 | 0,386 | +8,5 % | 16,17 |
| Lineal | 0,416 | 0,533 | −6,9 % | 14,47 |
| MLP | 0,380 | 0,572 | −4,9 % | 12,28 |
| KAN | 0,368 | 0,597 | −4,2 % | 12,40 |

**Prueba** (`04_Evaluacion_interpretacion.ipynb`, semilla elegida con la validación cruzada)

| Modelo | L² relativo | KGE | Sesgo anual | Error de los dos picos (días) | RECM (mm) |
| --- | --- | --- | --- | --- | --- |
| Satélite sin corregir | 0,445 | 0,401 | +2,3 % | 22,95 | 13,94 |
| Lineal | 0,393 | 0,456 | −14,3 % | 14,39 | 13,08 |
| MLP | 0,356 | 0,510 | −9,5 % | 15,93 | 11,23 |
| KAN | 0,363 | 0,493 | −6,6 % | 14,67 | 10,95 |
