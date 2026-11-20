# proyecto_final_analiticaDeDatos — Corrección funcional de la lluvia satelital en el Valle del Cauca

<img width="100" height="30" alt="uao-color" src="https://github.com/user-attachments/assets/b3112555-4eb4-4302-9a71-b88e845d8065" />

**Universidad Autónoma de Occidente**<br>
**Materia:** Analítica de Datos<br>
**Estudiantes:** Mariana Cifuentes, John David Franco, Estefanía Hernández, Natalia Paredes<br>
**Profesor:** Carlos Andrés Ferro Sánchez y Johann Alexis Ospina Galindez

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
├── data/                               # Datos crudos (estaciones, catálogos, CHIRPS, ONI, cartografía)
├── Curvas/                             # Curvas exportadas (coeficientes, bases, metadatos)
├── Resultados_Modelos1/                # Pesos de las tres cabezas × 3 semillas
├── Tablas/ 
├── Figuras/
├── predicciones_estaciones_reservadas.csv  # Predicciones de las estaciones de prueba (KAN)
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
   - **Sorteo:** dentro de cada tercil de altitud, 15 % de las estaciones a prueba, elegidas entre las que tienen datos desde 1998; el resto queda en entrenamiento.
   - **Años de cada grupo:** entrenamiento 1981–1997 y prueba 1998–2003, que incluye La Niña de 1998–2001 y El Niño de 2002–2003.
   - **Validación cruzada:** las estaciones de entrenamiento con datos hasta 1997 se reparten en 5 grupos dentro de cada tercil de altitud; cada estación queda en un solo grupo (`01_particion_estaciones.csv`, columna `grupo_vc`). La regresión funcional y la red usan los mismos grupos.
   - **Quedan fuera del modelo:** las estaciones sin coordenadas o sin altitud. Sí se usan en el EDA y en las comparaciones.
3. **Preprocesamiento:**
   - errores aislados de más de 400 mm;
   - totales mensuales registrados en un solo día;
   - imputación de días faltantes por IDW con las K estaciones de entrenamiento más cercanas (K elegido en cada zona por validación cruzada dejando una estación fuera, `01_preprocesamiento_imputacion_seleccion_K.csv`) que tienen dato ese día (peso 1/distancia, distancia euclidiana con las coordenadas), solo en estación-año con menos del 20 % de días faltantes; las vecinas son de la misma zona, valle (hasta 1239 m) o montaña, como en Russi y Larrahondo (2017) (`01_preprocesamiento_imputacion.csv`);
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
  - **MLP:** una capa oculta con Tanh (63 neuronas con q = 8, para igualar los parámetros de la KAN);
  - **KAN:** una capa oculta de 4 nodos y 12 intervalos, con φ(x) = w_b·SiLU(x) + w_s·Σ c_i B_i(x).
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


## Parte 6 - Predicciones de las estaciones de prueba

`predicciones_estaciones_reservadas.csv` tiene las predicciones de la KAN seleccionada para las estaciones de prueba (curvas estación-año de 1998–2003).

| Columna | Contenido |
| --- | --- |
| `station_id` | Código de la estación |
| `year` | Año |
| `day` | Día del año en que empieza la péntada (1, 6, …, 26, 32, …) |
| `value_mm` | Lluvia total predicha en la péntada (mm) |

## Parte 7 - Resultados de referencia

**Partición**

| Grupo | Estaciones sorteadas | Curvas estación-año |
| --- | --- | --- |
| Entrenamiento (1981–1995) | 135 | 1307 |
| Validación (1996–2000) | 25 | 106 |
| Prueba (2001–2003) | 24 | 60 |


**Configuración elegida con validación:** q = 8 y M = 5. λ es 1 en la lineal, 10⁴ en la MLP y 10² en la KAN (`Tablas/08_red_configuracion_elegida.csv`).
