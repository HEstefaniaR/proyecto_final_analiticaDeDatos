####################################################
# Proyecto 3
# Por: Natalia Paredes, Mariana Cifuentes, John Franco
# y Estefanía Hernández
####################################################

##### Librerias #####
#install.packages("fda")
#install.packages("readxl")
#install.packages("tidyverse")
library(fda)
library(readxl)
library(terra)
library(geodata)
library(tidyverse)

set.seed(2026)

##### Directorio de trabajo #####
rm(list = ls())
dir.create("Figuras", showWarnings = FALSE)
dir.create("Resultados", showWarnings = FALSE)

##### Carga de datos #####
# Precipitacion diaria en tierra
p_tierra_1980_1990 <- read.csv("data/P_1980_1990.csv")
p_tierra_1991_2003 <- read.csv("data/P_1991_2003.csv")

p_tierra <- rbind(p_tierra_1980_1990, p_tierra_1991_2003)
p_tierra <- p_tierra[p_tierra$YEAR >= 1981, ]
p_tierra$fecha <- as.Date(sprintf("%04d-%02d-%02d", p_tierra$YEAR, p_tierra$MONTH, p_tierra$DAY))
p_tierra <- p_tierra[order(p_tierra$STATION_ID, p_tierra$fecha), ]

# Verificar duplicados estacion-fecha
sum(duplicated(p_tierra[, c("STATION_ID", "fecha")]))

# Separar datos
p_prueba <- p_tierra[p_tierra$YEAR >= 2001, ]
p_tierra <- p_tierra[p_tierra$YEAR <= 2000, ]

# Coordenadas y altitud de las estaciones
coord <- read.csv("data/CVC_Wheather_Stations_Long_Lat.csv", fileEncoding = "latin1")
coord$STATION_ID <- as.character(coord$CODIGO)

alt <- read_excel("data/CatalogoEstacionesElevacion.xlsx")
alt$STATION_ID <- as.character(alt$codigo)

info_estaciones <- merge(data.frame(STATION_ID = coord$STATION_ID, longitud = coord$Lon, latitud = coord$lat),
                         data.frame(STATION_ID = alt$STATION_ID, altitud = alt$alt),
                         by = "STATION_ID", all = TRUE)

# Borde del Valle
departamentos_col <- gadm(country = "COL", level = 1, path = tempdir())
valle <- departamentos_col[departamentos_col$NAME_1 == "Valle del Cauca", ]


##### Valores imposibles y errores de digitacion #####
# Valores negativos 
sum(p_tierra$VALUE < 0, na.rm = TRUE)
sum(p_prueba$VALUE < 0, na.rm = TRUE)

summary(p_tierra)
quantile(p_tierra$VALUE, probs = c(0.5, 0.9, 0.95, 0.99, 0.999, 0.99999), na.rm = TRUE)
tapply(p_tierra$VALUE, p_tierra$YEAR, summary)

# Valores mayores a 300 mm (1984: valor 859375)
mayores_300 <- p_tierra[p_tierra$VALUE > 300, c("STATION_ID", "YEAR", "MONTH", "DAY", "VALUE")]
mayores_300 <- mayores_300[order(-mayores_300$YEAR), ]
mayores_300

# Valores menores al percentil 99.9 (120 mm)
hist(p_tierra$VALUE[p_tierra$VALUE < 120],
     main = "Precipitación diaria", xlab = "Precipitación (mm)", ylab = "Frecuencia")

# Outliers en cada estacion y año
casos <- unique(mayores_300[, c("STATION_ID", "YEAR")])
casos <- casos[order(casos$YEAR, casos$STATION_ID), ]

par(mfrow = c(2, 3))
for (i in 1:nrow(casos)) {
  d <- p_tierra[p_tierra$STATION_ID == casos$STATION_ID[i] & p_tierra$YEAR == casos$YEAR[i], ]
  plot(d$fecha, d$VALUE, type = "b", pch = 16, cex = 0.4, lwd = 0.5, col = "steelblue",
       ylim = c(0, min(max(d$VALUE), 1000)),
       xlab = "Fecha", ylab = "Precipitación (mm)",
       main = paste("Estación", casos$STATION_ID[i], "-", casos$YEAR[i]))
  points(d$fecha[d$VALUE > 200], d$VALUE[d$VALUE > 200], col = "red", pch = 19)
}
par(mfrow = c(1, 1))

# Eliminar los valores aislados > 400 mm que son 5 veces mayores que sus vecinos
# y no caen el ultimo dia del mes
p_tierra$fin_mes <- format(p_tierra$fecha + 1, "%d") == "01"

n <- nrow(p_tierra)
anterior  <- c(NA, p_tierra$VALUE[-n])
siguiente <- c(p_tierra$VALUE[-1], NA)
vecino_max <- pmax(anterior, siguiente, na.rm = TRUE)

error <- which(p_tierra$VALUE > 400 & p_tierra$VALUE > 5 * vecino_max & !p_tierra$fin_mes)

data.frame(p_tierra[error, c("STATION_ID", "YEAR", "MONTH", "DAY")],
           anterior = anterior[error], VALUE = p_tierra$VALUE[error], siguiente = siguiente[error])

p_tierra$VALUE[error] <- NA


##### Acumulados mensuales #####
p_m <- p_tierra

# Resumen de registros por estacion y mes
resumen_mensual <- p_m %>%
  group_by(STATION_ID, YEAR, MONTH) %>%
  summarize(n_registros = n(),
            n_validos = sum(!is.na(VALUE)),
            n_ceros = sum(VALUE == 0, na.rm = TRUE),
            n_positivos = sum(VALUE > 0, na.rm = TRUE),
            total_observado = sum(VALUE, na.rm = TRUE),
            valor_fin_mes = if (any(fin_mes)) VALUE[fin_mes] else NA_real_,
            .groups = "drop")

resumen_mensual$proporcion_fin_mes <- resumen_mensual$valor_fin_mes / resumen_mensual$total_observado

# Candidatos: lluvia del ultimo dia >= 100 mm y >= 80 % del total del mes
candidatos_M <- resumen_mensual[which(resumen_mensual$valor_fin_mes >= 100 &
                                        resumen_mensual$proporcion_fin_mes >= 0.80), ]
nrow(candidatos_M)
print(candidatos_M, n = 40)

# Meses candidatos por estacion y año
print(count(candidatos_M, STATION_ID, YEAR, sort = TRUE), n = 30)

# Casos conocidos en detalle
print(resumen_mensual[resumen_mensual$STATION_ID == 5421200201 & resumen_mensual$YEAR == 1981, ], n = 12)
print(resumen_mensual[resumen_mensual$STATION_ID == 5300000201 & resumen_mensual$YEAR == 1986, ], n = 12)

# Meses con un unico registro, ubicado en el ultimo dia del mes
mensuales_M <- resumen_mensual[resumen_mensual$n_registros == 1 & resumen_mensual$n_validos == 1 &
                                 !is.na(resumen_mensual$valor_fin_mes), ]
nrow(mensuales_M)
print(mensuales_M, n = 40)

# El registro de fin de mes de los candidatos pasa a NA
clave_candidatos <- paste(candidatos_M$STATION_ID, candidatos_M$YEAR, candidatos_M$MONTH)
acumulados <- which(p_m$fin_mes & !is.na(p_m$VALUE) &
                      paste(p_m$STATION_ID, p_m$YEAR, p_m$MONTH) %in% clave_candidatos)
length(acumulados)

p_m_depurado <- p_m
p_m_depurado$VALUE[acumulados] <- NA


##### Ceros falsos #####
# Meses con muchos ceros y un pico al final del mes
diagnostico_ceros_M <- resumen_mensual[which(resumen_mensual$n_registros >= 25 & resumen_mensual$n_ceros >= 20 &
                                               resumen_mensual$n_positivos <= 3 & resumen_mensual$valor_fin_mes >= 30 &
                                               resumen_mensual$proporcion_fin_mes >= 0.80), ]
diagnostico_ceros_M

# Porcentaje de ceros de cada mes
resumen_mensual$pct_ceros <- 100 * resumen_mensual$n_ceros / resumen_mensual$n_validos

# Comparar con otras estaciones en el mismo mes y año, y con su historico 1981-1995
clave_ceros <- paste(diagnostico_ceros_M$STATION_ID, diagnostico_ceros_M$YEAR, diagnostico_ceros_M$MONTH)
referencia <- resumen_mensual[resumen_mensual$n_registros >= 25 & resumen_mensual$n_validos >= 25 &
                                !paste(resumen_mensual$STATION_ID, resumen_mensual$YEAR, resumen_mensual$MONTH) %in% clave_ceros, ]

diagnostico_ceros_M$pct_ceros <- round(100 * diagnostico_ceros_M$n_ceros / diagnostico_ceros_M$n_registros, 1)
diagnostico_ceros_M$mediana_mismo_anio <- NA
diagnostico_ceros_M$mediana_historica <- NA

for (i in 1:nrow(diagnostico_ceros_M)) {
  otras <- referencia[referencia$MONTH == diagnostico_ceros_M$MONTH[i] &
                        referencia$STATION_ID != diagnostico_ceros_M$STATION_ID[i], ]
  diagnostico_ceros_M$mediana_mismo_anio[i] <- median(otras$pct_ceros[otras$YEAR == diagnostico_ceros_M$YEAR[i]])
  diagnostico_ceros_M$mediana_historica[i] <- median(otras$pct_ceros[otras$YEAR <= 1995 &
                                                                       otras$YEAR != diagnostico_ceros_M$YEAR[i]])
}

diagnostico_ceros_M$diferencia_mismo_anio <- round(diagnostico_ceros_M$pct_ceros - diagnostico_ceros_M$mediana_mismo_anio, 1)
diagnostico_ceros_M$diferencia_historica <- round(diagnostico_ceros_M$pct_ceros - diagnostico_ceros_M$mediana_historica, 1)
print(diagnostico_ceros_M[order(-diagnostico_ceros_M$diferencia_mismo_anio), ], width = Inf)

# Meses completos con todos los valores en cero
resumen_mensual$dias_mes <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)[resumen_mensual$MONTH]
resumen_mensual$dias_mes[resumen_mensual$MONTH == 2 & resumen_mensual$YEAR %% 4 == 0] <- 29

todo_cero <- resumen_mensual[resumen_mensual$n_registros == resumen_mensual$dias_mes &
                               resumen_mensual$n_validos == resumen_mensual$dias_mes &
                               resumen_mensual$n_ceros == resumen_mensual$dias_mes, ]
nrow(todo_cero)

# Comparar con el historico de la estacion (1981-1995) y con la region el mismo mes
todo_cero$mediana_historica <- NA
todo_cero$mediana_regional <- NA

for (i in 1:nrow(todo_cero)) {
  historico <- resumen_mensual[resumen_mensual$STATION_ID == todo_cero$STATION_ID[i] &
                                 resumen_mensual$MONTH == todo_cero$MONTH[i] &
                                 resumen_mensual$YEAR <= 1995 & resumen_mensual$YEAR != todo_cero$YEAR[i] &
                                 resumen_mensual$n_validos >= 25, ]
  region <- resumen_mensual[resumen_mensual$YEAR == todo_cero$YEAR[i] & resumen_mensual$MONTH == todo_cero$MONTH[i] &
                              resumen_mensual$n_validos >= 25, ]
  todo_cero$mediana_historica[i] <- median(historico$pct_ceros)
  todo_cero$mediana_regional[i] <- median(region$pct_ceros)
}

todo_cero$diferencia_historica <- round(100 - todo_cero$mediana_historica, 1)
todo_cero$diferencia_regional <- round(100 - todo_cero$mediana_regional, 1)
print(todo_cero[order(-todo_cero$diferencia_regional), c("STATION_ID", "YEAR", "MONTH", "mediana_historica",
                                                         "mediana_regional", "diferencia_historica",
                                                         "diferencia_regional")], n = 30)

# Meses seleccionados: todos sus valores pasan a NA
meses_ceros <- data.frame(
  STATION_ID = c(2623300101, 2611400103, 2612710106, 2622900301, 2610000105,
                 2612000113, 2612000202, 2612700203, 2612720101, 2622600201),
  YEAR = c(1987, 1987, 1987, 1992, 1986, 1997, 1990, 1990, 1992, 1994),
  MONTH = c(11, 9, 7, 8, 12, 1, 2, 9, 10, 9))

ceros_falsos <- which(!is.na(p_m_depurado$VALUE) &
                        paste(p_m_depurado$STATION_ID, p_m_depurado$YEAR, p_m_depurado$MONTH) %in%
                        paste(meses_ceros$STATION_ID, meses_ceros$YEAR, meses_ceros$MONTH))
length(ceros_falsos)

p_m_depurado$VALUE[ceros_falsos] <- NA


##### Trazabilidad #####
trazabilidad <- rbind(data.frame(p_m[acumulados, c("STATION_ID", "fecha", "VALUE")],
                                 motivo = "Precipitacion final >= 100 mm y >= 80% del total mensual"),
                      data.frame(p_m[ceros_falsos, c("STATION_ID", "fecha", "VALUE")],
                                 motivo = "Mes con ceros sospechosos"))
trazabilidad

comparacion_M <- data.frame(
  etapa = c("Antes de depurar acumulados mensuales y ceros sospechosos",
            "Despues de depurar acumulados mensuales y ceros sospechosos"),
  registros = c(nrow(p_m), nrow(p_m_depurado)),
  validos = c(sum(!is.na(p_m$VALUE)), sum(!is.na(p_m_depurado$VALUE))),
  estaciones = c(length(unique(p_m$STATION_ID[!is.na(p_m$VALUE)])),
                 length(unique(p_m_depurado$STATION_ID[!is.na(p_m_depurado$VALUE)]))),
  pares_estacion_anio = c(nrow(unique(p_m[!is.na(p_m$VALUE), c("STATION_ID", "YEAR")])),
                          nrow(unique(p_m_depurado[!is.na(p_m_depurado$VALUE), c("STATION_ID", "YEAR")]))))
comparacion_M


##### Dias faltantes y años incompletos #####
# Racha mas larga de dias faltantes (como en el script de ANOVA funcional)
racha_max <- function(x) {
  r <- rle(is.na(x))
  if (any(r$values)) max(r$lengths[r$values]) else 0
}

# Dias sin dato al inicio o al final del año
borde_max <- function(x) {
  obs <- which(!is.na(x))
  if (length(obs) == 0) return(length(x))
  max(obs[1] - 1, length(x) - obs[length(obs)])
}

# Calendario completo de cada estacion-año
fechas <- data.frame(fecha = seq(as.Date("1981-01-01"), as.Date("2000-12-31"), by = "day"))
fechas$YEAR <- as.numeric(format(fechas$fecha, "%Y"))

calendario <- merge(unique(p_m_depurado[, c("STATION_ID", "YEAR")]), fechas, by = "YEAR")
calendario <- merge(calendario, p_m_depurado[, c("STATION_ID", "fecha", "VALUE")],
                    by = c("STATION_ID", "fecha"), all.x = TRUE)
calendario <- calendario[order(calendario$STATION_ID, calendario$fecha), ]

# Cobertura, hueco mas largo y dias faltantes en los bordes
cobertura <- calendario %>%
  group_by(STATION_ID, YEAR) %>%
  summarize(dias = n(),
            dias_validos = sum(!is.na(VALUE)),
            hueco = racha_max(VALUE),
            borde = borde_max(VALUE),
            .groups = "drop")

cobertura$cobertura <- cobertura$dias_validos / cobertura$dias
summary(cobertura)

# Comparar umbrales de cobertura en entrenamiento (1981-1995)
train <- cobertura[cobertura$YEAR <= 1995, ]
for (u in c(0.85, 0.90, 0.95)) {
  cat("Umbral", u, "- solo cobertura:", sum(train$cobertura >= u),
      "- con hueco y borde <= 7:", sum(train$cobertura >= u & train$hueco <= 7 & train$borde <= 7), "\n")
}

# Criterio: cobertura >= 95 %, hueco mas largo y borde de maximo 7 dias
cobertura$apto <- cobertura$cobertura >= 0.95 & cobertura$hueco <= 7 & cobertura$borde <= 7

cat("Pares estacion-año antes:", nrow(cobertura), "| aptos:", sum(cobertura$apto), "\n")
cat("Estaciones antes:", length(unique(cobertura$STATION_ID)),
    "| despues:", length(unique(cobertura$STATION_ID[cobertura$apto])), "\n")


##### Conjunto de prueba #####
# Se aplican las mismas reglas con los mismos umbrales.
# Los ceros falsos no se aplican porque esos meses se eligieron caso por caso.

# Errores aislados
p_prueba$fin_mes <- format(p_prueba$fecha + 1, "%d") == "01"
n <- nrow(p_prueba)
anterior  <- c(NA, p_prueba$VALUE[-n])
siguiente <- c(p_prueba$VALUE[-1], NA)
vecino_max <- pmax(anterior, siguiente, na.rm = TRUE)

error_prueba <- which(p_prueba$VALUE > 400 & p_prueba$VALUE > 5 * vecino_max & !p_prueba$fin_mes)
p_prueba[error_prueba, c("STATION_ID", "YEAR", "MONTH", "DAY", "VALUE")]
p_prueba$VALUE[error_prueba] <- NA

# Acumulados mensuales
resumen_prueba <- p_prueba %>%
  group_by(STATION_ID, YEAR, MONTH) %>%
  summarize(total_observado = sum(VALUE, na.rm = TRUE),
            valor_fin_mes = if (any(fin_mes)) VALUE[fin_mes] else NA_real_,
            .groups = "drop")

candidatos_prueba <- resumen_prueba[which(resumen_prueba$valor_fin_mes >= 100 &
                                            resumen_prueba$valor_fin_mes / resumen_prueba$total_observado >= 0.80), ]
acumulados_prueba <- which(p_prueba$fin_mes & !is.na(p_prueba$VALUE) &
                             paste(p_prueba$STATION_ID, p_prueba$YEAR, p_prueba$MONTH) %in%
                             paste(candidatos_prueba$STATION_ID, candidatos_prueba$YEAR, candidatos_prueba$MONTH))
length(acumulados_prueba)
p_prueba$VALUE[acumulados_prueba] <- NA

# Cobertura
fechas_prueba <- data.frame(fecha = seq(as.Date("2001-01-01"), as.Date("2003-12-31"), by = "day"))
fechas_prueba$YEAR <- as.numeric(format(fechas_prueba$fecha, "%Y"))

calendario_prueba <- merge(unique(p_prueba[, c("STATION_ID", "YEAR")]), fechas_prueba, by = "YEAR")
calendario_prueba <- merge(calendario_prueba, p_prueba[, c("STATION_ID", "fecha", "VALUE")],
                           by = c("STATION_ID", "fecha"), all.x = TRUE)
calendario_prueba <- calendario_prueba[order(calendario_prueba$STATION_ID, calendario_prueba$fecha), ]

cobertura_prueba <- calendario_prueba %>%
  group_by(STATION_ID, YEAR) %>%
  summarize(cobertura = mean(!is.na(VALUE)),
            hueco = racha_max(VALUE),
            borde = borde_max(VALUE),
            .groups = "drop")

cobertura_prueba$apto <- cobertura_prueba$cobertura >= 0.95 & cobertura_prueba$hueco <= 7 & cobertura_prueba$borde <= 7

cat("Prueba - errores:", length(error_prueba), "| acumulados:", length(acumulados_prueba),
    "| pares antes:", nrow(cobertura_prueba), "| pares aptos:", sum(cobertura_prueba$apto),
    "| estaciones aptas:", length(unique(cobertura_prueba$STATION_ID[cobertura_prueba$apto])), "\n")


##### Pentadas de las estaciones #####
# Datos diarios de los pares aptos (con los dias faltantes como NA)
pares_aptos <- paste(cobertura$STATION_ID, cobertura$YEAR)[cobertura$apto]
datos_tierra <- calendario[paste(calendario$STATION_ID, calendario$YEAR) %in% pares_aptos, ]

# Pentada de cada dia: 6 por mes, la sexta va del dia 26 al fin de mes
datos_tierra$MONTH <- as.numeric(format(datos_tierra$fecha, "%m"))
datos_tierra$DAY <- as.numeric(format(datos_tierra$fecha, "%d"))
datos_tierra$pentada <- (datos_tierra$MONTH - 1) * 6 + pmin(6, (datos_tierra$DAY - 1) %/% 5 + 1)
datos_tierra$STATION_ID <- as.character(datos_tierra$STATION_ID)
datos_tierra$curva <- paste(datos_tierra$STATION_ID, datos_tierra$YEAR, sep = "_")

# Lluvia por pentada: si falta algun dia, la pentada queda NA
tierra_pentadal <- datos_tierra %>%
  group_by(curva, STATION_ID, YEAR, pentada) %>%
  summarize(lluvia = sum(VALUE), .groups = "drop")

# Matriz terrestre: filas = pentadas, columnas = curvas estacion-año
Y_terrestre_total <- tapply(tierra_pentadal$lluvia, list(tierra_pentadal$pentada, tierra_pentadal$curva), sum)
dim(Y_terrestre_total)


##### Satelite y emparejamiento estacion-pixel #####
# Imagenes pentadales de 1981 a 2000
archivos <- list.files("data/ImagenesSatelitePentadal", pattern = "ppt.*\\.tif$", full.names = TRUE)
anio_archivo <- as.numeric(substr(basename(archivos), 4, 7))
archivos <- archivos[anio_archivo <= 2000]
table(substr(basename(archivos), 4, 7))   # 72 imagenes por año

p_sat <- rast(archivos)

# Estaciones con coordenadas y el pixel que les corresponde
curvas <- unique(tierra_pentadal[, c("curva", "STATION_ID", "YEAR")])
curvas <- merge(curvas, info_estaciones, by = "STATION_ID", all.x = TRUE)
curvas <- curvas[order(curvas$curva), ]

estaciones <- unique(curvas[!is.na(curvas$longitud), c("STATION_ID", "longitud", "latitud")])
puntos <- vect(estaciones, geom = c("longitud", "latitud"), crs = "EPSG:4326")
estaciones$pixel <- cellFromXY(p_sat, crds(puntos))

# Lluvia satelital en el pixel de cada estacion (filas = estaciones, columnas = imagenes)
lluvia_sat <- as.matrix(terra::extract(p_sat, puntos, ID = FALSE))
rownames(lluvia_sat) <- estaciones$STATION_ID

# Estado del emparejamiento de cada curva
curvas$pixel <- estaciones$pixel[match(curvas$STATION_ID, estaciones$STATION_ID)]
curvas$estado <- "Emparejado"
curvas$estado[is.na(curvas$longitud)] <- "Sin coordenadas"
curvas$estado[!is.na(curvas$longitud) & is.na(curvas$pixel)] <- "Fuera del raster"
for (i in which(curvas$estado == "Emparejado")) {
  columnas <- substr(colnames(lluvia_sat), 4, 7) == curvas$YEAR[i]
  if (any(is.na(lluvia_sat[curvas$STATION_ID[i], columnas]))) curvas$estado[i] <- "Pixel incompleto"
}
table(curvas$estado)

# Estaciones que comparten pixel
table(estaciones$pixel)[table(estaciones$pixel) > 1]

# Matrices alineadas: filas = pentadas, columnas = curvas emparejadas
metadatos <- curvas[curvas$estado == "Emparejado", ]

Y_est <- Y_terrestre_total[, metadatos$curva]
X_sat <- matrix(NA, 72, nrow(metadatos), dimnames = list(1:72, metadatos$curva))
for (i in 1:nrow(metadatos)) {
  columnas <- substr(colnames(lluvia_sat), 4, 7) == metadatos$YEAR[i]
  X_sat[, i] <- lluvia_sat[metadatos$STATION_ID[i], columnas]
}

dim(Y_est)
dim(X_sat)
sum(is.na(Y_est))
sum(!is.na(metadatos$altitud))


##### Matrices y particion #####
tiempos <- 1:72

# Transponer: filas = curvas estacion-año, columnas = pentadas
Y_matriz <- t(Y_est)
X_matriz <- t(X_sat)
colnames(Y_matriz) <- paste0("P", sprintf("%02d", tiempos))
colnames(X_matriz) <- paste0("P", sprintf("%02d", tiempos))

# Entrenamiento (1981-1995) y validacion (1996-2000)
Y_train <- Y_matriz[metadatos$YEAR <= 1995, ]
X_train <- X_matriz[metadatos$YEAR <= 1995, ]
Y_val <- Y_matriz[metadatos$YEAR >= 1996, ]
X_val <- X_matriz[metadatos$YEAR >= 1996, ]

nrow(Y_matriz)
nrow(Y_train)
nrow(Y_val)

# Guardar las matrices
saveRDS(list(Y_matriz = Y_matriz, X_matriz = X_matriz, metadatos = metadatos,
             Y_train = Y_train, X_train = X_train, Y_val = Y_val, X_val = X_val),
        "Resultados/matrices_funcionales.rds")

# Pentadas faltantes por curva
faltantes_curva <- rowSums(is.na(Y_matriz))
table(faltantes_curva)
faltantes_curva[faltantes_curva > 0]


##### Base B-spline y rango de la matriz #####
K_evaluar <- c(6, 8, 10, 12, 15, 20, 25, 30)

# Base B-spline cubica (orden 4) evaluada en las pentadas
construir_bspline <- function(t, K) {
  BasePrec <- create.bspline.basis(rangeval = c(1, 72), nbasis = K, norder = 4)
  eval.basis(t, BasePrec)
}

# Particion de la unidad
Phi_ejemplo <- construir_bspline(tiempos, 12)
dim(Phi_ejemplo)
range(rowSums(Phi_ejemplo))

# Rango de Phi en las pentadas observadas de cada curva
for (K in K_evaluar) {
  Phi <- construir_bspline(tiempos, K)
  rangos <- c()
  for (i in 1:nrow(Y_matriz)) {
    rangos[i] <- qr(Phi[!is.na(Y_matriz[i, ]), ], tol = 1e-9)$rank
  }
  cat("K =", K, "- curvas con rango completo:", sum(rangos == K), "de", nrow(Y_matriz), "\n")
}

# Curva con mas pentadas faltantes
i <- which.max(faltantes_curva)
rownames(Y_matriz)[i]
for (K in c(12, 15, 20, 25, 30)) {
  Phi_i <- construir_bspline(tiempos, K)[!is.na(Y_matriz[i, ]), ]
  cat("K =", K, "- rango:", qr(Phi_i, tol = 1e-9)$rank, "- condicion:", kappa(Phi_i, exact = TRUE), "\n")
}


##### Construccion por minimos cuadrados ordinarios (MCO) #####
# c = (Phi'Phi)^-1 Phi'y,  C = (Phi'Phi)^-1 Phi' (y2cMap),  S = Phi C
# Cada curva se ajusta solo con sus pentadas observadas
ajustar_MCO <- function(M, K) {
  Phi <- construir_bspline(tiempos, K)
  coeficientes <- matrix(NA, nrow(M), K, dimnames = list(rownames(M), NULL))
  predicciones <- matrix(NA, nrow(M), 72, dimnames = dimnames(M))
  residuos <- matrix(NA, nrow(M), 72, dimnames = dimnames(M))
  estado <- rep("ok", nrow(M))
  SSE <- rep(NA, nrow(M))
  gl <- rep(NA, nrow(M))
  
  for (i in 1:nrow(M)) {
    obs <- which(!is.na(M[i, ]))
    P <- Phi[obs, ]
    if (length(obs) <= K || kappa(P, exact = TRUE) > 1e4) {
      estado[i] <- "mal_condicionada"
      next
    }
    C <- solve(t(P) %*% P) %*% t(P)
    coeficientes[i, ] <- C %*% M[i, obs]
    predicciones[i, ] <- Phi %*% coeficientes[i, ]
    residuos[i, obs] <- M[i, obs] - predicciones[i, obs]
    SSE[i] <- sum(residuos[i, obs]^2)
    gl[i] <- sum(diag(P %*% C))
  }
  list(Phi = Phi, coeficientes = coeficientes, predicciones = predicciones,
       residuos = residuos, estado = estado, SSE = SSE, gl = gl)
}

# Ajuste de las dos fuentes para cada K
fuentes <- list(terrestre = Y_train, satelital = X_train)
resultados_MCO <- list()
for (fuente in names(fuentes)) {
  resultados_MCO[[fuente]] <- list()
  for (K in K_evaluar) {
    resultados_MCO[[fuente]][[as.character(K)]] <- ajustar_MCO(fuentes[[fuente]], K)
  }
}

# Curvas ajustadas por fuente y K
for (fuente in names(fuentes)) {
  for (K in K_evaluar) {
    a <- resultados_MCO[[fuente]][[as.character(K)]]
    cat(fuente, "- K =", K, "- curvas ajustadas:", sum(a$estado == "ok"), "\n")
  }
}

# Verificacion de la matriz suavizadora con una curva completa (K = 12)
i <- which(rowSums(is.na(Y_train)) == 0)[1]
Phi <- construir_bspline(tiempos, 12)
y <- Y_train[i, ]

y2cMap <- solve(t(Phi) %*% Phi) %*% t(Phi)
S <- Phi %*% y2cMap
c_est <- y2cMap %*% y

max(abs(S - t(S)))         # simetrica
max(abs(S %*% S - S))      # idempotente
sum(diag(S))               # traza = K

# Comparacion con smooth.basis
BasePrec <- create.bspline.basis(c(1, 72), nbasis = 12, norder = 4)
ajuste_fda <- smooth.basis(tiempos, y, BasePrec)
max(abs(c_est - ajuste_fda$fd$coefs))

y_sat <- X_train[which(rowSums(is.na(X_train)) == 0)[1], ]
max(abs(y2cMap %*% y_sat - smooth.basis(tiempos, y_sat, BasePrec)$fd$coefs))

# Grafica: efecto del numero de bases
par(mfrow = c(1, 2))
for (fuente in names(fuentes)) {
  plot(tiempos, fuentes[[fuente]][i, ], pch = 16, cex = 0.5,
       xlab = "Pentada", ylab = "Precipitacion (mm)", main = paste("Ajuste MCO -", fuente))
  lines(tiempos, resultados_MCO[[fuente]][["6"]]$predicciones[i, ], col = 2, lwd = 2)
  lines(tiempos, resultados_MCO[[fuente]][["12"]]$predicciones[i, ], col = 3, lwd = 2)
  lines(tiempos, resultados_MCO[[fuente]][["20"]]$predicciones[i, ], col = 4, lwd = 2)
  legend("topright", legend = c("K = 6", "K = 12", "K = 20"), col = 2:4, lty = 1, bty = "n")
}
par(mfrow = c(1, 1))


##### Diagnostico de residuos de MCO #####
# MCO supone residuos independientes y con varianza constante (Var[y] = sigma^2 I).
# Si los residuos muestran heterocedasticidad o correlacion, MCO no es adecuado
# y se justifica usar MCP o MCG.

# Varianza residual en temporada lluviosa frente al resto del año
lluviosas <- c(13:30, 49:66)
for (fuente in names(fuentes)) {
  for (K in c(6, 12, 20)) {
    R <- resultados_MCO[[fuente]][[as.character(K)]]$residuos
    v_lluviosa <- apply(R[, lluviosas], 1, var, na.rm = TRUE)
    v_otra <- apply(R[, -lluviosas], 1, var, na.rm = TRUE)
    razon <- (v_lluviosa / v_otra)[which(v_otra > 0)]
    cat(fuente, "- K =", K, "- razon mediana de varianzas:", round(median(razon), 3),
        "- % de curvas con mas varianza en temporada lluviosa:", round(100 * mean(razon > 1), 1), "\n")
  }
}

# Graficas de varianza y correlacion residual (K = 12)
par(mfrow = c(2, 2))
for (fuente in names(fuentes)) {
  R <- resultados_MCO[[fuente]][["12"]]$residuos
  plot(tiempos, apply(R, 2, var, na.rm = TRUE), type = "l", lwd = 2,
       xlab = "Pentada", ylab = "Varianza residual", main = paste("Varianza residual -", fuente))
  contour(tiempos, tiempos, cor(R, use = "pairwise.complete.obs"), nlevels = 10,
          xlab = "Pentada s", ylab = "Pentada t", main = paste("Correlacion residual -", fuente))
}
par(mfrow = c(1, 1))

# Prueba de independencia de los residuos de MCO
autocorrelacion_residuos <- function(datos, Phi) {
  completas <- datos[rowSums(is.na(datos)) == 0, ]
  S <- Phi %*% solve(t(Phi) %*% Phi) %*% t(Phi)
  r1 <- c()
  for (i in 1:nrow(completas)) {
    residuos <- completas[i, ] - S %*% completas[i, ]
    r1[i] <- acf(residuos, lag.max = 1, plot = FALSE)$acf[2]
  }
  r1
}

r1_estaciones <- autocorrelacion_residuos(Y_train, eval.basis(tiempos, create.fourier.basis(c(0, 72), nbasis = 5)))
r1_satelite <- autocorrelacion_residuos(X_train, eval.basis(tiempos, create.fourier.basis(c(0, 72), nbasis = 9)))

# Porcentaje de curvas con autocorrelacion significativa (se espera 5 % si hay independencia)
limite <- 1.96 / sqrt(72)
mean(abs(r1_estaciones) > limite) * 100
mean(abs(r1_satelite) > limite) * 100

# Autocorrelacion de rezago 1 de todas las curvas
par(mfrow = c(1, 2))
hist(r1_estaciones, breaks = 30, main = "Estaciones", xlab = "Autocorrelacion de rezago 1")
abline(v = c(-limite, limite), lty = 2, col = "red")
hist(r1_satelite, breaks = 30, main = "Satelite", xlab = "Autocorrelacion de rezago 1")
abline(v = c(-limite, limite), lty = 2, col = "red")
par(mfrow = c(1, 1))


# Heterocedasticidad de los residuos de MCO
residuos_MCO <- function(datos, Phi) {
  residuos <- matrix(NA, nrow(datos), 72)
  for (i in 1:nrow(datos)) {
    obs <- which(!is.na(datos[i, ]))
    P <- Phi[obs, ]
    residuos[i, obs] <- datos[i, obs] - P %*% solve(t(P) %*% P) %*% t(P) %*% datos[i, obs]
  }
  residuos
}

residuos_estaciones <- residuos_MCO(Y_train, eval.basis(tiempos, create.fourier.basis(c(0, 72), nbasis = 5)))
residuos_satelite <- residuos_MCO(X_train, eval.basis(tiempos, create.fourier.basis(c(0, 72), nbasis = 9)))

# Residuos de todas las curvas y su varianza en cada pentada
par(mfrow = c(2, 2))
matplot(tiempos, t(residuos_estaciones), type = "l", lty = 1, col = "grey", xlab = "Pentada", ylab = "Residuo",
        main = "Residuos MCO - estaciones")
matplot(tiempos, t(residuos_satelite), type = "l", lty = 1, col = "grey", xlab = "Pentada", ylab = "Residuo",
        main = "Residuos MCO - satelite")
plot(tiempos, apply(residuos_estaciones, 2, var, na.rm = TRUE), type = "b", pch = 19,
     xlab = "Pentada", ylab = "Varianza residual", main = "Varianza por pentada - estaciones")
plot(tiempos, apply(residuos_satelite, 2, var, na.rm = TRUE), type = "b", pch = 19,
     xlab = "Pentada", ylab = "Varianza residual", main = "Varianza por pentada - satelite")
par(mfrow = c(1, 1))

# Cuantas veces es mayor la varianza en la pentada mas variable que en la menos variable
varianza_estaciones <- apply(residuos_estaciones, 2, var, na.rm = TRUE)
varianza_satelite <- apply(residuos_satelite, 2, var, na.rm = TRUE)
max(varianza_estaciones) / min(varianza_estaciones)
max(varianza_satelite) / min(varianza_satelite)


# Precipitaciones negativas en los ajustes MCO
for (fuente in names(fuentes)) {
  for (K in K_evaluar) {
    P <- resultados_MCO[[fuente]][[as.character(K)]]$predicciones
    cat(fuente, "- K =", K, "- valores negativos:", sum(P < -1e-8, na.rm = TRUE),
        "- curvas afectadas:", sum(rowSums(P < -1e-8, na.rm = TRUE) > 0),
        "- minimo:", round(min(P, na.rm = TRUE), 2), "\n")
  }
}


##### Comparacion MCO, MCP y MCG #####
# MCP: W diagonal (varianzas distintas). MCG: W incluye covarianzas.
# Sigma se estima con los residuos de MCO del 70 % de las estaciones
# y se evalua en el 30 % restante, reservando una de cada 5 pentadas.
# La muestra usa la semilla fijada al inicio (set.seed(2026)).
estaciones_train <- sub("_[0-9]{4}$", "", rownames(Y_train))
est_ref <- sample(unique(estaciones_train), floor(0.7 * length(unique(estaciones_train))))
idx_ref <- estaciones_train %in% est_ref
pentadas_prueba <- seq(5, 72, by = 5)

# Minimos cuadrados con matriz Sigma (Sigma = NULL es MCO)
ajustar_metodo <- function(y, Phi, obs, Sigma = NULL) {
  P <- Phi[obs, ]
  z <- y[obs]
  if (!is.null(Sigma)) {
    L <- t(chol(Sigma[obs, obs]))
    z <- forwardsolve(L, z)
    P <- forwardsolve(L, P)
  }
  if (nrow(P) <= ncol(P) || kappa(P, exact = TRUE) > 1e4) return(NULL)
  coef <- solve(t(P) %*% P, t(P) %*% z)
  list(coef = as.vector(coef), pred = as.vector(Phi %*% coef))
}

# MCG con coeficientes c >= 0: como las B-splines son no negativas, la curva no es negativa
ajustar_MCG_positivo <- function(y, Phi, obs, Sigma) {
  K <- ncol(Phi)
  if (length(obs) <= K) return(list(estado = "observaciones_insuficientes"))
  L <- t(chol(Sigma[obs, obs]))
  Pw <- forwardsolve(L, Phi[obs, ])
  zw <- forwardsolve(L, y[obs])
  if (kappa(Pw, exact = TRUE) > 1e4) return(list(estado = "mal_condicionada"))
  
  escala <- max(1, sum(zw^2))
  objetivo <- function(c) sum((zw - Pw %*% c)^2) / escala
  gradiente <- function(c) as.vector(2 * t(Pw) %*% (Pw %*% c - zw) / escala)
  inicio <- pmax(0, as.vector(qr.solve(Pw, zw)))
  
  ajuste <- optim(inicio, objetivo, gradiente, method = "L-BFGS-B",
                  lower = rep(0, K), control = list(maxit = 3000))
  if (ajuste$convergence != 0) {
    ajuste <- nlminb(inicio, objetivo, gradiente, lower = rep(0, K),
                     control = list(iter.max = 3000, eval.max = 5000))
  }
  if (ajuste$convergence != 0) return(list(estado = "sin_convergencia"))
  
  coef <- pmax(0, ajuste$par)
  g <- gradiente(coef)
  g[coef <= 1e-8 & g > 0] <- 0
  if (max(abs(g)) > 1e-4) return(list(estado = "optimalidad_no_verificada"))
  list(estado = "ok", coef = coef, pred = as.vector(Phi %*% coef))
}

# ECM en las pentadas reservadas de las estaciones de evaluacion
comparacion_metodos <- data.frame()
sigmas <- list()
for (fuente in names(fuentes)) {
  M <- fuentes[[fuente]]
  sigmas[[fuente]] <- list()
  for (K in K_evaluar) {
    a <- resultados_MCO[[fuente]][[as.character(K)]]
    
    # Covarianza de los residuos MCO de las curvas completas de referencia
    R_ref <- a$residuos[idx_ref & rowSums(is.na(a$residuos)) == 0, ]
    Sigma <- cov(R_ref)
    eps <- 0.01 * median(diag(Sigma))
    Sigma_MCP <- diag(diag(Sigma)) + eps * diag(72)
    Sigma_MCG <- 0.5 * Sigma + 0.5 * diag(diag(Sigma)) + eps * diag(72)
    sigmas[[fuente]][[as.character(K)]] <- Sigma_MCG
    
    ecm <- c()
    negativos <- c(MCO = 0, MCP = 0, MCG = 0)
    for (i in which(!idx_ref)) {
      y <- M[i, ]
      prueba <- intersect(which(!is.na(y)), pentadas_prueba)
      obs <- setdiff(which(!is.na(y)), prueba)
      mco <- ajustar_metodo(y, a$Phi, obs)
      mcp <- ajustar_metodo(y, a$Phi, obs, Sigma_MCP)
      mcg <- ajustar_metodo(y, a$Phi, obs, Sigma_MCG)
      if (length(prueba) == 0 || is.null(mco) || is.null(mcp) || is.null(mcg)) next
      ecm <- rbind(ecm, c(MCO = mean((y[prueba] - mco$pred[prueba])^2),
                          MCP = mean((y[prueba] - mcp$pred[prueba])^2),
                          MCG = mean((y[prueba] - mcg$pred[prueba])^2)))
      negativos <- negativos + c(sum(mco$pred < -1e-8), sum(mcp$pred < -1e-8), sum(mcg$pred < -1e-8))
    }
    comparacion_metodos <- rbind(comparacion_metodos,
                                 data.frame(fuente = fuente, K = K, curvas = nrow(ecm),
                                            ECM_MCO = mean(ecm[, "MCO"]), ECM_MCP = mean(ecm[, "MCP"]),
                                            ECM_MCG = mean(ecm[, "MCG"]),
                                            neg_MCO = negativos[1], neg_MCP = negativos[2], neg_MCG = negativos[3]))
  }
}
comparacion_metodos$mejora_MCP <- 100 * (comparacion_metodos$ECM_MCO - comparacion_metodos$ECM_MCP) / comparacion_metodos$ECM_MCO
comparacion_metodos$mejora_MCG <- 100 * (comparacion_metodos$ECM_MCO - comparacion_metodos$ECM_MCG) / comparacion_metodos$ECM_MCO
rownames(comparacion_metodos) <- NULL
comparacion_metodos

# Propiedades de la matriz suavizadora de MCG (terrestre, K = 12)
Phi <- construir_bspline(tiempos, 12)
W <- solve(sigmas$terrestre[["12"]])
S_MCG <- Phi %*% solve(t(Phi) %*% W %*% Phi) %*% t(Phi) %*% W
max(abs(S_MCG %*% S_MCG - S_MCG))   # idempotente
sum(diag(S_MCG))                    # traza = K
max(abs(S_MCG - t(S_MCG)))          # no es simetrica


##### Construccion con MCG restringido (curvas no negativas) #####
verificacion_final <- data.frame()
comparacion_validacion <- data.frame()
for (fuente in names(fuentes)) {
  M <- fuentes[[fuente]]
  for (K in K_evaluar) {
    Phi <- construir_bspline(tiempos, K)
    Sigma <- sigmas[[fuente]][[as.character(K)]]
    
    # Todas las curvas
    coeficientes <- matrix(NA, nrow(M), K)
    predicciones <- matrix(NA, nrow(M), 72)
    for (i in 1:nrow(M)) {
      positivo <- ajustar_MCG_positivo(M[i, ], Phi, which(!is.na(M[i, ])), Sigma)
      if (positivo$estado == "ok") {
        coeficientes[i, ] <- positivo$coef
        predicciones[i, ] <- positivo$pred
      }
    }
    verificacion_final <- rbind(verificacion_final,
                                data.frame(fuente = fuente, K = K,
                                           ajustadas = sum(!is.na(coeficientes[, 1])),
                                           pendientes = sum(is.na(coeficientes[, 1])),
                                           coef_negativos = sum(coeficientes < -1e-8, na.rm = TRUE),
                                           pred_negativas = sum(predicciones < -1e-8, na.rm = TRUE),
                                           minimo = min(predicciones, na.rm = TRUE)))
    
    # Validacion en las pentadas reservadas: MCG frente a MCG restringido
    ecm <- c()
    for (i in which(!idx_ref)) {
      y <- M[i, ]
      prueba <- intersect(which(!is.na(y)), pentadas_prueba)
      obs <- setdiff(which(!is.na(y)), prueba)
      if (length(prueba) == 0) next
      mcg <- ajustar_metodo(y, Phi, obs, Sigma)
      positivo <- ajustar_MCG_positivo(y, Phi, obs, Sigma)
      if (is.null(mcg) || positivo$estado != "ok") next
      ecm <- rbind(ecm, c(mean((y[prueba] - mcg$pred[prueba])^2),
                          mean((y[prueba] - positivo$pred[prueba])^2)))
    }
    comparacion_validacion <- rbind(comparacion_validacion,
                                    data.frame(fuente = fuente, K = K, curvas = nrow(ecm),
                                               ECM_MCG = mean(ecm[, 1]), ECM_restringido = mean(ecm[, 2])))
  }
}
verificacion_final
comparacion_validacion$mejora_pct <- 100 * (comparacion_validacion$ECM_MCG - comparacion_validacion$ECM_restringido) /
  comparacion_validacion$ECM_MCG
comparacion_validacion

# Grafica: MCG frente a MCG restringido en una curva (K = 12)
Phi <- construir_bspline(tiempos, 12)
par(mfrow = c(1, 2))
for (fuente in names(fuentes)) {
  y <- fuentes[[fuente]][1, ]
  obs <- which(!is.na(y))
  mcg <- ajustar_metodo(y, Phi, obs, sigmas[[fuente]][["12"]])
  positivo <- ajustar_MCG_positivo(y, Phi, obs, sigmas[[fuente]][["12"]])
  plot(tiempos, y, pch = 16, cex = 0.5, xlab = "Pentada", ylab = "Precipitacion (mm)",
       main = paste("MCG vs MCG restringido -", fuente))
  lines(tiempos, mcg$pred, col = "blue", lwd = 2)
  lines(tiempos, positivo$pred, col = "red", lwd = 2)
  abline(h = 0, lty = 2)
  legend("topright", c("MCG", "MCG restringido"), col = c("blue", "red"), lty = 1, bty = "n")
}
par(mfrow = c(1, 1))

##### Seleccion de la base, el numero de bases y el metodo #####
criterios <- function(datos, Phi, pesos = rep(1, 72), R = 0, lambda = 0) {
  ECM <- c()
  VC <- c()
  VCG <- c()
  for (i in 1:nrow(datos)) {
    obs <- which(!is.na(datos[i, ]))
    y <- datos[i, obs]
    P <- Phi[obs, ]
    PW <- t(P * pesos[obs])                  # Phi' W
    S <- P %*% solve(PW %*% P + lambda * R) %*% PW
    error <- y - S %*% y
    n <- length(obs)
    ECM[i] <- mean(error^2)
    VC[i] <- mean((error / (1 - diag(S)))^2)
    VCG[i] <- n * sum(error^2) / (n - sum(diag(S)))^2
  }
  c(ECM = mean(ECM), VC = mean(VC), VCG = sum(VCG))
}

# Pesos de MCP
pesos_MCP <- function(datos, Phi) {
  residuos <- matrix(NA, nrow(datos), 72)
  for (i in 1:nrow(datos)) {
    obs <- which(!is.na(datos[i, ]))
    P <- Phi[obs, ]
    residuos[i, obs] <- datos[i, obs] - P %*% solve(t(P) %*% P) %*% t(P) %*% datos[i, obs]
  }
  1 / colMeans(residuos^2, na.rm = TRUE)
}

crear_base <- function(nombre, K) {
  if (nombre == "B-spline") return(create.bspline.basis(c(0, 72), nbasis = K, norder = 4))
  if (nombre == "Fourier") return(create.fourier.basis(c(0, 72), nbasis = K))
  create.monomial.basis(c(0, 72), nbasis = K)
}

# B-spline cubica
valores_K <- list("B-spline" = 4:15, "Fourier" = seq(3, 15, 2), "Monomial" = 1:4)

# MCO y MCP para cada base y K
tabla <- c()
for (nombre in names(valores_K)) {
  for (K in valores_K[[nombre]]) {
    Phi <- eval.basis(tiempos, crear_base(nombre, K))
    tabla <- rbind(tabla,
                   data.frame(fuente = "estaciones", base = nombre, metodo = "MCO", parametro = K,
                              t(criterios(Y_train, Phi))),
                   data.frame(fuente = "estaciones", base = nombre, metodo = "MCP", parametro = K,
                              t(criterios(Y_train, Phi, pesos_MCP(Y_train, Phi)))),
                   data.frame(fuente = "satelite", base = nombre, metodo = "MCO", parametro = K,
                              t(criterios(X_train, Phi))),
                   data.frame(fuente = "satelite", base = nombre, metodo = "MCP", parametro = K,
                              t(criterios(X_train, Phi, pesos_MCP(X_train, Phi)))))
  }
}

# Penalizacion por rugosidad con una base amplia; lambda se elige con la VC (parametro = log10(lambda))
# B-spline: K = 72 y penalizacion de la segunda derivada
# Fourier: K = 71 y aceleracion armonica (no penaliza el ciclo anual)
base_bspline <- create.bspline.basis(c(0, 72), nbasis = 72, norder = 4)
base_fourier <- create.fourier.basis(c(0, 72), nbasis = 71)
penalizadas <- list("B-spline" = list(Phi = eval.basis(tiempos, base_bspline),
                                      R = eval.penalty(base_bspline, 2)),
                    "Fourier" = list(Phi = eval.basis(tiempos, base_fourier),
                                     R = eval.penalty(base_fourier, vec2Lfd(c(0, (2 * pi / 72)^2, 0), c(0, 72)))))

for (nombre in names(penalizadas)) {
  Phi <- penalizadas[[nombre]]$Phi
  R <- penalizadas[[nombre]]$R
  for (log_lambda in seq(-1, 5, 0.5)) {
    tabla <- rbind(tabla,
                   data.frame(fuente = "estaciones", base = nombre, metodo = "Penalizacion", parametro = log_lambda,
                              t(criterios(Y_train, Phi, R = R, lambda = 10^log_lambda))),
                   data.frame(fuente = "satelite", base = nombre, metodo = "Penalizacion", parametro = log_lambda,
                              t(criterios(X_train, Phi, R = R, lambda = 10^log_lambda))))
  }
}
tabla

# Tabla de ablacion: la mejor configuracion (menor VC) de cada fuente, base y metodo
ablacion <- tabla[order(tabla$fuente, tabla$VC), ]
ablacion <- ablacion[!duplicated(ablacion[, c("fuente", "base", "metodo")]), ]
ablacion

# Graficas de VC: MCO y MCP para cada base (estaciones arriba, satelite abajo)
par(mfrow = c(2, 3))
for (f in c("estaciones", "satelite")) {
  for (nombre in names(valores_K)) {
    mco <- tabla[tabla$fuente == f & tabla$base == nombre & tabla$metodo == "MCO", ]
    mcp <- tabla[tabla$fuente == f & tabla$base == nombre & tabla$metodo == "MCP", ]
    plot(mco$parametro, mco$VC, type = "b", pch = 19, xlab = "K", ylab = "VC",
         main = paste(nombre, "-", f), ylim = range(mco$VC, mcp$VC))
    lines(mcp$parametro, mcp$VC, type = "b", pch = 17, col = "red")
    legend("topright", c("MCO", "MCP"), col = c("black", "red"), pch = c(19, 17), bty = "n")
  }
}
par(mfrow = c(1, 1))

# Graficas de VC de la penalizacion
par(mfrow = c(1, 2))
for (f in c("estaciones", "satelite")) {
  b <- tabla[tabla$fuente == f & tabla$base == "B-spline" & tabla$metodo == "Penalizacion", ]
  fo <- tabla[tabla$fuente == f & tabla$base == "Fourier" & tabla$metodo == "Penalizacion", ]
  plot(b$parametro, b$VC, type = "b", pch = 19, xlab = expression(log[10](lambda)), ylab = "VC",
       main = paste("Penalizacion -", f), ylim = range(b$VC, fo$VC))
  lines(fo$parametro, fo$VC, type = "b", pch = 17, col = "red")
  legend("topright", c("B-spline", "Fourier"), col = c("black", "red"), pch = c(19, 17), bty = "n")
}
par(mfrow = c(1, 1))


# Validacion cruzada: el atajo con la diagonal de S da lo mismo que dejar uno por fuera
Phi <- eval.basis(tiempos, crear_base("B-spline", 6))
y <- Y_train[which(rowSums(is.na(Y_train)) == 0)[1], ]

S <- Phi %*% solve(t(Phi) %*% Phi) %*% t(Phi)
VC_atajo <- mean(((y - S %*% y) / (1 - diag(S)))^2)

prediccion <- c()
for (j in 1:72) {
  coeficientes <- solve(t(Phi[-j, ]) %*% Phi[-j, ]) %*% t(Phi[-j, ]) %*% y[-j]
  prediccion[j] <- Phi[j, ] %*% coeficientes
}
VC_dejar_uno <- mean((y - prediccion)^2)

VC_atajo
VC_dejar_uno

# Valores negativos en las curvas ajustadas 
ajuste_MCP <- function(datos, Phi, pesos) {
  ajustes <- matrix(NA, nrow(datos), 72)
  for (i in 1:nrow(datos)) {
    obs <- which(!is.na(datos[i, ]))
    P <- Phi[obs, ]
    ajustes[i, ] <- Phi %*% solve(t(P * pesos[obs]) %*% P) %*% t(P * pesos[obs]) %*% datos[i, obs]
  }
  ajustes
}

Phi_estaciones <- eval.basis(tiempos, create.fourier.basis(c(0, 72), nbasis = 5))
Phi_satelite <- eval.basis(tiempos, create.fourier.basis(c(0, 72), nbasis = 9))
ajustes_estaciones <- ajuste_MCP(Y_train, Phi_estaciones, 1 / apply(residuos_estaciones^2, 2, mean, na.rm = TRUE))
ajustes_satelite <- ajuste_MCP(X_train, Phi_satelite, 1 / apply(residuos_satelite^2, 2, mean, na.rm = TRUE))

# Porcentaje de valores negativos, curvas afectadas y valor minimo
mean(ajustes_estaciones < 0) * 100
sum(rowSums(ajustes_estaciones < 0) > 0)
min(ajustes_estaciones)
mean(ajustes_satelite < 0) * 100
sum(rowSums(ajustes_satelite < 0) > 0)
min(ajustes_satelite)

# La diferencia entre penalizar y no penalizar es minima (menos del 1 % en la VC).
# Fourier es la base con mejor rendimiento en todas las comparaciones.
# MCO se rechaza: la varianza de los residuos cambia a lo largo del año
# (mayor en las temporadas lluviosas). La autocorrelacion de los residuos es debil.
# Se escoge MCP, que corrige la heterocedasticidad. MCG no se usa porque la Sigma
# estimada con los residuos no es invertible.
# Se escoge K = 5 para estaciones y K = 9 para satelite (menor VC con Fourier y MCP).

# El 0.95% de los datos de estacion tiene valores menores a cero.




