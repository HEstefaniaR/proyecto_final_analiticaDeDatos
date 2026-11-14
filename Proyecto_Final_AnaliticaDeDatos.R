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

# Conteo de registros validos, estaciones y pares estacion-año en cada paso del preprocesamiento
contar <- function(datos, paso) {
  validos <- datos[!is.na(datos$VALUE), ]
  data.frame(paso = paso, registros = nrow(validos),
             estaciones = length(unique(validos$STATION_ID)),
             pares_estacion_anio = nrow(unique(validos[, c("STATION_ID", "YEAR")])))
}

trazabilidad_pasos <- contar(p_tierra, "Datos 1981-2003")


# Separar datos
p_prueba <- p_tierra[p_tierra$YEAR >= 2001, ]
p_tierra <- p_tierra[p_tierra$YEAR <= 2000, ]

trazabilidad_pasos <- rbind(trazabilidad_pasos, contar(p_tierra, "Sin el conjunto de prueba (1981-2000)"))
trazabilidad_prueba <- contar(p_prueba, "Datos 2001-2003")

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

trazabilidad_pasos <- rbind(trazabilidad_pasos, contar(p_tierra, "Errores aislados"))

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

trazabilidad_pasos <- rbind(trazabilidad_pasos, contar(p_m_depurado, "Acumulados mensuales"))

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

trazabilidad_pasos <- rbind(trazabilidad_pasos, contar(p_m_depurado, "Ceros falsos"))

##### Trazabilidad #####
trazabilidad <- rbind(data.frame(p_m[acumulados, c("STATION_ID", "fecha", "VALUE")],
                                 motivo = "Precipitacion final >= 100 mm y >= 80% del total mensual"),
                      data.frame(p_m[ceros_falsos, c("STATION_ID", "fecha", "VALUE")],
                                 motivo = "Mes con ceros sospechosos"))
trazabilidad

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

trazabilidad_pasos <- rbind(trazabilidad_pasos, 
                            contar(p_m_depurado[paste(p_m_depurado$STATION_ID, 
                                                      p_m_depurado$YEAR) %in% 
                                                  paste(cobertura$STATION_ID, cobertura$YEAR)[cobertura$apto], ], "Cobertura minima"))

cat("Pares estacion-año antes:", nrow(cobertura), "| aptos:", sum(cobertura$apto), "\n")
cat("Estaciones antes:", length(unique(cobertura$STATION_ID)),
    "| despues:", length(unique(cobertura$STATION_ID[cobertura$apto])), "\n")


##### Tabla ONI #####

tabla_oni <- read_excel("data/ONI_1981_2003.xlsx")
names(tabla_oni)[-1] <- c("DEF", "EFM", "FMA", "MAM", "AMJ", "MJJ",
                          "JJA", "JAS", "ASO", "SON", "OND", "NDE")

# Ventanas trimestrales ordenadas por su mes central. Ejemplo: DEF = enero, NDE = diciembre
ventanas_oni <- names(tabla_oni)[-1]
oni_temporal <- data.frame(
  YEAR = rep(tabla_oni$Year, each = length(ventanas_oni)),
  mes_central = rep(seq_along(ventanas_oni), times = nrow(tabla_oni)),
  ventana = rep(ventanas_oni, times = nrow(tabla_oni)),
  ONI = as.vector(t(as.matrix(tabla_oni[, -1])))
) %>%
  arrange(YEAR, mes_central) %>%
  mutate(fase_umbral = case_when(
    ONI >= 0.5 ~ "El Niño",
    ONI <= -0.5 ~ "La Niña",
    TRUE ~ "Neutral"
  ))

# Rachas continuas entre años, sin reiniciar el conteo en enero
rachas_oni <- rle(oni_temporal$fase_umbral)
oni_temporal$racha <- rep(seq_along(rachas_oni$lengths), rachas_oni$lengths)
oni_temporal$ventanas_racha <- rep(rachas_oni$lengths, rachas_oni$lengths)

# Criterio histórico NOAA - al menos 5 ventanas consecutivas del mismo signo
oni_temporal <- oni_temporal %>%
  mutate(fase_enso = ifelse(ventanas_racha >= 5, fase_umbral, "Neutral"))

# Inicio y fin expresados como ventanas trimestrales
episodios_enso <- oni_temporal %>%
  filter(fase_enso != "Neutral") %>%
  group_by(racha, fase_enso) %>%
  summarise(inicio = paste(first(ventana), first(YEAR)),
            fin = paste(last(ventana), last(YEAR)),
            ventanas = n(),
            .groups = "drop") %>%
  mutate(episodio = row_number()) %>%
  select(episodio, fase_enso, inicio, fin, ventanas)

# Fases dentro de cada año, asignadas por el mes central de cada ventana
enso_anual <- oni_temporal %>%
  group_by(YEAR) %>%
  summarise(ventanas_nino = sum(fase_enso == "El Niño"),
            ventanas_nina = sum(fase_enso == "La Niña"),
            ventanas_neutral = sum(fase_enso == "Neutral"),
            secuencia_fases = paste(rle(fase_enso)$values, collapse = "-> "),
            cambio_fase = n_distinct(fase_enso) > 1,
            nino_y_nina = ventanas_nino > 0 & ventanas_nina > 0,
            .groups = "drop")

anios_cambio_fase <- enso_anual %>% filter(cambio_fase)

# Resúmenes anuales de ONI para comparar como entradas del modelo
oni_anual <- oni_temporal %>%
  group_by(YEAR) %>%
  summarise(oni_promedio = mean(ONI),
            oni_max_intensidad = ONI[which.max(abs(ONI))],
            .groups = "drop")
print(oni_anual, n=Inf)

# Fase ENSO por año y mes central para evaluar tramos de las curvas
enso_evaluacion <- oni_temporal %>%
  select(YEAR, mes_central, ventana, ONI, fase_enso)

head(enso_evaluacion, n=12)

write.csv(oni_anual, "data/oni_anual.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(enso_evaluacion, "data/enso_evaluacion.csv", row.names = FALSE, fileEncoding = "UTF-8")

print(episodios_enso, n = Inf)
print(enso_anual, n = Inf, width = Inf)


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

trazabilidad_prueba <- rbind(trazabilidad_prueba, contar(p_prueba, "Errores aislados"))

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

trazabilidad_prueba <- rbind(trazabilidad_prueba, 
                             contar(p_prueba, "Acumulados mensuales"))

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

trazabilidad_prueba <- rbind(trazabilidad_prueba, 
                             contar(p_prueba[paste(p_prueba$STATION_ID, 
                                                   p_prueba$YEAR) %in% 
                                               paste(cobertura_prueba$STATION_ID, 
                                                     cobertura_prueba$YEAR)[cobertura_prueba$apto], ], "Cobertura minima"))

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

# Trazabilidad
trazabilidad_pasos <- rbind(trazabilidad_pasos,
                            contar(datos_tierra[datos_tierra$curva %in% metadatos$curva, ], "Emparejamiento estacion-pixel"))
trazabilidad_pasos
trazabilidad_prueba

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

# Residuos de MCO de cada curva
residuos_MCO <- function(datos, Phi) {
  residuos <- matrix(NA, nrow(datos), 72)
  for (i in 1:nrow(datos)) {
    obs <- which(!is.na(datos[i, ]))
    P <- Phi[obs, ]
    residuos[i, obs] <- datos[i, obs] - P %*% solve(t(P) %*% P) %*% t(P) %*% datos[i, obs]
  }
  residuos
}

# Pesos de MCP
pesos_MCP <- function(datos, Phi) {
  1 / colMeans(residuos_MCO(datos, Phi)^2, na.rm = TRUE)
}

crear_base <- function(nombre, K) {
  if (nombre == "B-spline") return(create.bspline.basis(c(0, 72), nbasis = K, norder = 4))
  if (nombre == "Fourier") return(create.fourier.basis(c(0, 72), nbasis = K))
  create.monomial.basis(c(0, 72), nbasis = K)
}

# B-spline cubica: K de 4 a 15 (con K mayor algunas curvas con NA interpolan, s_ii = 1)
# Fourier: K impar (K = 2M + 1).  Monomial: K de 1 a 4 (con mas, Phi'Phi no se puede invertir)
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

# Penalizacion por rugosidad con una base amplia
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

# Tabla de ablacion
ablacion <- tabla[order(tabla$fuente, tabla$VC), ]
ablacion <- ablacion[!duplicated(ablacion[, c("fuente", "base", "metodo")]), ]
ablacion

# Graficas de VC
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

# Base elegida
# ECM, VC y VCG frente a K
estaciones_MCP <- tabla[tabla$fuente == "estaciones" & tabla$base == "B-spline" & tabla$metodo == "MCP", ]
satelite_MCP <- tabla[tabla$fuente == "satelite" & tabla$base == "B-spline" & tabla$metodo == "MCP", ]
K_estaciones <- estaciones_MCP$parametro[which.min(estaciones_MCP$VC)]
K_satelite <- satelite_MCP$parametro[which.min(satelite_MCP$VC)]
K_estaciones
K_satelite

par(mfrow = c(2, 3))
for (criterio in c("ECM", "VC", "VCG")) {
  plot(estaciones_MCP$parametro, estaciones_MCP[[criterio]], type = "b", pch = 19, xlab = "K", ylab = criterio,
       main = paste(criterio, "- estaciones"))
  abline(v = K_estaciones, lty = 3, col = "red")
}
for (criterio in c("ECM", "VC", "VCG")) {
  plot(satelite_MCP$parametro, satelite_MCP[[criterio]], type = "b", pch = 19, xlab = "K", ylab = criterio,
       main = paste(criterio, "- satelite"))
  abline(v = K_satelite, lty = 3, col = "red")
}
par(mfrow = c(1, 1))

Phi_estaciones <- eval.basis(tiempos, crear_base("B-spline", K_estaciones))
Phi_satelite <- eval.basis(tiempos, crear_base("B-spline", K_satelite))

# Validacion cruzada
y <- Y_train[which(rowSums(is.na(Y_train)) == 0)[1], ]
S <- Phi_estaciones %*% solve(t(Phi_estaciones) %*% Phi_estaciones) %*% t(Phi_estaciones)
VC_atajo <- mean(((y - S %*% y) / (1 - diag(S)))^2)

prediccion <- c()
for (j in 1:72) {
  coeficientes <- solve(t(Phi_estaciones[-j, ]) %*% Phi_estaciones[-j, ]) %*% t(Phi_estaciones[-j, ]) %*% y[-j]
  prediccion[j] <- Phi_estaciones[j, ] %*% coeficientes
}
VC_dejar_uno <- mean((y - prediccion)^2)

VC_atajo
VC_dejar_uno

# Matriz suavizadora de MCP
pesos_estaciones <- pesos_MCP(Y_train, Phi_estaciones)
pesos_satelite <- pesos_MCP(X_train, Phi_satelite)

S_estaciones <- Phi_estaciones %*% solve(t(Phi_estaciones * pesos_estaciones) %*% Phi_estaciones) %*% t(Phi_estaciones * pesos_estaciones)
max(abs(S_estaciones %*% S_estaciones - S_estaciones))
sum(diag(S_estaciones))

S_satelite <- Phi_satelite %*% solve(t(Phi_satelite * pesos_satelite) %*% Phi_satelite) %*% t(Phi_satelite * pesos_satelite)
max(abs(S_satelite %*% S_satelite - S_satelite))
sum(diag(S_satelite))


##### Residuos de MCO #####
# MCO supone residuos independientes y con varianza constante 
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

r1_estaciones <- autocorrelacion_residuos(Y_train, Phi_estaciones)
r1_satelite <- autocorrelacion_residuos(X_train, Phi_satelite)

limite <- 1.96 / sqrt(72)
mean(abs(r1_estaciones) > limite) * 100
mean(abs(r1_satelite) > limite) * 100

par(mfrow = c(1, 2))
hist(r1_estaciones, breaks = 30, main = "Estaciones", xlab = "Autocorrelacion de rezago 1")
abline(v = c(-limite, limite), lty = 2, col = "red")
hist(r1_satelite, breaks = 30, main = "Satelite", xlab = "Autocorrelacion de rezago 1")
abline(v = c(-limite, limite), lty = 2, col = "red")
par(mfrow = c(1, 1))

# Heterocedasticidad: varianza de los residuos en cada pentada
residuos_estaciones <- residuos_MCO(Y_train, Phi_estaciones)
residuos_satelite <- residuos_MCO(X_train, Phi_satelite)
varianza_estaciones <- apply(residuos_estaciones, 2, var, na.rm = TRUE)
varianza_satelite <- apply(residuos_satelite, 2, var, na.rm = TRUE)

par(mfrow = c(2, 2))
matplot(tiempos, t(residuos_estaciones), type = "l", lty = 1, col = "grey", xlab = "Pentada", ylab = "Residuo",
        main = "Residuos MCO - estaciones")
matplot(tiempos, t(residuos_satelite), type = "l", lty = 1, col = "grey", xlab = "Pentada", ylab = "Residuo",
        main = "Residuos MCO - satelite")
plot(tiempos, varianza_estaciones, type = "b", pch = 19, xlab = "Pentada", ylab = "Varianza residual",
     main = "Varianza por pentada - estaciones")
plot(tiempos, varianza_satelite, type = "b", pch = 19, xlab = "Pentada", ylab = "Varianza residual",
     main = "Varianza por pentada - satelite")
par(mfrow = c(1, 1))

max(varianza_estaciones) / min(varianza_estaciones)
max(varianza_satelite) / min(varianza_satelite)

##### Valores negativos en las curvas ajustadas con MCP #####
ajuste_MCP <- function(datos, Phi, pesos) {
  ajustes <- matrix(NA, nrow(datos), 72)
  for (i in 1:nrow(datos)) {
    obs <- which(!is.na(datos[i, ]))
    P <- Phi[obs, ]
    ajustes[i, ] <- Phi %*% solve(t(P * pesos[obs]) %*% P) %*% t(P * pesos[obs]) %*% datos[i, obs]
  }
  ajustes
}

ajustes_estaciones <- ajuste_MCP(Y_train, Phi_estaciones, pesos_estaciones)
ajustes_satelite <- ajuste_MCP(X_train, Phi_satelite, pesos_satelite)

# Porcentaje de valores negativos, curvas afectadas y valor minimo
mean(ajustes_estaciones < 0) * 100
sum(rowSums(ajustes_estaciones < 0) > 0)
min(ajustes_estaciones)
mean(ajustes_satelite < 0) * 100
sum(rowSums(ajustes_satelite < 0) > 0)
min(ajustes_satelite)

# Fourier tiene una VC un poco menor, pero se usa B-spline porque el codificador de la red
# exige curvas y pesos en B-splines (z = A W c).
# La diferencia entre penalizar y no penalizar es minima (menos del 1 % en la VC).
# MCO se rechaza: la varianza de los residuos cambia a lo largo del año (mayor en las
# temporadas lluviosas). La autocorrelacion de los residuos es debil.
# Se escoge MCP


##### Varianza y bandas #####
# Varianza residual de referencia MCO
residuos_referencia <- list(terrestre = residuos_estaciones, satelital = residuos_satelite)
bases_referencia <- list(terrestre = Phi_estaciones, satelital = Phi_satelite)
varianza_residual_MCO <- list()

# SCE / (n - K) por curva con rango completo y grados de libertad positivos
for (fuente in names(fuentes)) {
  datos <- fuentes[[fuente]]
  R_mco <- residuos_referencia[[fuente]]
  Phi_ref <- bases_referencia[[fuente]]
  n_obs <- rowSums(!is.na(datos))
  K_ref <- ncol(Phi_ref)
  gl_residual <- n_obs - K_ref
  SSE <- rowSums(R_mco^2, na.rm = TRUE)
  aptas <- gl_residual > 0 & rowSums(is.finite(R_mco)) == n_obs
  for (i in which(aptas)) {
    obs <- !is.na(datos[i, ])
    aptas[i] <- qr(Phi_ref[obs, , drop = FALSE], tol = 1e-9)$rank == K_ref
  }
  
  gl_residual[!aptas] <- NA_integer_
  SSE[!aptas] <- NA_real_
  sigma2 <- SSE / gl_residual
  varianza_residual_MCO[[fuente]] <- data.frame(curva = rownames(datos), n_obs, K = K_ref,
                                                gl_residual, SSE, sigma2)
  cat(fuente, "- estimaciones de sigma^2:", sum(aptas), "de", length(aptas), "curvas\n")
  print(summary(sigma2))
}


# Mapa de coeficientes y covarianza del MCP
mapas_MCP <- list()
for (fuente in names(fuentes)) {
  datos_mcp <- fuentes[[fuente]]
  Phi_mcp <- bases_referencia[[fuente]]
  R_mcp <- residuos_referencia[[fuente]]
  
  # Diagonal de Sigma estimada con el segundo momento residual por pentada
  varianzas_mcp <- apply(R_mcp^2, 2, mean, na.rm = TRUE)
  n_residuos_mcp <- colSums(is.finite(R_mcp))
  pesos_mcp <- 1 / varianzas_mcp
  pesos_validos <- is.finite(pesos_mcp) & pesos_mcp > 0
  varianzas_mcp[!pesos_validos] <- NA_real_
  pesos_mcp[!pesos_validos] <- NA_real_
  Sigma_mcp <- diag(varianzas_mcp)
  
  curvas_mcp <- vector("list", nrow(datos_mcp))
  names(curvas_mcp) <- rownames(datos_mcp)
  for (i in seq_len(nrow(datos_mcp))) {
    obs <- which(!is.na(datos_mcp[i, ]))
    if (length(obs) < ncol(Phi_mcp) || any(!pesos_validos[obs])) next
    P <- Phi_mcp[obs, , drop = FALSE]
    if (qr(P, tol = 1e-9)$rank < ncol(Phi_mcp)) next
    
    # C = (P' W P)^-1 P' W sobre las pentadas observadas
    PW <- t(P * pesos_mcp[obs])
    C_mcp <- solve(PW %*% P, PW)
    curvas_mcp[[i]] <- list(obs = obs, C = C_mcp, varianzas_obs = varianzas_mcp[obs])
  }
  
  mapas_MCP[[fuente]] <- list(Phi = Phi_mcp, pesos = pesos_mcp, Sigma = Sigma_mcp,
                              n_residuos = n_residuos_mcp, curvas = curvas_mcp)
  cat(fuente, "- mapas C:", sum(lengths(curvas_mcp) > 0), "de", length(curvas_mcp),
      "curvas | pentadas sin peso:", sum(!pesos_validos), "\n")
}


# Distribucion de la varianza residual de referencia por fuente
varianzas_grafico <- lapply(varianza_residual_MCO, function(tabla) {
  tabla$sigma2[is.finite(tabla$sigma2) & tabla$sigma2 > 0]
})
if (any(lengths(varianzas_grafico) > 0)) {
  boxplot(varianzas_grafico, names = c("Estaciones", "Satélite"), log = "y",
          col = c("lightblue", "moccasin"),
          main = "Varianza residual por fuente", ylab = "Varianza (mm², escala log)")
}

# Comparacion de las varianzas para la misma curva estacion-año
varianzas_pares <- merge(varianza_residual_MCO$terrestre[, c("curva", "sigma2")],
                         varianza_residual_MCO$satelital[, c("curva", "sigma2")],
                         by = "curva", suffixes = c("_terrestre", "_satelital"))
pares_positivos <- with(varianzas_pares, is.finite(sigma2_terrestre) & sigma2_terrestre > 0 &
                          is.finite(sigma2_satelital) & sigma2_satelital > 0)
varianzas_pares <- varianzas_pares[pares_positivos, ]

if (nrow(varianzas_pares) > 0) {
  limites_varianza <- range(varianzas_pares$sigma2_terrestre, varianzas_pares$sigma2_satelital)
  plot(varianzas_pares$sigma2_terrestre, varianzas_pares$sigma2_satelital,
       log = "xy", xlim = limites_varianza, ylim = limites_varianza, pch = 16, cex = 0.7,
       col = "steelblue",
       main = "Varianza residual por par estación-satélite",
       xlab = "Estaciones (mm², escala log)", ylab = "Satélite (mm², escala log)")
  abline(0, 1, lty = 2, col = "darkorange")
}


# Varianza de la curva ajustada con MCP en cada pentada
varianza_curvas_MCP <- list()
for (fuente in names(mapas_MCP)) {
  mapa <- mapas_MCP[[fuente]]
  datos_mcp <- fuentes[[fuente]]
  varianzas_curvas <- matrix(NA_real_, nrow(datos_mcp), ncol(datos_mcp),
                             dimnames = dimnames(datos_mcp))
  
  for (i in which(lengths(mapa$curvas) > 0)) {
    curva <- mapa$curvas[[i]]
    Phi_C <- mapa$Phi %*% curva$C
    # Diagonal de Phi C Sigma C' Phi' cuando Sigma es diagonal
    varianzas_curvas[i, ] <- Phi_C^2 %*% curva$varianzas_obs
  }
  
  varianza_curvas_MCP[[fuente]] <- varianzas_curvas
}


# Gráficos Q-Q de los residuos MCP en las péntadas 2, 40 y 60
residuos_MCP <- list(terrestre = Y_train - ajustes_estaciones,
                     satelital = X_train - ajustes_satelite)
colores_bandas <- c(terrestre = "steelblue", satelital = "darkorange")
etiquetas_bandas <- c(terrestre = "Estación", satelital = "Satélite")
pentadas_QQ <- c(2, 40, 60)
par_QQ <- par(mfrow = c(2, 3))
for (fuente in names(residuos_MCP)) {
  for (j in pentadas_QQ) {
    residuos_pentada <- residuos_MCP[[fuente]][, j]
    qqnorm(residuos_pentada, pch = 16, cex = 0.5,
           main = paste(etiquetas_bandas[[fuente]], "- Péntada", j),
           xlab = "Cuantiles normales", ylab = "Residuos MCP (mm)")
    qqline(residuos_pentada, col = colores_bandas[[fuente]], lwd = 2)
  }
}
par(par_QQ)


# Comparación de pares estación-año para las bandas
tabla_pares_bandas <- data.frame(
  curva = rownames(Y_train),
  observadas_tierra = rowSums(!is.na(Y_train)),
  observadas_satelite = rowSums(!is.na(X_train)),
  acumulado_tierra = rowSums(ajustes_estaciones),
  acumulado_satelite = rowSums(ajustes_satelite),
  RMSE_tierra = sqrt(rowMeans((Y_train - ajustes_estaciones)^2, na.rm = TRUE))
)

# Distancia al percentil 50 del acumulado ajustado en las dos fuentes
tabla_pares_bandas$percentil_tierra <- 100 * (rank(tabla_pares_bandas$acumulado_tierra) - 0.5) /
  nrow(tabla_pares_bandas)
tabla_pares_bandas$percentil_satelite <- 100 * (rank(tabla_pares_bandas$acumulado_satelite) - 0.5) /
  nrow(tabla_pares_bandas)
tabla_pares_bandas$distancia_mediana <- abs(tabla_pares_bandas$percentil_tierra - 50) +
  abs(tabla_pares_bandas$percentil_satelite - 50)

orden_central <- order(tabla_pares_bandas$distancia_mediana, tabla_pares_bandas$curva)
orden_RMSE <- order(-tabla_pares_bandas$RMSE_tierra, tabla_pares_bandas$curva)
cat("Pares más cercanos a las dos medianas\n")
print(head(tabla_pares_bandas[orden_central, ], 5), row.names = FALSE)
cat("Pares con mayor RMSE terrestre (mm)\n")
print(head(tabla_pares_bandas[orden_RMSE, ], 5), row.names = FALSE)


# Bandas punto a punto para los pares elegidos manualmente
curvas_bandas <- c("2610000111_1989", "5421200110_1982")
ajustes_bandas <- list(terrestre = ajustes_estaciones, satelital = ajustes_satelite)

# Multiplicador de referencia para visualizar la incertidumbre del ajuste
multiplicador_bandas <- 2

bandas_MCP <- list()
for (curva_bandas in curvas_bandas) {
  par_bandas <- match(curva_bandas, tabla_pares_bandas$curva)
  bandas_par <- list()
  for (fuente in names(ajustes_bandas)) {
    ajuste_banda <- ajustes_bandas[[fuente]][par_bandas, ]
    EE_banda <- sqrt(varianza_curvas_MCP[[fuente]][par_bandas, ])
    bandas_par[[fuente]] <- data.frame(
      observado = fuentes[[fuente]][par_bandas, ],
      ajuste = ajuste_banda,
      inferior = ajuste_banda - multiplicador_bandas * EE_banda,
      superior = ajuste_banda + multiplicador_bandas * EE_banda
    )
  }
  bandas_MCP[[curva_bandas]] <- bandas_par
  limites_bandas <- range(bandas_par$terrestre, bandas_par$satelital, na.rm = TRUE)
  par_bandas_anterior <- par(mfrow = c(1, 2))
  for (fuente in names(bandas_par)) {
    banda <- bandas_par[[fuente]]
    color <- colores_bandas[[fuente]]
    plot(tiempos, banda$observado, pch = 16, cex = 0.65, ylim = limites_bandas,
         xlab = "Péntada", ylab = "Precipitación (mm)",
         main = paste(etiquetas_bandas[[fuente]], curva_bandas))
    lines(tiempos, banda$ajuste, col = color, lwd = 2)
    lines(tiempos, banda$inferior, col = color, lty = 2)
    lines(tiempos, banda$superior, col = color, lty = 2)
    legend("topright", c("Observado", "Ajuste MCP", paste("Referencia ±", multiplicador_bandas, "EE")),
           col = c("black", color, color), pch = c(16, NA, NA),
           lty = c(NA, 1, 2), bty = "n", cex = 0.8)
  }
  par(par_bandas_anterior)
}

