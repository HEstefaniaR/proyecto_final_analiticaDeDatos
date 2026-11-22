#################### Proyecto 3 ####################
# Por: Natalia Paredes, Mariana Cifuentes, John Franco
# y Estefanía Hernández
####################################################

#################### Librerias ####################
#install.packages("fda")
#install.packages("readxl")
#install.packages("tidyverse")
#install.packages("terra")
#install.packages("geodata")
#install.packages("fda.usc")
library(fda)
library(readxl)
library(terra)
library(geodata)
library(tidyverse)

set.seed(2026)
#################### Directorio de trabajo ####################
rm(list = ls())

# Carpetas
dir.create("Tablas", showWarnings = FALSE)
dir.create("Figuras", showWarnings = FALSE)
dir.create("Curvas", showWarnings = FALSE)

guardar_tabla <- function(tabla, nombre) {
  write.csv(tabla, file.path("Tablas", nombre), row.names = FALSE, fileEncoding = "UTF-8")
}

guardar_figura <- function(nombre, ancho = 10, alto = 6) {
  dev.copy(png, file.path("Figuras", nombre), width = ancho, height = alto, units = "in", res = 150)
  dev.off()
}

#################### Carga de datos ####################
# Precipitacion diaria en tierra
p_tierra_1980_1990 <- read.csv("data/P_1980_1990.csv")
p_tierra_1991_2003 <- read.csv("data/P_1991_2003.csv")

p_tierra <- rbind(p_tierra_1980_1990, p_tierra_1991_2003)
p_tierra <- p_tierra[p_tierra$YEAR >= 1981, ]
p_tierra$fecha <- as.Date(sprintf("%04d-%02d-%02d", p_tierra$YEAR, p_tierra$MONTH, p_tierra$DAY))
p_tierra <- p_tierra[order(p_tierra$STATION_ID, p_tierra$fecha), ]

# Verificar duplicados
sum(duplicated(p_tierra[, c("STATION_ID", "fecha")]))

# Dias faltantes por año 
faltantes_por_anio <- p_tierra %>%
  group_by(STATION_ID) %>%
  complete(fecha = seq(min(fecha), max(fecha), by = "day")) %>%
  ungroup() %>%
  mutate(YEAR = as.integer(format(fecha, "%Y"))) %>%
  group_by(YEAR) %>%
  summarize(total = n(), faltantes = sum(is.na(VALUE)),
            porcentaje_faltantes = round(100 * mean(is.na(VALUE)), 1))
faltantes_por_anio
guardar_tabla(faltantes_por_anio, "01_preprocesamiento_faltantes_por_anio.csv")

# Dias faltantes por año y estacion
faltantes_por_anio_estacion <- p_tierra %>%
  group_by(STATION_ID) %>%
  complete(fecha = seq(min(fecha), max(fecha), by = "day")) %>%
  ungroup() %>%
  mutate(YEAR = as.integer(format(fecha, "%Y"))) %>%
  group_by(STATION_ID, YEAR) %>%
  summarize(total = n(), faltantes = sum(is.na(VALUE)),
            porcentaje_faltantes = round(100 * mean(is.na(VALUE)), 1),
            .groups = "drop")
guardar_tabla(faltantes_por_anio_estacion, "01_preprocesamiento_faltantes_por_anio_estacion.csv")

# Conteo Faltantes año-estacion
conteo_por_anio <- faltantes_por_anio_estacion %>%
  group_by(YEAR) %>%
  summarize(estaciones_20_o_mas = sum(porcentaje_faltantes >= 20),
            estaciones_menos_20 = sum(porcentaje_faltantes < 20),
            .groups = "drop")
conteo_por_anio
guardar_tabla(conteo_por_anio, "01_preprocesamiento_conteo_20pct_por_anio.csv")

contar <- function(datos, paso) {
  validos <- datos[!is.na(datos$VALUE), ]
  data.frame(paso = paso, registros = nrow(validos),
             estaciones = length(unique(validos$STATION_ID)),
             pares_estacion_anio = nrow(unique(validos[, c("STATION_ID", "YEAR")])))
}

trazabilidad_pasos <- contar(p_tierra, "Datos 1981-2003")

# Coordenadas y altitud de las estaciones
coord <- read.csv("data/CVC_Wheather_Stations_Long_Lat.csv", fileEncoding = "latin1")
coord$STATION_ID <- as.character(coord$CODIGO)

alt <- read_excel("data/CatalogoEstacionesElevacion.xlsx")
alt$STATION_ID <- as.character(alt$codigo)

info_estaciones <- merge(data.frame(STATION_ID = coord$STATION_ID, longitud = coord$Lon, latitud = coord$lat),
                         data.frame(STATION_ID = alt$STATION_ID, altitud = alt$alt),
                         by = "STATION_ID", all = TRUE)

# Imputar altitud faltante con satelite
elevacion_col <- elevation_30s(country = "COL", path = tempdir())
con_coordenadas <- !is.na(info_estaciones$longitud) & !is.na(info_estaciones$latitud)
info_estaciones$altitud_DEM <- NA
info_estaciones$altitud_DEM[con_coordenadas] <-
  terra::extract(elevacion_col, as.matrix(info_estaciones[con_coordenadas, c("longitud", "latitud")]))[, 1]

# Comparacion con las estaciones que tienen altitud
ambas <- !is.na(info_estaciones$altitud) & !is.na(info_estaciones$altitud_DEM)
comparacion_altitud <- data.frame(
  estaciones = sum(ambas),
  diferencia_media_m = mean(info_estaciones$altitud_DEM[ambas] - info_estaciones$altitud[ambas]),
  diferencia_absoluta_media_m = mean(abs(info_estaciones$altitud_DEM[ambas] - info_estaciones$altitud[ambas])),
  diferencia_absoluta_mediana_m = median(abs(info_estaciones$altitud_DEM[ambas] - info_estaciones$altitud[ambas])),
  correlacion = cor(info_estaciones$altitud_DEM[ambas], info_estaciones$altitud[ambas]))
comparacion_altitud

plot(info_estaciones$altitud[ambas], info_estaciones$altitud_DEM[ambas], pch = 16, cex = 0.6,
     xlab = "Altitud datos (m)", ylab = "Altitud satélite (m)", main = "Altitud: de los datos frente a satélite")
abline(0, 1, lty = 2, col = "red")
guardar_figura("01_altitud_imputada_comparacion.png", 6, 6)

# Imputar altitud si es menor a 2000 m
usar_dem <- is.na(info_estaciones$altitud) & !is.na(info_estaciones$altitud_DEM) & info_estaciones$altitud_DEM < 2000
info_estaciones$altitud_completa <- ifelse(usar_dem, info_estaciones$altitud_DEM, info_estaciones$altitud)
info_estaciones$fuente_altitud <- ifelse(!is.na(info_estaciones$altitud), "Datos", ifelse(usar_dem, "Imputada", NA))
table(info_estaciones$fuente_altitud, useNA = "ifany")
guardar_tabla(comparacion_altitud, "01_altitud_imputada_comparacion.csv")

# Borde del Valle
departamentos_col <- gadm(country = "COL", level = 1, path = tempdir())
valle <- departamentos_col[departamentos_col$NAME_1 == "Valle del Cauca", ]



#################### Particion por estaciones ####################
# Entrenamiento 1981-1997 y prueba 1998-2003
# se asegura que las estaciones del mismo pixel queden en la misma particion
estaciones_particion <- info_estaciones[info_estaciones$STATION_ID %in% p_tierra$STATION_ID &
  !is.na(info_estaciones$longitud) & !is.na(info_estaciones$latitud) &
  !is.na(info_estaciones$altitud_completa), c("STATION_ID", "longitud", "latitud", "altitud_completa")]
estaciones_particion$STATION_ID <- as.character(estaciones_particion$STATION_ID)
raster_pixel <- rast(list.files("data/ImagenesSatelitePentadal", pattern = "ppt.*\\.tif$", full.names = TRUE)[1])
estaciones_particion$pixel <- cellFromXY(raster_pixel, as.matrix(estaciones_particion[, c("longitud", "latitud")]))

estaciones_particion$bloque <- ifelse(is.na(estaciones_particion$pixel), paste0("estacion_", estaciones_particion$STATION_ID),
                                      paste0("pixel_", estaciones_particion$pixel))
cortes_altitud <- quantile(estaciones_particion$altitud_completa, c(0, 1/3, 2/3, 1))
estaciones_particion$tercil_altitud <- cut(estaciones_particion$altitud_completa, cortes_altitud,
  labels = c("Baja", "Media", "Alta"), include.lowest = TRUE)

bloques <- estaciones_particion %>%
  group_by(bloque) %>%
  summarize(altitud = median(altitud_completa),
            con_prueba = any(STATION_ID %in% as.character(p_tierra$STATION_ID[p_tierra$YEAR >= 1998])),
            con_entrenamiento = any(STATION_ID %in% as.character(p_tierra$STATION_ID[p_tierra$YEAR <= 1997])),
            .groups = "drop")
bloques$tercil_altitud <- cut(bloques$altitud, cortes_altitud, labels = c("Baja", "Media", "Alta"), include.lowest = TRUE)

set.seed(2026)
bloques$particion <- "entrenamiento"
for (tercil in levels(bloques$tercil_altitud)) {
  en_tercil <- bloques$tercil_altitud == tercil
  candidatas <- which(en_tercil & bloques$con_prueba)
  n_reservadas <- min(length(candidatas), round(0.15 * sum(en_tercil)))
  if (n_reservadas > 0) bloques$particion[sample(candidatas, n_reservadas)] <- "prueba"
}
bloques$grupo_vc <- NA
for (tercil in levels(bloques$tercil_altitud)) {
  filas <- which(bloques$tercil_altitud == tercil & bloques$particion == "entrenamiento" & bloques$con_entrenamiento)
  bloques$grupo_vc[filas] <- sample(rep(1:5, length.out = length(filas)))
}
table(bloques$tercil_altitud, bloques$particion)
table(bloques$tercil_altitud, bloques$grupo_vc)

estaciones_particion$particion <- bloques$particion[match(estaciones_particion$bloque, bloques$bloque)]
con_entrenamiento <- estaciones_particion$STATION_ID %in% as.character(p_tierra$STATION_ID[p_tierra$YEAR <= 1997])
estaciones_particion$grupo_vc <- ifelse(con_entrenamiento,
                                        bloques$grupo_vc[match(estaciones_particion$bloque, bloques$bloque)], NA)
table(estaciones_particion$tercil_altitud, estaciones_particion$particion)

con_pixel <- estaciones_particion[!is.na(estaciones_particion$pixel), ]
stopifnot(all(tapply(con_pixel$particion, con_pixel$pixel, function(x) length(unique(x))) == 1),
          all(tapply(con_pixel$grupo_vc, con_pixel$pixel, function(x) length(unique(na.omit(x))) <= 1)))
grupo_vc <- function(estacion) {
  estaciones_particion$grupo_vc[match(as.character(estacion), estaciones_particion$STATION_ID)]
}
guardar_tabla(estaciones_particion[order(estaciones_particion$particion, estaciones_particion$STATION_ID), ],
  "01_particion_estaciones.csv")
particion_curva <- function(estacion, anio) {
  grupo <- estaciones_particion$particion[match(as.character(estacion), estaciones_particion$STATION_ID)]
  grupo[is.na(grupo)] <- "sin_uso"
  ifelse(grupo == "entrenamiento" & anio >= 1981 & anio <= 1997, "entrenamiento",
    ifelse(grupo == "prueba" & anio >= 1998 & anio <= 2003, "prueba", "sin_uso"))
}
estaciones_train <- estaciones_particion$STATION_ID[estaciones_particion$particion == "entrenamiento"]
estaciones_test <- estaciones_particion$STATION_ID[estaciones_particion$particion == "prueba"]
stopifnot(length(intersect(estaciones_train, estaciones_test)) == 0)
stopifnot(length(estaciones_train) >= 5, length(estaciones_test) > 0)
es_entrenamiento <- function(estacion) {
  grupo <- estaciones_particion$particion[match(as.character(estacion), estaciones_particion$STATION_ID)]
  !is.na(grupo) & grupo == "entrenamiento"
}

#################### Preprocesamiento de las estaciones ####################
# Las reglas se revisan con las estaciones de entrenamiento (1981-1997) y se aplican a todos los datos
p_revision <- p_tierra[es_entrenamiento(p_tierra$STATION_ID) & p_tierra$YEAR <= 1997, ]

# Valores negativos 
sum(p_tierra$VALUE < 0, na.rm = TRUE)

summary(p_revision)
quantile(p_revision$VALUE, probs = c(0.5, 0.9, 0.95, 0.99, 0.999, 0.99999), na.rm = TRUE)
tapply(p_revision$VALUE, p_revision$YEAR, summary)

# Valores mayores a 300 mm
mayores_300 <- p_revision[which(p_revision$VALUE > 300), c("STATION_ID", "YEAR", "MONTH", "DAY", "VALUE")]
mayores_300 <- mayores_300[order(-mayores_300$YEAR), ]

n_mayores <- nrow(mayores_300)
n_total <- sum(!is.na(p_revision$VALUE))
round(100 * n_mayores / n_total, 4)

# Histograma sin los valores extremos
hist(p_revision$VALUE[p_revision$VALUE < 130],
     main = "Precipitación diaria", xlab = "Precipitación (mm)", ylab = "Frecuencia")
guardar_figura("01_preprocesamiento_histograma_diaria.png", 7, 5)

# Outliers en cada estacion y año
casos <- unique(mayores_300[, c("STATION_ID", "YEAR")])
casos <- casos[order(casos$YEAR, casos$STATION_ID), ]

# Una imagen cada 6 casos
png("Figuras/01_preprocesamiento_casos_mayores_300mm_%02d.png", width = 12, height = 7, units = "in", res = 150)
par(mfrow = c(2, 3))
for (i in 1:nrow(casos)) {
  d <- p_revision[p_revision$STATION_ID == casos$STATION_ID[i] & p_revision$YEAR == casos$YEAR[i], ]
  plot(d$fecha, d$VALUE, type = "b", pch = 16, cex = 0.4, lwd = 0.5, col = "steelblue",
       ylim = c(0, min(max(d$VALUE), 1000)),
       xlab = "Fecha", ylab = "Precipitación (mm)",
       main = paste("Estación", casos$STATION_ID[i], "-", casos$YEAR[i]))
  points(d$fecha[d$VALUE > 200], d$VALUE[d$VALUE > 200], col = "red", pch = 19)
}
dev.off()

# Errores aislados: mas de 400 mm, 5 veces sus vecinos y no a fin de mes
p_tierra$fin_mes <- format(p_tierra$fecha + 1, "%d") == "01"

n <- nrow(p_tierra)
anterior  <- c(NA, p_tierra$VALUE[-n])
siguiente <- c(p_tierra$VALUE[-1], NA)
vecino_max <- pmax(anterior, siguiente, na.rm = TRUE)

error <- which(p_tierra$VALUE > 400 & p_tierra$VALUE > 5 * vecino_max & !p_tierra$fin_mes)

errores_aislados <- data.frame(p_tierra[error, c("STATION_ID", "YEAR", "MONTH", "DAY")],
                               anterior = anterior[error], VALUE = p_tierra$VALUE[error], siguiente = siguiente[error])
errores_aislados
guardar_tabla(errores_aislados, "01_preprocesamiento_errores_aislados.csv")

# Grafico valores altos de precipitacion
png("Figuras/01_errores_aislados_anio_completo.png", width = 10, height = 3 * nrow(errores_aislados),
    units = "in", res = 150)
par(mfrow = c(nrow(errores_aislados), 1), mar = c(4, 4, 3, 1))
for (i in seq_len(nrow(errores_aislados))) {
  s <- errores_aislados$STATION_ID[i]
  f <- p_tierra$fecha[error[i]]
  anio <- errores_aislados$YEAR[i]
  datos <- p_tierra[p_tierra$STATION_ID == s & p_tierra$YEAR == anio, ]
  valor_error <- errores_aislados$VALUE[i]
  
  y_resto <- datos$VALUE[datos$fecha != f]
  ymax <- max(c(y_resto, 50), na.rm = TRUE) * 1.2
  
  plot(datos$fecha, ifelse(datos$fecha == f, NA, datos$VALUE), type = "h",
       ylim = c(0, ymax), xlab = "Fecha", ylab = "Precipitacion (mm)",
       main = paste0("Estacion ", s, " - año ", anio))
  points(f, ymax, pch = 17, col = "red", cex = 1.5)
  text(f, ymax, labels = paste(format(f, "%d-%m"), ":", valor_error, "mm"), pos = 1, col = "red", cex = 0.8)
}
dev.off()
par(mfrow = c(1, 1))
guardar_figura("01_errores_aislados_anio_completo.png", 10, 3 * nrow(errores_aislados))

p_tierra$VALUE[error] <- NA

trazabilidad_pasos <- rbind(trazabilidad_pasos, contar(p_tierra, "Errores aislados"))

# Acumulados mensuales
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

# Candidatos: ultimo dia con 100 mm o mas y al menos el 80 % del mes sin dato
candidatos_M <- resumen_mensual[which(resumen_mensual$valor_fin_mes >= 100 &
                                        resumen_mensual$proporcion_fin_mes >= 0.80), ]
nrow(candidatos_M)
print(candidatos_M, n = 40)
guardar_tabla(candidatos_M, "01_preprocesamiento_acumulados_mensuales.csv")

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

# El fin de mes de los candidatos pasa a NA
clave_candidatos <- paste(candidatos_M$STATION_ID, candidatos_M$YEAR, candidatos_M$MONTH)
acumulados <- which(p_m$fin_mes & !is.na(p_m$VALUE) &
                      paste(p_m$STATION_ID, p_m$YEAR, p_m$MONTH) %in% clave_candidatos)
length(acumulados)

p_m_depurado <- p_m
p_m_depurado$VALUE[acumulados] <- NA

trazabilidad_pasos <- rbind(trazabilidad_pasos, contar(p_m_depurado, "Acumulados mensuales"))

# Trazabilidad
trazabilidad <- rbind(data.frame(p_m[acumulados, c("STATION_ID", "fecha", "VALUE")],
                                 motivo = "Precipitacion final >= 100 mm y >= 80% del total mensual"))
trazabilidad
guardar_tabla(trazabilidad, "01_preprocesamiento_registros_eliminados.csv")

# Dias faltantes y años incompletos
# Racha mas larga de dias faltantes
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
fechas <- data.frame(fecha = seq(as.Date("1981-01-01"), as.Date("2003-12-31"), by = "day"))
fechas$YEAR <- as.numeric(format(fechas$fecha, "%Y"))

calendario <- merge(unique(p_m_depurado[, c("STATION_ID", "YEAR")]), fechas, by = "YEAR")
calendario <- merge(calendario, p_m_depurado[, c("STATION_ID", "fecha", "VALUE")],
                    by = c("STATION_ID", "fecha"), all.x = TRUE)
calendario <- calendario[order(calendario$STATION_ID, calendario$fecha), ]


# Imputacion por IDW con las K estaciones mas cercanas
# Estacion-año con menos del 20 % de dias faltantes; las vecinas son estaciones de entrenamiento con dato ese dia
# y de la misma zona: valle hasta 1239 m y montaña por encima (Russi y Larrahondo, 2017)
LIMITE_VALLE <- 1239
info_imp <- info_estaciones[!is.na(info_estaciones$longitud) & !is.na(info_estaciones$latitud) &
                              !is.na(info_estaciones$altitud_completa), ]
info_imp$zona_imputacion <- ifelse(info_imp$altitud_completa <= LIMITE_VALLE, "Valle", "Montaña")
puntos_imp <- vect(info_imp[, c("STATION_ID", "longitud", "latitud")], geom = c("longitud", "latitud"), crs = "EPSG:4326")
paisajes_imp <- vect("data/Paisajes_IGAC_hidrografia_IDEAM.gpkg", layer = "paisajes_IGAC")
hidrografia_imp <- vect("data/Paisajes_IGAC_hidrografia_IDEAM.gpkg", layer = "hidrografia_IDEAM")
paisaje_punto_imp <- terra::extract(paisajes_imp, project(puntos_imp, crs(paisajes_imp)))
cuenca_punto_imp <- terra::extract(hidrografia_imp, project(puntos_imp, crs(hidrografia_imp)))
stopifnot(identical(as.integer(paisaje_punto_imp$id.y), seq_len(nrow(info_imp))),
          identical(as.integer(cuenca_punto_imp$id.y), seq_len(nrow(info_imp))))
info_imp$zona_imputacion[grepl("^Planicie", paisaje_punto_imp$PAISAJE) & cuenca_punto_imp$nom_ah %in% "Pacifico"] <- "Costa"
table(info_imp$zona_imputacion)
lat0 <- mean(info_imp$latitud)
info_imp$x <- info_imp$longitud * 111.32 * cos(lat0 * pi / 180)
info_imp$y <- info_imp$latitud * 110.57

faltante <- calendario %>%
  group_by(STATION_ID, YEAR) %>%
  summarize(pct_faltante = 100 * mean(is.na(VALUE)), .groups = "drop")
imputables <- faltante[faltante$pct_faltante > 0 & faltante$pct_faltante < 20, ]

ancho <- calendario %>%
  dplyr::select(STATION_ID, fecha, VALUE) %>%
  pivot_wider(names_from = STATION_ID, values_from = VALUE) %>%
  arrange(fecha)
donantes <- intersect(as.character(info_imp$STATION_ID[es_entrenamiento(info_imp$STATION_ID) &
                                                         info_imp$zona_imputacion != "Costa"]), names(ancho))
zona_donante <- info_imp$zona_imputacion[match(donantes, info_imp$STATION_ID)]

# Estimacion IDW de todos los dias de una estacion para cada K de la lista (peso 1 / distancia)
estimar_idw <- function(s, K_lista) {
  i_s <- match(s, info_imp$STATION_ID)
  candidatas <- setdiff(donantes[zona_donante == info_imp$zona_imputacion[i_s]], s)
  j <- match(candidatas, info_imp$STATION_ID)
  d <- pmax(sqrt((info_imp$x[j] - info_imp$x[i_s])^2 + (info_imp$y[j] - info_imp$y[i_s])^2), 0.1)
  orden <- order(d)
  candidatas <- candidatas[orden]
  d <- d[orden]
  suma <- peso <- n_vecinos <- numeric(nrow(ancho))
  estimado <- matrix(NA_real_, nrow(ancho), length(K_lista))
  for (k in seq_along(candidatas)) {
    v <- ancho[[candidatas[k]]]
    usar <- !is.na(v) & n_vecinos < max(K_lista)
    suma[usar] <- suma[usar] + v[usar] / d[k]
    peso[usar] <- peso[usar] + 1 / d[k]
    n_vecinos[usar] <- n_vecinos[usar] + 1
    for (m in seq_along(K_lista)) {
      llega <- usar & n_vecinos == K_lista[m]
      estimado[llega, m] <- suma[llega] / peso[llega]
    }
    if (all(n_vecinos >= max(K_lista))) break
  }
  estimado
}

# Eleccion de K por validacion cruzada
# El K elegido se aplica a las de prueba
K_lista <- c(1, 3, 5, 7, 10, 15, 20, 25, 30)
dias_entrenamiento <- as.numeric(format(ancho$fecha, "%Y")) <= 1997
errores_K <- data.frame()
for (s in donantes) {
  observado <- ancho[[s]]
  estimado <- estimar_idw(s, K_lista)
  usar <- dias_entrenamiento & !is.na(observado) & rowSums(is.na(estimado)) == 0
  errores_K <- rbind(errores_K, data.frame(
    zona = info_imp$zona_imputacion[match(s, info_imp$STATION_ID)], K = K_lista,
    suma_error2 = colSums((estimado[usar, , drop = FALSE] - observado[usar])^2, na.rm = TRUE),
    n = colSums(!is.na(estimado[usar, , drop = FALSE]))))
}
seleccion_K <- errores_K %>%
  group_by(zona, K) %>%
  summarize(RECM_mm = sqrt(sum(suma_error2) / sum(n)), dias = sum(n), .groups = "drop")
seleccion_K
guardar_tabla(seleccion_K, "01_preprocesamiento_imputacion_seleccion_K.csv")
K_zona <- seleccion_K %>% group_by(zona) %>% slice_min(RECM_mm, n = 1)
K_zona

plot(K_lista, seleccion_K$RECM_mm[seleccion_K$zona == "Valle"], type = "b", pch = 19,
     ylim = range(seleccion_K$RECM_mm), xlab = "K (estaciones vecinas)", ylab = "RECM de validacion cruzada (mm)",
     main = "Eleccion de K para la imputacion IDW")
lines(K_lista, seleccion_K$RECM_mm[seleccion_K$zona == "Montaña"], type = "b", pch = 17, col = "red")
legend("topright", c("Valle", "Montaña"), col = c("black", "red"), pch = c(19, 17), bty = "n")
guardar_figura("01_preprocesamiento_imputacion_seleccion_K.png", 7, 5)

# Imputacion con el K elegido en cada zona
calendario$imputado <- FALSE
sin_costa <- as.character(info_imp$STATION_ID[info_imp$zona_imputacion != "Costa"])
for (s in intersect(as.character(unique(imputables$STATION_ID)), sin_costa)) {
  K_s <- K_zona$K[K_zona$zona == info_imp$zona_imputacion[match(s, info_imp$STATION_ID)]]
  estimado <- estimar_idw(s, K_s)[, 1]
  anios <- imputables$YEAR[as.character(imputables$STATION_ID) == s]
  filas <- which(as.character(calendario$STATION_ID) == s & is.na(calendario$VALUE) & calendario$YEAR %in% anios)
  valor <- estimado[match(calendario$fecha[filas], ancho$fecha)]
  calendario$VALUE[filas] <- valor
  calendario$imputado[filas] <- !is.na(valor)
}
calendario$zona_imputacion <- info_imp$zona_imputacion[match(as.character(calendario$STATION_ID), info_imp$STATION_ID)]
en_imputables <- paste(calendario$STATION_ID, calendario$YEAR) %in% paste(imputables$STATION_ID, imputables$YEAR)
imputacion <- calendario[en_imputables & !is.na(calendario$zona_imputacion), ] %>%
  group_by(zona = zona_imputacion) %>%
  summarize(K = K_zona$K[match(first(zona), K_zona$zona)],
            estaciones = n_distinct(STATION_ID),
            estacion_anio = n_distinct(paste(STATION_ID, YEAR)),
            dias_imputados = sum(imputado),
            dias_sin_vecinas = sum(is.na(VALUE)), .groups = "drop")
imputacion
guardar_tabla(imputacion, "01_preprocesamiento_imputacion.csv")
trazabilidad_pasos <- rbind(trazabilidad_pasos, contar(calendario, "Imputacion IDW"))

# Cobertura, hueco mas largo y faltantes en los bordes
cobertura <- calendario %>%
  group_by(STATION_ID, YEAR) %>%
  summarize(dias = n(),
            dias_validos = sum(!is.na(VALUE)),
            hueco = racha_max(VALUE),
            borde = borde_max(VALUE),
            .groups = "drop")

cobertura$cobertura <- cobertura$dias_validos / cobertura$dias
summary(cobertura)

# Umbrales de cobertura (entrenamiento)
train <- cobertura[cobertura$YEAR <= 1997 & es_entrenamiento(cobertura$STATION_ID), ]
umbrales_cobertura <- data.frame()
for (u in c(0.85, 0.90, 0.95)) {
  umbrales_cobertura <- rbind(umbrales_cobertura,
                              data.frame(umbral = u, solo_cobertura = sum(train$cobertura >= u),
                                         con_hueco_y_borde_7 = sum(train$cobertura >= u & train$hueco <= 7 &
                                                                     train$borde <= 7)))
}
umbrales_cobertura
guardar_tabla(umbrales_cobertura, "01_preprocesamiento_umbrales_cobertura.csv")

# Criterio: 95 % de cobertura y huecos de maximo 7 dias
cobertura$apto <- cobertura$cobertura >= 0.95 & cobertura$hueco <= 7 & cobertura$borde <= 7

trazabilidad_pasos <- rbind(trazabilidad_pasos, 
                            contar(p_m_depurado[paste(p_m_depurado$STATION_ID, 
                                                      p_m_depurado$YEAR) %in% 
                                                  paste(cobertura$STATION_ID, cobertura$YEAR)[cobertura$apto], ], "Cobertura minima"))

cat("Pares estacion-año antes:", nrow(cobertura), "| aptos:", sum(cobertura$apto), "\n")
cat("Estaciones antes:", length(unique(cobertura$STATION_ID)),
    "| despues:", length(unique(cobertura$STATION_ID[cobertura$apto])), "\n")

# Pentadas de las estaciones
# Datos diarios de los pares aptos
pares_aptos <- paste(cobertura$STATION_ID, cobertura$YEAR)[cobertura$apto]
datos_tierra <- calendario[paste(calendario$STATION_ID, calendario$YEAR) %in% pares_aptos, ]

datos_tierra$MONTH <- as.numeric(format(datos_tierra$fecha, "%m"))
datos_tierra$DAY <- as.numeric(format(datos_tierra$fecha, "%d"))
datos_tierra$pentada <- (datos_tierra$MONTH - 1) * 6 + pmin(6, (datos_tierra$DAY - 1) %/% 5 + 1)
datos_tierra$STATION_ID <- as.character(datos_tierra$STATION_ID)
datos_tierra$curva <- paste(datos_tierra$STATION_ID, datos_tierra$YEAR, sep = "_")

# Lluvia por pentada
tierra_pentadal <- datos_tierra %>%
  group_by(curva, STATION_ID, YEAR, pentada) %>%
  summarize(lluvia = sum(VALUE), .groups = "drop")

# Matriz terrestre (pentadas por curvas)
Y_terrestre_total <- tapply(tierra_pentadal$lluvia, list(tierra_pentadal$pentada, tierra_pentadal$curva), sum)
dim(Y_terrestre_total)


#################### Fases ENSO (ONI) ####################
tabla_oni <- read_excel("data/ONI_1981_2003.xlsx")
names(tabla_oni)[-1] <- c("DEF", "EFM", "FMA", "MAM", "AMJ", "MJJ",
                          "JJA", "JAS", "ASO", "SON", "OND", "NDE")

# Ventanas ordenadas por su mes central (DEF es enero)
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

# Rachas continuas entre años
rachas_oni <- rle(oni_temporal$fase_umbral)
oni_temporal$racha <- rep(seq_along(rachas_oni$lengths), rachas_oni$lengths)
oni_temporal$ventanas_racha <- rep(rachas_oni$lengths, rachas_oni$lengths)

# Criterio NOAA: 5 ventanas seguidas del mismo signo
oni_temporal <- oni_temporal %>%
  mutate(fase_enso = ifelse(ventanas_racha >= 5, fase_umbral, "Neutral"))

# Episodios ENSO
episodios_enso <- oni_temporal %>%
  filter(fase_enso != "Neutral") %>%
  group_by(racha, fase_enso) %>%
  summarise(inicio = paste(first(ventana), first(YEAR)),
            fin = paste(last(ventana), last(YEAR)),
            ventanas = n(),
            .groups = "drop") %>%
  mutate(episodio = row_number()) %>%
  select(episodio, fase_enso, inicio, fin, ventanas)

# Fases de cada año
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

# ONI anual para el modelo
oni_anual <- oni_temporal %>%
  group_by(YEAR) %>%
  summarise(oni_promedio = mean(ONI),
            oni_max_intensidad = ONI[which.max(abs(ONI))],
            .groups = "drop")
print(oni_anual, n=Inf)

# Fase por año y mes
enso_evaluacion <- oni_temporal %>%
  select(YEAR, mes_central, ventana, ONI, fase_enso)

head(enso_evaluacion, n=12)

write.csv(oni_anual, "data/oni_anual.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(enso_evaluacion, "data/enso_evaluacion.csv", row.names = FALSE, fileEncoding = "UTF-8")

print(episodios_enso, n = Inf)
print(enso_anual, n = Inf, width = Inf)

guardar_tabla(episodios_enso, "02_ONI_episodios.csv")
guardar_tabla(enso_anual, "02_ONI_fases_por_anio.csv")
guardar_tabla(oni_anual, "02_ONI_resumen_anual.csv")


#################### Emparejamiento estacion-pixel ####################
# Satelite y pixel de cada estacion
# Imagenes pentadales de 1981 a 2003
archivos <- list.files("data/ImagenesSatelitePentadal", pattern = "ppt.*\\.tif$", full.names = TRUE)
anio_archivo <- as.numeric(substr(basename(archivos), 4, 7))
archivos <- archivos[anio_archivo <= 2003]
# 72 imagenes por año
table(substr(basename(archivos), 4, 7))

p_sat <- rast(archivos)

# Estaciones con coordenadas y el pixel que les corresponde
curvas <- unique(tierra_pentadal[, c("curva", "STATION_ID", "YEAR")])
curvas <- merge(curvas, info_estaciones, by = "STATION_ID", all.x = TRUE)
curvas <- curvas[order(curvas$curva), ]

estaciones <- unique(curvas[!is.na(curvas$longitud), c("STATION_ID", "longitud", "latitud")])
puntos <- vect(estaciones, geom = c("longitud", "latitud"), crs = "EPSG:4326")
estaciones$pixel <- cellFromXY(p_sat, crds(puntos))

# Lluvia satelital en el pixel de cada estacion 
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
guardar_tabla(as.data.frame(table(estado = curvas$estado)), "03_emparejamiento_estado_curvas.csv")

table(estaciones$pixel)[table(estaciones$pixel) > 1]
en_raster <- estaciones[!is.na(estaciones$pixel), ]
pixeles_compartidos <- en_raster[en_raster$pixel %in% en_raster$pixel[duplicated(en_raster$pixel)], ]
pixeles_compartidos$particion <- estaciones_particion$particion[match(pixeles_compartidos$STATION_ID,
                                                                      estaciones_particion$STATION_ID)]
pixeles_compartidos$grupo_vc <- grupo_vc(pixeles_compartidos$STATION_ID)
pixeles_compartidos <- pixeles_compartidos[order(pixeles_compartidos$pixel), ]
guardar_tabla(pixeles_compartidos, "03_emparejamiento_pixeles_compartidos.csv")

# Curvas emparejadas
metadatos <- curvas[curvas$estado == "Emparejado", ]

# Trazabilidad
trazabilidad_pasos <- rbind(trazabilidad_pasos,
                            contar(datos_tierra[datos_tierra$curva %in% metadatos$curva, ], "Emparejamiento estacion-pixel"))
trazabilidad_pasos
guardar_tabla(trazabilidad_pasos, "01_preprocesamiento_trazabilidad.csv")

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

# Matrices y particion
tiempos <- 1:72

# Filas: curvas; columnas: pentadas
Y_matriz <- t(Y_est)
X_matriz <- t(X_sat)
colnames(Y_matriz) <- paste0("P", sprintf("%02d", tiempos))
colnames(X_matriz) <- paste0("P", sprintf("%02d", tiempos))

# Entrenamiento 1981-1997 y prueba 1998-2003
particion_metadatos <- particion_curva(metadatos$STATION_ID, metadatos$YEAR)
table(particion_metadatos)

Y_train <- Y_matriz[particion_metadatos == "entrenamiento", , drop = FALSE]
X_train <- X_matriz[particion_metadatos == "entrenamiento", , drop = FALSE]
Y_test <- Y_matriz[particion_metadatos == "prueba", , drop = FALSE]
X_test <- X_matriz[particion_metadatos == "prueba", , drop = FALSE]

# Verificar periodos
stopifnot(all(metadatos$YEAR[particion_metadatos == "entrenamiento"] >= 1981 &
              metadatos$YEAR[particion_metadatos == "entrenamiento"] <= 1997))
stopifnot(all(metadatos$YEAR[particion_metadatos == "prueba"] >= 1998 &
              metadatos$YEAR[particion_metadatos == "prueba"] <= 2003))

# Verificar que no se repitan estaciones
estaciones_train <- unique(as.character(metadatos$STATION_ID[particion_metadatos == "entrenamiento"]))
estaciones_test <- unique(as.character(metadatos$STATION_ID[particion_metadatos == "prueba"]))
stopifnot(length(intersect(estaciones_train, estaciones_test)) == 0)

# Resumen
cat("Curvas totales:", nrow(Y_matriz), "\n")
cat("Curvas entrenamiento:", nrow(Y_train), "\n")
cat("Curvas prueba:", nrow(Y_test), "\n")
cat("Estaciones entrenamiento:", length(estaciones_train), "\n")
cat("Estaciones prueba:", length(estaciones_test), "\n")
cat("Estaciones compartidas:", length(intersect(estaciones_train, estaciones_test)), "\n")

# Guardar las matrices
saveRDS(list(Y_matriz = Y_matriz, X_matriz = X_matriz, metadatos = metadatos,
             Y_train = Y_train, X_train = X_train, Y_test = Y_test, X_test = X_test),
        "Curvas/matrices_funcionales.rds")

# Pentadas faltantes por curva
faltantes_curva <- rowSums(is.na(Y_matriz))
table(faltantes_curva)
faltantes_curva[faltantes_curva > 0]

# Matrices emparejadas de 1981 a 2003 (la particion separa entrenamiento y prueba)
metadatos_completos <- metadatos
Y_completa <- Y_matriz
X_completa <- X_matriz


# Verificacion del emparejamiento estacion-pixel
par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))

# Sistema de referencia y cobertura
crs(p_sat, describe = TRUE)[, c("name", "code")]
same.crs(puntos, p_sat)
ext(p_sat)
res(p_sat)
# 23 años por 72 pentadas
nlyr(p_sat)
puntos_r <- project(puntos, crs(p_sat))
# Pixeles que cambian al reproyectar
sum(cellFromXY(p_sat, crds(puntos_r)) != estaciones$pixel, na.rm = TRUE)
range(estaciones$longitud)
range(estaciones$latitud)

# Estaciones fuera del raster
fuera <- unique(curvas[curvas$estado == "Fuera del raster", c("STATION_ID", "longitud", "latitud")])
nrow(fuera)
head(fuera)
plot(p_sat[[1]], main = "Estaciones sobre el raster")
lines(project(valle, crs(p_sat)), lwd = 2)
points(puntos_r, pch = 16, cex = 0.6)
points(vect(fuera, geom = c("longitud", "latitud"), crs = "EPSG:4326"), col = "red", pch = 17)
guardar_figura("03_emparejamiento_estaciones_raster.png", 7, 7)

# Estaciones sin coordenadas
sin <- unique(curvas$STATION_ID[curvas$estado == "Sin coordenadas"])
length(sin)
head(sin)
head(coord$STATION_ID)
head(alt$STATION_ID)
# Longitud de los codigos de estacion
table(nchar(sin))
table(nchar(coord$STATION_ID))
c(en_coord = sum(sin %in% coord$STATION_ID), en_alt = sum(sin %in% alt$STATION_ID))
sum(is.na(coord$Lon) | is.na(coord$lat))

# Orden cronologico de las imagenes
nm <- basename(archivos)
anio <- as.integer(substr(nm, 4, 7))
pent <- as.integer(substr(nm, 8, 9))

# 72 pentadas por año y en orden
verificacion <- tapply(pent, anio, function(p) identical(as.integer(p), 1:72))
verificacion
all(verificacion)
table(anio)

# Correlacion con rezagos: el maximo debe estar en 0
lags <- -3:3
cor_lag <- sapply(lags, function(l) {
  a <- if (l >= 0) 1:(72 - l) else (1 - l):72
  median(suppressWarnings(sapply(1:nrow(Y_matriz), function(i)
    cor(Y_matriz[i, a], X_matriz[i, a + l], use = "complete.obs"))), na.rm = TRUE)
})
correlacion_rezagos <- data.frame(rezago = lags, correlacion_mediana = round(cor_lag, 3))
correlacion_rezagos
guardar_tabla(correlacion_rezagos, "03_emparejamiento_correlacion_rezagos.csv")

# Escala de las dos fuentes y pixeles compartidos
c(media_estaciones = mean(Y_matriz, na.rm = TRUE), media_satelite = mean(X_matriz, na.rm = TRUE))
c(minimo_satelite = min(X_matriz, na.rm = TRUE), maximo_satelite = max(X_matriz, na.rm = TRUE))
# Estaciones por pixel
table(table(estaciones$pixel))


#################### Descriptivo de los datos emparejados ####################
# Pentadas de 1981 a 2003 ya limpias e imputadas, con todos los datos
stopifnot(identical(rownames(Y_matriz), metadatos$curva), identical(rownames(X_matriz), metadatos$curva))
largo <- data.frame(curva = rep(metadatos$curva, times = 72),
                    STATION_ID = rep(metadatos$STATION_ID, times = 72),
                    YEAR = rep(metadatos$YEAR, times = 72),
                    altitud = rep(metadatos$altitud_completa, times = 72),
                    pentada = rep(tiempos, each = nrow(metadatos)),
                    estacion = as.vector(Y_matriz),
                    satelite = as.vector(X_matriz))

# Dias faltantes de las estaciones en las curvas emparejadas
dias_emparejados <- datos_tierra[datos_tierra$curva %in% metadatos$curva, ]
faltantes_dias <- rbind(
  data.frame(YEAR = "Todos", pct_dias_faltantes = 100 * mean(is.na(dias_emparejados$VALUE))),
  dias_emparejados %>%
    group_by(YEAR = as.character(YEAR)) %>%
    summarize(pct_dias_faltantes = 100 * mean(is.na(VALUE)), .groups = "drop"))

# Estadisticas descriptivas por pentada
describir <- function(x) {
  data.frame(pentadas = length(x),
             pct_faltantes = 100 * mean(is.na(x)),
             media = mean(x, na.rm = TRUE),
             DE = sd(x, na.rm = TRUE),
             CV_pct = 100 * sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE),
             minimo = min(x, na.rm = TRUE),
             Q1 = unname(quantile(x, 0.25, na.rm = TRUE)),
             mediana = median(x, na.rm = TRUE),
             Q3 = unname(quantile(x, 0.75, na.rm = TRUE)),
             maximo = max(x, na.rm = TRUE))
}

# Todos los años juntos y cada año
descriptivo <- data.frame()
for (anio in c("Todos", sort(unique(largo$YEAR)))) {
  filas <- if (anio == "Todos") rep(TRUE, nrow(largo)) else largo$YEAR == anio
  for (fuente in c("estacion", "satelite")) {
    descriptivo <- rbind(descriptivo,
                         data.frame(anio = anio, fuente = fuente,
                                    curvas = length(unique(largo$curva[filas])),
                                    estaciones = length(unique(largo$STATION_ID[filas])),
                                    describir(largo[[fuente]][filas])))
  }
}
descriptivo$pct_dias_faltantes <- ifelse(descriptivo$fuente == "estacion",
                                         faltantes_dias$pct_dias_faltantes[match(descriptivo$anio, faltantes_dias$YEAR)],
                                         NA)
descriptivo[descriptivo$anio == "Todos", ]
guardar_tabla(descriptivo, "03_descriptivo_estadisticas.csv")

# Coeficiente de variacion por año
cv_anual <- descriptivo[descriptivo$anio != "Todos", ]
plot(as.integer(cv_anual$anio[cv_anual$fuente == "estacion"]), cv_anual$CV_pct[cv_anual$fuente == "estacion"],
     type = "b", pch = 16, col = "darkgreen", ylim = range(cv_anual$CV_pct),
     xlab = "Año", ylab = "Coeficiente de variación (%)", main = "Coeficiente de variación por año")
lines(as.integer(cv_anual$anio[cv_anual$fuente == "satelite"]), cv_anual$CV_pct[cv_anual$fuente == "satelite"],
      type = "b", pch = 17, col = "steelblue")
legend("topright", c("Estación", "Satélite"), col = c("darkgreen", "steelblue"), pch = c(16, 17), lty = 1, bty = "n")
guardar_figura("03_descriptivo_CV_por_anio.png", 8, 5)

# Correlacion estacion-satelite por pentada y de la lluvia media de cada curva con la altitud
medias_curva <- largo %>%
  group_by(curva, YEAR) %>%
  summarize(estacion = mean(estacion, na.rm = TRUE), satelite = mean(satelite, na.rm = TRUE),
            altitud = first(altitud), .groups = "drop")
correlaciones <- function(p, m) {
  data.frame(estacion_satelite = cor(p$estacion, p$satelite, use = "complete.obs"),
             estacion_altitud = cor(m$estacion, m$altitud, use = "complete.obs"),
             satelite_altitud = cor(m$satelite, m$altitud, use = "complete.obs"))
}
correlacion <- data.frame(anio = "Todos", correlaciones(largo, medias_curva))
for (anio in sort(unique(largo$YEAR))) {
  correlacion <- rbind(correlacion,
                       data.frame(anio = as.character(anio),
                                  correlaciones(largo[largo$YEAR == anio, ], medias_curva[medias_curva$YEAR == anio, ])))
}
correlacion
guardar_tabla(correlacion, "03_descriptivo_correlacion.csv")

# Dispersion estacion frente a satelite por año
etiquetas_r <- data.frame(YEAR = as.integer(correlacion$anio[-1]),
                          texto = paste("r =", round(correlacion$estacion_satelite[-1], 2)))
grafico_dispersion <- ggplot(largo[!is.na(largo$estacion), ], aes(satelite, estacion)) +
  geom_point(size = 0.3, alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linetype = 2) +
  geom_text(data = etiquetas_r, aes(x = -Inf, y = Inf, label = texto), hjust = -0.1, vjust = 1.3, size = 3) +
  facet_wrap(~ YEAR, ncol = 6) +
  labs(x = "Satélite (mm por péntada)", y = "Estación (mm por péntada)",
       title = "Lluvia por péntada: estación frente a satélite por año") +
  theme_bw()
png("Figuras/03_descriptivo_dispersion_por_anio.png", width = 12, height = 9, units = "in", res = 150)
print(grafico_dispersion)
dev.off()

# Lluvia media de cada curva frente a la altitud
par(mfrow = c(1, 2))
plot(medias_curva$altitud, medias_curva$estacion, pch = 16, cex = 0.5, col = "darkgreen",
     xlab = "Altitud (m)", ylab = "Lluvia media (mm por péntada)", main = "Estación")
plot(medias_curva$altitud, medias_curva$satelite, pch = 16, cex = 0.5, col = "steelblue",
     xlab = "Altitud (m)", ylab = "Lluvia media (mm por péntada)", main = "Satélite")
guardar_figura("03_descriptivo_altitud.png", 10, 5)
par(mfrow = c(1, 1))


#################### Construccion de datos funcionales ####################
# Base B-spline y rango de la matriz
K_evaluar <- c(6, 8, 10, 12, 15, 20, 25, 30)

# Base B-spline cubica evaluada en las pentadas
construir_bspline <- function(t, K) {
  BasePrec <- create.bspline.basis(rangeval = c(0, 72), nbasis = K, norder = 4)
  eval.basis(t, BasePrec)
}

# Particion de la unidad
Phi_ejemplo <- construir_bspline(tiempos, 12)
dim(Phi_ejemplo)
range(rowSums(Phi_ejemplo))

# Rango de la base en las pentadas observadas
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


# Construccion por minimos cuadrados ordinarios (MCO)
# Cada curva usa solo sus pentadas observadas
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

# Seleccion de la base, el numero de bases y el metodo
criterios <- function(datos, Phi, pesos = rep(1, 72), R = 0, lambda = 0) {
  ECM <- c()
  VC <- c()
  VCG <- c()
  for (i in 1:nrow(datos)) {
    obs <- which(!is.na(datos[i, ]))
    y <- datos[i, obs]
    P <- Phi[obs, ]
    PW <- t(P * pesos[obs])
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

# B-spline: K de 4 a 15 (con mas K algunas curvas interpolan)
# Fourier: K impar. Monomial: K de 1 a 4
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

# Penalizacion por rugosidad
base_bspline <- create.bspline.basis(c(0, 72), nbasis = 72, norder = 4)
base_fourier <- create.fourier.basis(c(0, 72), nbasis = 71)
penalizadas <- list("B-spline" = list(Phi = eval.basis(tiempos, base_bspline),
                                      R = eval.penalty(base_bspline, 2)),
                    "Fourier" = list(Phi = eval.basis(tiempos, base_fourier),
                                     R = eval.penalty(base_fourier, vec2Lfd(c(0, (2 * pi / 72)^2, 0), c(0, 72)))))

for (nombre in names(penalizadas)) {
  Phi <- penalizadas[[nombre]]$Phi
  R <- penalizadas[[nombre]]$R
  for (log_lambda in seq(-1, 5, 0.25)) {
    tabla <- rbind(tabla,
                   data.frame(fuente = "estaciones", base = nombre, metodo = "Penalizacion", parametro = log_lambda,
                              t(criterios(Y_train, Phi, R = R, lambda = 10^log_lambda))),
                   data.frame(fuente = "satelite", base = nombre, metodo = "Penalizacion", parametro = log_lambda,
                              t(criterios(X_train, Phi, R = R, lambda = 10^log_lambda))))
  }
}
tabla
guardar_tabla(tabla, "04_construccion_criterios_todas_las_bases.csv")

# Mejor parametro de cada fuente, base y metodo
ablacion <- tabla[order(tabla$fuente, tabla$VCG), ]
ablacion <- ablacion[!duplicated(ablacion[, c("fuente", "base", "metodo")]), ]
ablacion
guardar_tabla(ablacion, "04_construccion_ablacion.csv")

# VCG de MCO y MCP
par(mfrow = c(2, 3))
for (f in c("estaciones", "satelite")) {
  for (nombre in names(valores_K)) {
    mco <- tabla[tabla$fuente == f & tabla$base == nombre & tabla$metodo == "MCO", ]
    mcp <- tabla[tabla$fuente == f & tabla$base == nombre & tabla$metodo == "MCP", ]
    plot(mco$parametro, mco$VCG, type = "b", pch = 19, xlab = "K", ylab = "VCG",
         main = paste(nombre, "-", f), ylim = range(mco$VCG, mcp$VCG))
    lines(mcp$parametro, mcp$VCG, type = "b", pch = 17, col = "red")
    legend("topright", c("MCO", "MCP"), col = c("black", "red"), pch = c(19, 17), bty = "n")
  }
}
guardar_figura("04_construccion_VCG_MCO_MCP.png", 12, 7)
par(mfrow = c(1, 1))

# VCG de la penalizacion
par(mfrow = c(1, 2))
for (f in c("estaciones", "satelite")) {
  b <- tabla[tabla$fuente == f & tabla$base == "B-spline" & tabla$metodo == "Penalizacion", ]
  fo <- tabla[tabla$fuente == f & tabla$base == "Fourier" & tabla$metodo == "Penalizacion", ]
  plot(b$parametro, b$VCG, type = "b", pch = 19, xlab = expression(log[10](lambda)), ylab = "VCG",
       main = paste("Penalizacion -", f), ylim = range(b$VCG, fo$VCG))
  lines(fo$parametro, fo$VCG, type = "b", pch = 17, col = "red")
  legend("topright", c("B-spline", "Fourier"), col = c("black", "red"), pch = c(19, 17), bty = "n")
}
guardar_figura("04_construccion_VCG_penalizacion.png", 10, 5)
par(mfrow = c(1, 1))


# Construccion elegida: B-spline cubica con MCP y K de menor VCG
estaciones_MCP <- tabla[tabla$fuente == "estaciones" & tabla$base == "B-spline" & tabla$metodo == "MCP", ]
satelite_MCP <- tabla[tabla$fuente == "satelite" & tabla$base == "B-spline" & tabla$metodo == "MCP", ]
K_estaciones <- estaciones_MCP$parametro[which.min(estaciones_MCP$VCG)]
K_satelite <- satelite_MCP$parametro[which.min(satelite_MCP$VCG)]
K_estaciones
K_satelite

# ECM, VC y VCG frente a K
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
guardar_figura("04_construccion_criterios_K.png", 12, 7)
par(mfrow = c(1, 1))

# Base de cada fuente y pesos de MCP (entrenamiento)
base_estaciones <- crear_base("B-spline", K_estaciones)
base_satelite <- crear_base("B-spline", K_satelite)
Phi_estaciones <- eval.basis(tiempos, base_estaciones)
Phi_satelite <- eval.basis(tiempos, base_satelite)
pesos_estaciones <- pesos_MCP(Y_train, Phi_estaciones)
pesos_satelite <- pesos_MCP(X_train, Phi_satelite)
guardar_tabla(data.frame(pentada = tiempos, estaciones = pesos_estaciones, satelite = pesos_satelite),
              "04_construccion_pesos_MCP.csv")

# Ajuste MCP de cada curva
ajuste_MCP <- function(datos, Phi, pesos) {
  coeficientes <- matrix(NA, nrow(datos), ncol(Phi), dimnames = list(rownames(datos), NULL))
  ajustes <- matrix(NA, nrow(datos), 72, dimnames = dimnames(datos))
  for (i in 1:nrow(datos)) {
    obs <- which(!is.na(datos[i, ]))
    P <- Phi[obs, ]
    PW <- t(P * pesos[obs])
    coeficientes[i, ] <- solve(PW %*% P) %*% PW %*% datos[i, obs]
    ajustes[i, ] <- Phi %*% coeficientes[i, ]
  }
  list(coeficientes = coeficientes, ajustes = ajustes, residuos = datos - ajustes)
}

construccion_estaciones <- ajuste_MCP(Y_train, Phi_estaciones, pesos_estaciones)
construccion_satelite <- ajuste_MCP(X_train, Phi_satelite, pesos_satelite)
ajustes_estaciones <- construccion_estaciones$ajustes
ajustes_satelite <- construccion_satelite$ajustes
residuos_estaciones <- construccion_estaciones$residuos
residuos_satelite <- construccion_satelite$residuos

# Curvas 1981-2003 con la base y los pesos elegidos
construccion_completa <- list(
  terrestre = ajuste_MCP(Y_completa, Phi_estaciones, pesos_estaciones),
  satelital = ajuste_MCP(X_completa, Phi_satelite, pesos_satelite))

# Parametros elegidos
parametros_construccion <- data.frame(
  fuente = c("estaciones", "satelite"), base = "B-spline", metodo = "MCP", orden = 4,
  K = c(K_estaciones, K_satelite),
  VCG = c(min(estaciones_MCP$VCG), min(satelite_MCP$VCG)))
parametros_construccion
guardar_tabla(parametros_construccion, "04_construccion_K_elegido.csv")

# Verificacion de la matriz suavizadora con una curva completa
i <- which(rowSums(is.na(Y_train)) == 0)[1]
verificacion_S <- data.frame()
for (f in c("estaciones", "satelite")) {
  if (f == "estaciones") {
    y <- Y_train[i, ]
    Phi <- Phi_estaciones
    base <- base_estaciones
    pesos <- pesos_estaciones
  } else {
    y <- X_train[i, ]
    Phi <- Phi_satelite
    base <- base_satelite
    pesos <- pesos_satelite
  }
  PW <- t(Phi * pesos)
  y2cMap <- solve(PW %*% Phi) %*% PW
  S <- Phi %*% y2cMap

  # VC por el atajo y dejando una pentada por fuera
  VC_atajo <- mean(((y - S %*% y) / (1 - diag(S)))^2)
  prediccion <- c()
  for (j in 1:72) {
    PW <- t(Phi[-j, ] * pesos[-j])
    coeficientes <- solve(PW %*% Phi[-j, ]) %*% PW %*% y[-j]
    prediccion[j] <- Phi[j, ] %*% coeficientes
  }
  VC_dejar_uno <- mean((y - prediccion)^2)

  # Comparacion con smooth.basis
  ajuste_fda <- smooth.basis(tiempos, y, base, wtvec = pesos)

  verificacion_S <- rbind(verificacion_S,
                          data.frame(fuente = f, curva = rownames(Y_train)[i],
                                     idempotencia = max(abs(S %*% S - S)),
                                     traza = sum(diag(S)), K = ncol(Phi),
                                     diferencia_coeficientes_fda = max(abs(y2cMap %*% y - ajuste_fda$fd$coefs)),
                                     VC_atajo = VC_atajo, VC_dejar_uno = VC_dejar_uno))
}
verificacion_S
guardar_tabla(verificacion_S, "04_construccion_verificacion_S.csv")
# Propiedades de la matriz suavizadora


# Residuos de la construccion elegida
# Autocorrelacion de rezago 1 en las curvas completas
completas_estaciones <- rowSums(is.na(Y_train)) == 0
completas_satelite <- rowSums(is.na(X_train)) == 0
r1_estaciones <- apply(residuos_estaciones[completas_estaciones, ], 1,
                       function(e) acf(e, lag.max = 1, plot = FALSE)$acf[2])
r1_satelite <- apply(residuos_satelite[completas_satelite, ], 1,
                     function(e) acf(e, lag.max = 1, plot = FALSE)$acf[2])

limite <- 1.96 / sqrt(72)
autocorrelacion <- data.frame(fuente = c("estaciones", "satelite"),
                              curvas = c(length(r1_estaciones), length(r1_satelite)),
                              mediana_r1 = c(median(r1_estaciones), median(r1_satelite)),
                              pct_fuera_de_banda = c(mean(abs(r1_estaciones) > limite) * 100,
                                                     mean(abs(r1_satelite) > limite) * 100))
autocorrelacion
guardar_tabla(autocorrelacion, "04_construccion_autocorrelacion_residuos.csv")

par(mfrow = c(1, 2))
hist(r1_estaciones, breaks = 30, main = "Estaciones", xlab = "Autocorrelacion de rezago 1")
abline(v = c(-limite, limite), lty = 2, col = "red")
hist(r1_satelite, breaks = 30, main = "Satelite", xlab = "Autocorrelacion de rezago 1")
abline(v = c(-limite, limite), lty = 2, col = "red")
guardar_figura("04_construccion_autocorrelacion_residuos.png", 10, 5)
par(mfrow = c(1, 1))

# Heterocedasticidad
varianza_estaciones <- apply(residuos_estaciones, 2, var, na.rm = TRUE)
varianza_satelite <- apply(residuos_satelite, 2, var, na.rm = TRUE)

par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
matplot(tiempos, t(residuos_estaciones), type = "l", lty = 1, col = "grey",
        xlab = "Pentada", ylab = "Residuo", main = "Residuos - estaciones")
matplot(tiempos, t(residuos_satelite), type = "l", lty = 1, col = "grey",
        xlab = "Pentada", ylab = "Residuo", main = "Residuos - satélite")
plot(tiempos, varianza_estaciones, type = "l",
     xlab = "Pentada", ylab = "Varianza", main = "Varianza residual - estaciones")
plot(tiempos, varianza_satelite, type = "l",
     xlab = "Pentada", ylab = "Varianza", main = "Varianza residual - satélite")
guardar_figura("04_construccion_residuos_por_pentada.png", 10, 8)
par(mfrow = c(1, 1))

max(varianza_estaciones) / min(varianza_estaciones)
max(varianza_satelite) / min(varianza_satelite)
# La varianza de los residuos cambia en el año, sobre todo en las temporadas de lluvia


# Valores negativos en las curvas ajustadas
valores_negativos <- data.frame(fuente = c("estaciones", "satelite"),
                                pct_valores = c(mean(ajustes_estaciones < 0), mean(ajustes_satelite < 0)) * 100,
                                curvas_afectadas = c(sum(rowSums(ajustes_estaciones < 0) > 0),
                                                     sum(rowSums(ajustes_satelite < 0) > 0)),
                                minimo_mm = c(min(ajustes_estaciones), min(ajustes_satelite)))
valores_negativos
guardar_tabla(valores_negativos, "04_construccion_valores_negativos.csv")


# Varianza residual y bandas
# Varianza residual de cada curva
varianza_residual <- data.frame(
  curva = rownames(Y_train),
  estaciones = rowSums(residuos_estaciones^2, na.rm = TRUE) /
    (rowSums(!is.na(Y_train)) - K_estaciones),
  satelite = rowSums(residuos_satelite^2, na.rm = TRUE) /
    (rowSums(!is.na(X_train)) - K_satelite))
summary(varianza_residual[, c("estaciones", "satelite")])

boxplot(varianza_residual[, c("estaciones", "satelite")], names = c("Estaciones", "Satélite"), log = "y",
        col = c("lightblue", "moccasin"),
        main = "Varianza residual por fuente", ylab = "Varianza (mm^2, escala log)")
guardar_figura("04_construccion_varianza_residual_fuentes.png", 6, 5)

# Comparacion de las varianzas para la misma curva estacion-año
limites_varianza <- range(varianza_residual$estaciones, varianza_residual$satelite)
plot(varianza_residual$estaciones, varianza_residual$satelite,
     log = "xy", xlim = limites_varianza, ylim = limites_varianza, pch = 16, cex = 0.7, col = "steelblue",
     main = "Varianza residual por par estación-satélite",
     xlab = "Estaciones (mm^2, escala log)", ylab = "Satélite (mm^2, escala log)")
abline(0, 1, lty = 2, col = "darkorange")
guardar_figura("04_construccion_varianza_residual_pares.png", 6, 6)

# Varianza residual por pentada
sigma2_estaciones <- colMeans(residuos_estaciones^2, na.rm = TRUE)
sigma2_satelite <- colMeans(residuos_satelite^2, na.rm = TRUE)

# Covarianza residual entre pentadas: promedio de e e^T sobre las curvas (solo pares de pentadas observadas)
covarianza_residual <- function(residuos) {
  observado <- !is.na(residuos)
  residuos[!observado] <- 0
  crossprod(residuos) / crossprod(observado)
}
Sigma_estaciones <- covarianza_residual(residuos_estaciones)
Sigma_satelite <- covarianza_residual(residuos_satelite)

par(mfrow = c(1, 2))
for (f in c("estaciones", "satelite")) {
  Sigma <- if (f == "estaciones") Sigma_estaciones else Sigma_satelite
  limite_z <- max(abs(Sigma))
  image(tiempos, tiempos, Sigma, zlim = c(-limite_z, limite_z), col = hcl.colors(41, "Blue-Red 3"),
        xlab = "Pentada", ylab = "Pentada", main = paste("Covarianza residual -", f))
}
guardar_figura("04_construccion_covarianza_residual.png", 11, 5)
par(mfrow = c(1, 1))

# Correlacion residual media segun la distancia entre pentadas
resumen_covarianza <- data.frame()
for (f in c("estaciones", "satelite")) {
  Sigma <- if (f == "estaciones") Sigma_estaciones else Sigma_satelite
  Rho <- Sigma / sqrt(diag(Sigma) %o% diag(Sigma))
  distancia <- abs(outer(tiempos, tiempos, "-"))
  resumen_covarianza <- rbind(resumen_covarianza, data.frame(
    fuente = f, varianza_media = mean(diag(Sigma)),
    correlacion_rezago_1 = mean(Rho[distancia == 1]), correlacion_rezago_2 = mean(Rho[distancia == 2]),
    correlacion_rezago_3_o_mas = mean(Rho[distancia >= 3])))
}
resumen_covarianza
guardar_tabla(resumen_covarianza, "04_construccion_covarianza_residual.csv")

# Varianza de la curva ajustada
varianza_ajuste <- function(y, Phi, pesos, sigma2) {
  obs <- which(!is.na(y))
  P <- Phi[obs, ]
  PW <- t(P * pesos[obs])
  PhiC <- Phi %*% solve(PW %*% P) %*% PW
  as.vector(PhiC^2 %*% sigma2[obs])
}

# Q-Q de los residuos
par(mfrow = c(2, 3))
for (j in c(2, 40, 60)) {
  qqnorm(residuos_estaciones[, j], pch = 16, cex = 0.5, main = paste("Estación - Péntada", j),
         xlab = "Cuantiles normales", ylab = "Residuos (mm)")
  qqline(residuos_estaciones[, j], col = "steelblue", lwd = 2)
}
for (j in c(2, 40, 60)) {
  qqnorm(residuos_satelite[, j], pch = 16, cex = 0.5, main = paste("Satélite - Péntada", j),
         xlab = "Cuantiles normales", ylab = "Residuos (mm)")
  qqline(residuos_satelite[, j], col = "darkorange", lwd = 2)
}
guardar_figura("04_construccion_QQ_residuos.png", 12, 7)
par(mfrow = c(1, 1))

# Pares estacion-año para elegir las curvas de las bandas
tabla_pares_bandas <- data.frame(
  curva = rownames(Y_train),
  acumulado_estacion = rowSums(ajustes_estaciones),
  acumulado_satelite = rowSums(ajustes_satelite),
  RMSE_estacion = sqrt(rowMeans(residuos_estaciones^2, na.rm = TRUE)))

# Distancia a la mediana del acumulado
percentil_estacion <- 100 * (rank(tabla_pares_bandas$acumulado_estacion) - 0.5) / nrow(tabla_pares_bandas)
percentil_satelite <- 100 * (rank(tabla_pares_bandas$acumulado_satelite) - 0.5) / nrow(tabla_pares_bandas)
tabla_pares_bandas$distancia_mediana <- abs(percentil_estacion - 50) + abs(percentil_satelite - 50)

# Pares cercanos a la mediana y con mayor RMSE
head(tabla_pares_bandas[order(tabla_pares_bandas$distancia_mediana), ], 5)
head(tabla_pares_bandas[order(-tabla_pares_bandas$RMSE_estacion), ], 5)

# Bandas punto a punto
curvas_bandas <- c("2610000111_1989", "5421200110_1982")
for (curva in curvas_bandas) {
  par(mfrow = c(1, 2))
  y <- Y_train[curva, ]
  x <- X_train[curva, ]
  EE_estacion <- sqrt(varianza_ajuste(y, Phi_estaciones, pesos_estaciones, sigma2_estaciones))
  EE_satelite <- sqrt(varianza_ajuste(x, Phi_satelite, pesos_satelite, sigma2_satelite))
  limites <- range(y, x, ajustes_estaciones[curva, ] + 2 * EE_estacion, ajustes_estaciones[curva, ] - 2 * EE_estacion,
                   ajustes_satelite[curva, ] + 2 * EE_satelite, ajustes_satelite[curva, ] - 2 * EE_satelite, na.rm = TRUE)
  
  plot(tiempos, y, pch = 16, cex = 0.65, ylim = limites, xlab = "Péntada", ylab = "Precipitación (mm)",
       main = paste("Estación", curva))
  lines(tiempos, ajustes_estaciones[curva, ], col = "steelblue", lwd = 2)
  lines(tiempos, ajustes_estaciones[curva, ] - 2 * EE_estacion, col = "steelblue", lty = 2)
  lines(tiempos, ajustes_estaciones[curva, ] + 2 * EE_estacion, col = "steelblue", lty = 2)
  legend("topright", c("Observado", "Ajuste MCP", "+- 2 EE"), col = c("black", "steelblue", "steelblue"),
         pch = c(16, NA, NA), lty = c(NA, 1, 2), bty = "n", cex = 0.8)
  
  plot(tiempos, x, pch = 16, cex = 0.65, ylim = limites, xlab = "Péntada", ylab = "Precipitación (mm)",
       main = paste("Satélite", curva))
  lines(tiempos, ajustes_satelite[curva, ], col = "darkorange", lwd = 2)
  lines(tiempos, ajustes_satelite[curva, ] - 2 * EE_satelite, col = "darkorange", lty = 2)
  lines(tiempos, ajustes_satelite[curva, ] + 2 * EE_satelite, col = "darkorange", lty = 2)
  legend("topright", c("Observado", "Ajuste MCP", "+- 2 EE"), col = c("black", "darkorange", "darkorange"),
         pch = c(16, NA, NA), lty = c(NA, 1, 2), bty = "n", cex = 0.8)
  guardar_figura(paste0("04_construccion_bandas_", curva, ".png"), 11, 5)
}
par(mfrow = c(1, 1))


# Sensibilidad del suavizado en la diferencia estacion - satelite
# Diferencias originales y funcionales (entrenamiento)
D_original <- Y_train - X_train
D_funcional <- ajustes_estaciones - ajustes_satelite

# Diferencia media por pentada
Media_original <- colMeans(D_original, na.rm = TRUE)
Media_funcional <- colMeans(D_funcional)

par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
plot(1:72, Media_original, type = "l", lwd = 2, col = "black",
     ylim = range(Media_original, Media_funcional),
     xlab = "Pentada", ylab = "Estacion - satelite (mm)",
     main = "Sensibilidad de la diferencia funcional")
lines(1:72, Media_funcional, col = "red", lwd = 2)
abline(h = 0, lty = 3)
legend("topright", c("Datos originales", "B-spline MCP"), col = c("black", "red"), lwd = 2, bty = "n")
guardar_figura("04_construccion_sensibilidad_diferencia.png", 8, 5)

# Efecto del suavizado y concordancia de signo
Cambio_suavizado <- Media_funcional - Media_original
mismo_signo <- sign(Media_original) == sign(Media_funcional)

resumen_sensibilidad <- data.frame(
  Periodo = c("Inicio (1-6)", "Interior (7-66)", "Final (67-72)", "Anio completo"),
  Diferencia_absoluta_media = c(mean(abs(Cambio_suavizado[1:6])), mean(abs(Cambio_suavizado[7:66])),
                                mean(abs(Cambio_suavizado[67:72])), mean(abs(Cambio_suavizado))),
  Concordancia_signo_pct = c(mean(mismo_signo[1:6]), mean(mismo_signo[7:66]),
                             mean(mismo_signo[67:72]), mean(mismo_signo)) * 100)
resumen_sensibilidad
guardar_tabla(resumen_sensibilidad, "04_construccion_sensibilidad_diferencia.csv")

# Pentadas con mayores diferencias por el suavizado
tabla_cambios <- data.frame(pentada = 1:72, original = Media_original, funcional = Media_funcional,
                            cambio = Cambio_suavizado, mismo_signo = mismo_signo)
head(tabla_cambios[order(-abs(tabla_cambios$cambio)), ], 10)


#################### Análisis exploratorio funcional ####################
# Curvas de 1981 a 2003
metadatos_eda <- metadatos_completos
fuentes_eda <- list(terrestre = Y_completa, satelital = X_completa)
bases_eda <- list(terrestre = base_estaciones, satelital = base_satelite)

# Coeficientes por columnas (una curva por columna)
coeficientes_eda <- list(terrestre = t(construccion_completa$terrestre$coeficientes),
                         satelital = t(construccion_completa$satelital$coeficientes))

# Cobertura anual de los pares con ajuste en ambas fuentes
pares_eda <- colnames(coeficientes_eda$terrestre)[colSums(!is.finite(coeficientes_eda$terrestre)) == 0 &
                                                    colSums(!is.finite(coeficientes_eda$satelital)) == 0]
metadatos_eda <- metadatos_eda[match(pares_eda, metadatos_eda$curva), ]
cobertura_eda <- cobertura[, c("YEAR", "apto")] %>%
  group_by(YEAR) %>%
  summarize(pares_aptos = sum(apto), .groups = "drop")
cobertura_eda$pares_eda <- as.integer(table(factor(metadatos_eda$YEAR, levels = cobertura_eda$YEAR)))
as.data.frame(cobertura_eda)
guardar_tabla(as.data.frame(cobertura_eda), "05_exploratorio_cobertura_por_anio.csv")

# Pares estación-año, estaciones y píxeles incluidos en el EDA
resumen_cobertura_eda <- data.frame(
  pares_estacion_anio = nrow(metadatos_eda),
  estaciones = n_distinct(metadatos_eda$STATION_ID),
  pixeles = n_distinct(metadatos_eda$pixel)
)
resumen_cobertura_eda
guardar_tabla(resumen_cobertura_eda, "05_exploratorio_resumen_cobertura.csv")

# Curvas emparejadas en la malla del EDA
pentadas_eda <- seq(min(tiempos), max(tiempos), by = 0.5)
curvas_eda_fd <- list()
valores_eda <- list()
for (fuente in names(fuentes_eda)) {
  fuentes_eda[[fuente]] <- fuentes_eda[[fuente]][pares_eda, , drop = FALSE]
  coeficientes_eda[[fuente]] <- coeficientes_eda[[fuente]][, pares_eda, drop = FALSE]
  curvas_eda_fd[[fuente]] <- fd(coeficientes_eda[[fuente]], bases_eda[[fuente]],
                                list("Péntada", pares_eda, "Precipitación (mm)"))
  valores_eda[[fuente]] <- eval.fd(pentadas_eda, curvas_eda_fd[[fuente]])
}

# Precipitación observada por péntada
limites_observados_eda <- range(unlist(fuentes_eda), finite = TRUE)
colores_curvas_eda <- c(terrestre = "steelblue", satelital = "darkorange")
par_curvas_eda <- par(mfrow = c(1, 2), bg = "white")
dev.control(displaylist = "enable")
matplot(tiempos, t(fuentes_eda$satelital), type = "o", pch = 1, cex = 0.35, lty = 1,
        col = adjustcolor(colores_curvas_eda[["satelital"]], alpha.f = 0.12),
        ylim = limites_observados_eda, xlab = "Péntada", ylab = "Precipitación observada (mm)",
        main = "Satélite")
matplot(tiempos, t(fuentes_eda$terrestre), type = "o", pch = 1, cex = 0.35, lty = 1,
        col = adjustcolor(colores_curvas_eda[["terrestre"]], alpha.f = 0.12),
        ylim = limites_observados_eda, xlab = "Péntada", ylab = "Precipitación observada (mm)",
        main = "Estaciones")
guardar_figura("05_exploratorio_observados.png", 12, 5)
# Los valores observados de estaciones alcanzan picos más altos que los del satélite

# Curvas ajustadas de precipitación por péntada
limites_curvas_eda <- range(unlist(valores_eda))
matplot(pentadas_eda, valores_eda$satelital, type = "l", lty = 1,
        col = adjustcolor(colores_curvas_eda[["satelital"]], alpha.f = 0.12),
        ylim = limites_curvas_eda, xlab = "Péntada", ylab = "Precipitación ajustada (mm)",
        main = "Satélite")
matplot(pentadas_eda, valores_eda$terrestre, type = "l", lty = 1,
        col = adjustcolor(colores_curvas_eda[["terrestre"]], alpha.f = 0.12),
        ylim = limites_curvas_eda, xlab = "Péntada", ylab = "Precipitación ajustada (mm)",
        main = "Estaciones")
guardar_figura("05_exploratorio_curvas_ajustadas.png", 12, 5)
par(par_curvas_eda)
# Las curvas de estaciones abarcan un rango más amplio, con algunos ajustes por debajo de 0 mm

# Media funcional y comprobación
media_eda <- data.frame(pentada = pentadas_eda)
coeficientes_medios_eda <- list()
media_eda_fd <- list()
error_media_eda <- numeric()
for (fuente in names(curvas_eda_fd)) {
  coefs <- coeficientes_eda[[fuente]]
  base_eda <- bases_eda[[fuente]]
  
  coeficientes_medios_eda[[fuente]] <- rowMeans(coefs)
  Phi_media_eda <- eval.basis(media_eda$pentada, base_eda)
  media_eda[[fuente]] <- as.vector(Phi_media_eda %*% coeficientes_medios_eda[[fuente]])
  
  media_eda_fd[[fuente]] <- mean.fd(curvas_eda_fd[[fuente]])
  media_libreria <- eval.fd(media_eda$pentada, media_eda_fd[[fuente]])
  error_media_eda[fuente] <- max(abs(media_eda[[fuente]] - media_libreria))
}
# Diferencia absoluta máxima entre la media manual y mean.fd() (mm)
error_media_eda
guardar_tabla(data.frame(fuente = names(error_media_eda), diferencia_mm = error_media_eda),
              "05_exploratorio_verificacion_media.csv")

# Media funcional de la precipitación por fuente
plot(media_eda$pentada, media_eda$terrestre, type = "l", lwd = 2, col = "steelblue",
     ylim = range(media_eda$terrestre, media_eda$satelital),
     xlab = "Péntada", ylab = "Precipitación media (mm)", main = "Media funcional de la precipitación")
lines(media_eda$pentada, media_eda$satelital, lwd = 2, lty = 2, col = "darkorange")
legend("topleft", c("Estaciones", "Satélite"), col = c("steelblue", "darkorange"),
       lty = c(1, 2), lwd = 2, bty = "n")
guardar_figura("05_exploratorio_media.png", 8, 5)
# Las dos fuentes muestran picos de lluvia en abril-mayo y en noviembre, con un valle en julio-agosto
# El pico de noviembre es el más alto del año en ambas fuentes
# Los dos máximos separados por un periodo menos lluvioso muestran un patrón bimodal en la media
# La media resume los pares estación-año disponibles; las estaciones con más años aportan más curvas


# Varianza y desviación estándar entre curvas estación-año
varianza_eda <- data.frame(pentada = media_eda$pentada)
desviacion_eda <- data.frame(pentada = media_eda$pentada)
error_varianza_eda <- numeric()
error_desviacion_eda <- numeric()
resumen_dispersion_eda <- data.frame()
for (fuente in names(curvas_eda_fd)) {
  n_curvas_eda <- ncol(valores_eda[[fuente]])
  diferencias_eda <- valores_eda[[fuente]] - media_eda[[fuente]]
  varianza_eda[[fuente]] <- apply(diferencias_eda^2, 1, sum) / (n_curvas_eda - 1)
  desviacion_eda[[fuente]] <- sqrt(varianza_eda[[fuente]])
  
  # Comprobación con la varianza y la desviación estándar de R en cada punto
  varianza_R <- apply(valores_eda[[fuente]], 1, var)
  desviacion_R <- apply(valores_eda[[fuente]], 1, sd)
  error_varianza_eda[fuente] <- max(abs(varianza_eda[[fuente]] - varianza_R))
  error_desviacion_eda[fuente] <- max(abs(desviacion_eda[[fuente]] - desviacion_R))
  
  minimo <- which.min(desviacion_eda[[fuente]])
  maximo <- which.max(desviacion_eda[[fuente]])
  resumen_dispersion_eda <- rbind(resumen_dispersion_eda,
                                  data.frame(fuente = fuente, curvas = n_curvas_eda,
                                             pentada_min = media_eda$pentada[minimo], DE_min_mm = desviacion_eda[[fuente]][minimo],
                                             pentada_max = media_eda$pentada[maximo], DE_max_mm = desviacion_eda[[fuente]][maximo]))
}
# Diferencia máxima entre la varianza manual y var() de R (mm^2)
error_varianza_eda
# Diferencia máxima entre la desviación estándar manual y sd() de R (mm)
error_desviacion_eda
# Resumen de la dispersión por fuente
resumen_dispersion_eda
guardar_tabla(resumen_dispersion_eda, "05_exploratorio_resumen_dispersion.csv")

# Dispersión de la precipitación por fuente
par_dispersion_eda <- par(mfrow = c(1, 2))
plot(varianza_eda$pentada, varianza_eda$terrestre, type = "l", lwd = 2, col = "steelblue",
     ylim = c(0, max(varianza_eda$terrestre, varianza_eda$satelital)),
     xlab = "Péntada", ylab = expression("Varianza (" * mm^2 * ")"), main = "Varianza funcional")
lines(varianza_eda$pentada, varianza_eda$satelital, lwd = 2, lty = 2, col = "darkorange")
legend("topleft", c("Estaciones", "Satélite"), col = c("steelblue", "darkorange"),
       lty = c(1, 2), lwd = 2, bty = "n")
# En toda la gráfica, las curvas de estaciones se dispersan más alrededor de su media que las del satélite
# La dispersión reúne diferencias entre lugares y entre años en una misma péntada
# Este cálculo no separa las curvas por altitud ni por fase ENSO, por lo que no identifica
# qué parte de la dispersión corresponde a diferencias entre estaciones o entre fases
plot(desviacion_eda$pentada, desviacion_eda$terrestre, type = "l", lwd = 2, col = "steelblue",
     ylim = c(0, max(desviacion_eda$terrestre, desviacion_eda$satelital)),
     xlab = "Péntada", ylab = "Desviación estándar (mm)", main = "Desviación estándar funcional")
lines(desviacion_eda$pentada, desviacion_eda$satelital, lwd = 2, lty = 2, col = "darkorange")
legend("topleft", c("Estaciones", "Satélite"), col = c("steelblue", "darkorange"),
       lty = c(1, 2), lwd = 2, bty = "n")
# La desviación estándar muestra la misma dispersión en mm y aumenta hacia noviembre (en ambas fuentes)
# Ese aumento coincide con el segundo pico de la media y señala mayor diferencia entre curvas en esa temporada
# El máximo del satélite es 16,15 mm en la péntada 62; el de estaciones es 26,20 mm en la primera péntada
# Entonces, dos picos de lluvia media no implican que la mayor dispersión ocurra siempre en esos dos picos
guardar_figura("05_exploratorio_varianza_desviacion.png", 12, 5)
par(par_dispersion_eda)

# Curvas de cada año con su media y su desviación estándar, un panel por año (Ospina Galíndez, 2016)
# El título da la fase ENSO del año según el ONI promedio: El Niño desde 0,5, La Niña desde -0,5 y Neutral entre ambos
# Todos los paneles de una fuente tienen la misma escala, hasta su valor máximo y como mucho 350 mm;
# en estaciones solo sale del eje el inicio de una curva de 1982 (5421200110), con un borde de 528 mm
anios_eda <- sort(unique(metadatos_eda$YEAR))
oni_anio_eda <- oni_anual$oni_promedio[match(anios_eda, oni_anual$YEAR)]
fase_anio_eda <- ifelse(oni_anio_eda >= 0.5, "El Niño", ifelse(oni_anio_eda <= -0.5, "La Niña", "Neutral"))
nombres_fuente_eda <- c(terrestre = "Estaciones", satelital = "Satélite")
for (fuente in names(valores_eda)) {
  limite_y <- c(min(valores_eda[[fuente]]), min(max(valores_eda[[fuente]]), 350))
  png(paste0("Figuras/05_exploratorio_curvas_por_anio_", fuente, ".png"), width = 18, height = 11, units = "in", res = 150)
  par(mfrow = c(4, 6), mar = c(3.5, 3.5, 2.5, 0.5), mgp = c(2, 0.6, 0), oma = c(0, 0, 3, 0))
  for (k in seq_along(anios_eda)) {
    curvas_anio <- valores_eda[[fuente]][, metadatos_eda$YEAR == anios_eda[k], drop = FALSE]
    matplot(pentadas_eda, curvas_anio, type = "l", lty = 1, col = "grey85", ylim = limite_y,
            xlab = "Péntada", ylab = "Precipitación (mm)",
            main = paste0(anios_eda[k], ": ", fase_anio_eda[k], " (ONI ",
                          format(round(oni_anio_eda[k], 2), nsmall = 2, decimal.mark = ","), ")"))
    lines(pentadas_eda, rowMeans(curvas_anio), lwd = 2.5, col = "red")
    lines(pentadas_eda, apply(curvas_anio, 1, sd), lwd = 2.5, lty = 2, col = "darkblue")
  }
  plot.new()
  legend("center", c("Curvas del año", "Media", "Desviación estándar"), col = c("grey85", "red", "darkblue"),
         lty = c(1, 1, 2), lwd = c(1, 2.5, 2.5), bty = "n", cex = 1.5)
  mtext(paste0(nombres_fuente_eda[fuente], ": curvas por año, media y desviación estándar",
               if (max(valores_eda[[fuente]]) > 350) " (eje hasta 350 mm)" else ""), outer = TRUE, cex = 1.4)
  dev.off()
}
# En el satélite el régimen bimodal se repite todos los años, con picos cerca de las péntadas 25 y 60 y poca dispersión
# En estaciones la media también tiene dos temporadas, pero es más plana: la desviación estándar es parecida a la media
# porque cada año unas pocas estaciones llueven mucho más que el resto
# La fase ENSO no cambia la forma bimodal; las diferencias entre años están en la altura y la anchura de las temporadas,
# como en 1999 (La Niña), cuando la primera temporada es más ancha y la lluvia media de estaciones es alta desde enero


# Covarianza y correlación funcionales de cada fuente
covarianza_eda <- list()
correlacion_eda <- list()
comprobacion_cov_cor_eda <- data.frame()
for (fuente in names(curvas_eda_fd)) {
  coefs_centrados_eda <- coeficientes_eda[[fuente]] - coeficientes_medios_eda[[fuente]]
  covarianza_coefs_eda <- coefs_centrados_eda %*% t(coefs_centrados_eda) /
    (ncol(coefs_centrados_eda) - 1)
  Phi_eda <- eval.basis(pentadas_eda, bases_eda[[fuente]])
  covarianza_eda[[fuente]] <- Phi_eda %*% covarianza_coefs_eda %*% t(Phi_eda)
  correlacion_eda[[fuente]] <- covarianza_eda[[fuente]] /
    outer(desviacion_eda[[fuente]], desviacion_eda[[fuente]])
  
  covarianza_R <- eval.bifd(pentadas_eda, pentadas_eda, var.fd(curvas_eda_fd[[fuente]]))
  correlacion_R <- cor.fd(pentadas_eda, curvas_eda_fd[[fuente]])
  comprobacion_cov_cor_eda <- rbind(comprobacion_cov_cor_eda,
                                    data.frame(fuente = fuente,
                                               diferencia_covarianza_mm2 = max(abs(covarianza_eda[[fuente]] - covarianza_R)),
                                               diferencia_correlacion = max(abs(correlacion_eda[[fuente]] - correlacion_R)),
                                               diferencia_diagonal_mm2 = max(abs(diag(covarianza_eda[[fuente]]) - varianza_eda[[fuente]]))))
}
# Diferencias frente a las funciones de R y frente a la varianza en la diagonal
comprobacion_cov_cor_eda
guardar_tabla(comprobacion_cov_cor_eda, "05_exploratorio_verificacion_covarianza_correlacion.csv")

# Curvas de nivel de covarianza y correlación por fuente
niveles_cov_eda <- pretty(range(unlist(covarianza_eda)), n = 12)
niveles_cor_eda <- seq(-1, 1, by = 0.1)
par_cov_cor_eda <- par(mfrow = c(1, 2))
contour(pentadas_eda, pentadas_eda, covarianza_eda$terrestre, levels = niveles_cov_eda,
        xlab = "Péntada s", ylab = "Péntada t", main = "Covarianza de estaciones (mm^2)")
# Cada línea une cruces con igual covarianza; en la diagonal se lee la varianza de cada péntada
# Las estaciones muestran covarianzas altas en los cruces entre péntadas cercanas a noviembre
contour(pentadas_eda, pentadas_eda, covarianza_eda$satelital, levels = niveles_cov_eda,
        xlab = "Péntada s", ylab = "Péntada t", main = "Covarianza de satélite (mm^2)")
# El satélite también muestra covarianzas altas alrededor de noviembre, en una escala de dispersión menor
guardar_figura("05_exploratorio_covarianza.png", 12, 6)
contour(pentadas_eda, pentadas_eda, correlacion_eda$terrestre, levels = niveles_cor_eda,
        xlab = "Péntada s", ylab = "Péntada t", main = "Correlación de estaciones")
# Las correlaciones más altas se ven cerca de la diagonal: los pares más lluviosos en una péntada tienden
# a seguir entre los más lluviosos en las péntadas cercanas
contour(pentadas_eda, pentadas_eda, correlacion_eda$satelital, levels = niveles_cor_eda,
        xlab = "Péntada s", ylab = "Péntada t", main = "Correlación de satélite")
# La franja de correlación alta cerca de la diagonal también aparece en satélite
guardar_figura("05_exploratorio_correlacion.png", 12, 6)
par(par_cov_cor_eda)

# Covariación entre finales de abril y comienzos de noviembre
covariacion_temporadas_eda <- data.frame()
for (fuente in names(covarianza_eda)) {
  s <- match(24, pentadas_eda)
  t <- match(62, pentadas_eda)
  covariacion_temporadas_eda <- rbind(covariacion_temporadas_eda,
                                      data.frame(fuente = fuente, pentada_s = 24, pentada_t = 62,
                                                 covarianza_mm2 = covarianza_eda[[fuente]][s, t],
                                                 correlacion = correlacion_eda[[fuente]][s, t]))
}
covariacion_temporadas_eda
guardar_tabla(covariacion_temporadas_eda, "05_exploratorio_covariacion_abril_noviembre.csv")
# Entre esas péntadas de las dos temporadas lluviosas, la covarianza es 327,53 mm^2 en estaciones
# y 162,62 mm^2 en satélite; las correlaciones son 0,79 y 0,84, respectivamente
# Los pares estación-año con lluvia por encima de la media en abril tienden a estar también por encima
# de la media de noviembre; la covarianza mayor de estaciones incluye su mayor dispersión
# La correlación se calculó reuniendo estaciones y años; no muestra si la relación se mantiene
# dentro de cada estación o al separar los años por fase ENSO


# Covarianza y correlación cruzadas entre satélite y estaciones
sat_centrado_eda <- valores_eda$satelital - media_eda$satelital
est_centrado_eda <- valores_eda$terrestre - media_eda$terrestre
covarianza_cruzada_eda <- tcrossprod(sat_centrado_eda, est_centrado_eda) / (length(pares_eda) - 1)
correlacion_cruzada_eda <- covarianza_cruzada_eda /
  outer(desviacion_eda$satelital, desviacion_eda$terrestre)
diagonal_cruzada_eda <- diag(correlacion_cruzada_eda)

covarianza_cruzada_R <- eval.bifd(pentadas_eda, pentadas_eda,
                                 var.fd(curvas_eda_fd$satelital, curvas_eda_fd$terrestre))
correlacion_cruzada_R <- cor.fd(pentadas_eda, curvas_eda_fd$satelital,
                               pentadas_eda, curvas_eda_fd$terrestre)

# Diferencias máximas frente a fda y resumen de la diagonal
c(covarianza_mm2 = max(abs(covarianza_cruzada_eda - covarianza_cruzada_R)),
  correlacion = max(abs(correlacion_cruzada_eda - correlacion_cruzada_R)))
summary(diagonal_cruzada_eda)

par_cor_cruzada_eda <- par(mfrow = c(1, 1), bg = "white")
contour(pentadas_eda, pentadas_eda, covarianza_cruzada_eda,
        xlab = "Péntada s: satélite", ylab = "Péntada t: estaciones",
        main = "Covarianza cruzada (mm^2)")
guardar_figura("05_exploratorio_covarianza_cruzada.png", 7, 6)

par(mfrow = c(1, 2))
image(pentadas_eda, pentadas_eda, correlacion_cruzada_eda,
      zlim = c(-1, 1), col = hcl.colors(40, "Blue-Red 3"),
      xlab = "Péntada s: satélite", ylab = "Péntada t: estaciones",
      main = "Correlación cruzada")
contour(pentadas_eda, pentadas_eda, correlacion_cruzada_eda, levels = seq(-1, 1, 0.2), add = TRUE)
abline(0, 1, lty = 3, col = "red")
plot(pentadas_eda, diagonal_cruzada_eda, type = "l", lwd = 2, col = "darkblue", ylim = c(-1, 1),
     xlab = "Péntada", ylab = "Correlación", main = "Diagonal de la correlación cruzada")
abline(h = 0, lty = 3)
guardar_figura("05_exploratorio_correlacion_cruzada.png", 12, 5)
par(par_cor_cruzada_eda)
# En la misma péntada, la correlación entre ambas fuentes varía de 0,31 a 0,68
# Los pares con más lluvia en una fuente tienden a mostrar más lluvia en la otra


# Primera derivada de las curvas de precipitación por péntada
derivadas_eda <- list()
derivada_media_eda <- data.frame(pentada = pentadas_eda)
error_derivada_eda <- numeric()
for (fuente in names(curvas_eda_fd)) {
  derivadas_eda[[fuente]] <- eval.fd(pentadas_eda, curvas_eda_fd[[fuente]], Lfdobj = 1)
  derivada_media_eda[[fuente]] <- as.vector(eval.fd(pentadas_eda, media_eda_fd[[fuente]], Lfdobj = 1))
  error_derivada_eda[fuente] <- max(abs(derivada_media_eda[[fuente]] - rowMeans(derivadas_eda[[fuente]])))
}
# Diferencia máxima entre derivar la media y promediar las derivadas
error_derivada_eda

# Tramos delimitados por los cambios de signo de la primera derivada
tramos_derivada_eda <- data.frame()
for (fuente in names(curvas_eda_fd)) {
  curva_media <- media_eda_fd[[fuente]]
  pendiente_media <- function(t) as.vector(eval.fd(t, curva_media, Lfdobj = 1))
  pendiente <- derivada_media_eda[[fuente]]
  cambios <- which(diff(sign(pendiente)) != 0)
  giros <- vapply(cambios, function(i) uniroot(pendiente_media, pentadas_eda[c(i, i + 1)])$root,
                  numeric(1))
  limites <- sort(unique(c(min(pentadas_eda), giros, max(pentadas_eda))))
  inicio <- head(limites, -1)
  fin <- tail(limites, -1)
  tramos_derivada_eda <- rbind(tramos_derivada_eda,
                               data.frame(fuente = fuente, inicio_pentada = inicio, fin_pentada = fin,
                                          signo_pendiente = sign(pendiente_media((inicio + fin) / 2)),
                                          lluvia_inicio_mm = as.vector(eval.fd(inicio, curva_media)),
                                          lluvia_fin_mm = as.vector(eval.fd(fin, curva_media))))
}
tramos_derivada_eda
guardar_tabla(tramos_derivada_eda, "05_exploratorio_tramos_derivada.csv")

# En los tramos de descenso de mitad de año, la media baja de 30,07 a 12,28 mm en estaciones
# y de 37,45 a 12,13 mm en satélite
# El valle conserva lluvia; corresponde al periodo menos lluvioso entre las dos temporadas del ciclo bimodal

# Inicio y fin según los cruces con la lluvia media del año
umbral_lluvia_eda <- colMeans(media_eda[pentadas_eda %in% tiempos, names(curvas_eda_fd)])
cruces_lluvia_eda <- data.frame()
temporadas_lluvia_eda <- data.frame()
for (fuente in names(curvas_eda_fd)) {
  anomalia_eda <- media_eda[[fuente]] - umbral_lluvia_eda[fuente]
  cambios_eda <- diff(anomalia_eda > 0)
  indices_cruce <- which(cambios_eda != 0)
  if (!length(indices_cruce)) next
  
  # Interpolación entre los puntos que encierran cada cruce
  pentada_cruce <- pentadas_eda[indices_cruce] -
    anomalia_eda[indices_cruce] * diff(pentadas_eda)[indices_cruce] / diff(anomalia_eda)[indices_cruce]
  pendiente_cruce <- as.vector(eval.fd(pentada_cruce, media_eda_fd[[fuente]], Lfdobj = 1))
  tipo_cruce <- rep("Indeterminado", length(indices_cruce))
  tipo_cruce[cambios_eda[indices_cruce] > 0 & pendiente_cruce > 0] <- "Inicio"
  tipo_cruce[cambios_eda[indices_cruce] < 0 & pendiente_cruce < 0] <- "Fin"
  
  # Tramos con inicio y fin dentro del año
  inicios <- which(head(tipo_cruce, -1) == "Inicio" & tail(tipo_cruce, -1) == "Fin")
  temporada_cruce <- rep(NA_integer_, length(tipo_cruce))
  temporada_cruce[inicios] <- temporada_cruce[inicios + 1] <- seq_along(inicios)
  cruces_lluvia_eda <- rbind(cruces_lluvia_eda,
                             data.frame(fuente = fuente, pentada = pentada_cruce, tipo = tipo_cruce,
                                        pendiente = pendiente_cruce, temporada = temporada_cruce))
  temporadas_lluvia_eda <- rbind(temporadas_lluvia_eda,
                                 data.frame(fuente = rep(fuente, length(inicios)), temporada = seq_along(inicios),
                                            inicio_pentada = pentada_cruce[inicios], fin_pentada = pentada_cruce[inicios + 1]))
}
umbral_lluvia_eda
cruces_lluvia_eda
temporadas_lluvia_eda
guardar_tabla(cruces_lluvia_eda, "05_exploratorio_cruces_umbral_lluvia.csv")
guardar_tabla(temporadas_lluvia_eda, "05_exploratorio_temporadas_lluvia.csv")
# Con umbrales de 23,69 mm en estaciones y 25,08 mm en satélite, los dos tramos completos abarcan
# aproximadamente marzo-junio y septiembre-diciembre, con los picos de abril-mayo y noviembre dentro de ellos
# Son periodos por encima del promedio de cada fuente; sus límites dependen del umbral y del suavizado
# El cruce ascendente del satélite cerca de la péntada 71 queda sin cierre dentro del año
# Estas fechas corresponden a la media de la red

# Primera derivada de la media de estaciones y satélite
plot(pentadas_eda, derivada_media_eda$terrestre, type = "l", lwd = 2, col = "steelblue",
     ylim = range(0, derivada_media_eda[-1]), xlab = "Péntada",
     ylab = "Cambio (mm por péntada)", main = "Primera derivada de la lluvia media")
lines(pentadas_eda, derivada_media_eda$satelital, lwd = 2, lty = 2, col = "darkorange")
abline(h = 0, lty = 3)
legend("topleft", c("Estaciones", "Satélite"), col = c("steelblue", "darkorange"),
       lty = c(1, 2), lwd = 2, bty = "n")
guardar_figura("05_exploratorio_derivada_media.png", 8, 5)
# Por encima de cero, la lluvia media aumenta; por debajo, disminuye
# En estaciones, la subida de comienzos de año empieza entre la cuarta y la quinta péntada
# Después del primer pico, la lluvia baja hasta cerca de la péntada 44 en estaciones y de la 41 en satélite
# Las dos fuentes vuelven a subir hacia el pico de noviembre y después bajan
# Los dos tramos principales de aumento describen la entrada hacia los máximos del ciclo bimodal
# La pendiente negativa después de cada máximo indica que disminuye la cantidad de lluvia por péntada,
# aunque todavía puede estar por encima del umbral de la temporada lluviosa


# Diferencia funcional por par estación-año
diferencia_curvas_eda <- valores_eda$terrestre - valores_eda$satelital
diferencia_media_eda <- rowMeans(diferencia_curvas_eda)
ee_diferencia_eda <- apply(diferencia_curvas_eda, 1, sd) / sqrt(ncol(diferencia_curvas_eda))
inferior_diferencia_eda <- diferencia_media_eda - 2 * ee_diferencia_eda
superior_diferencia_eda <- diferencia_media_eda + 2 * ee_diferencia_eda

par_diferencia_eda <- par(mfrow = c(1, 2), bg = "white")
matplot(pentadas_eda, diferencia_curvas_eda, type = "l", lty = 1,
        col = adjustcolor("grey40", alpha.f = 0.12),
        xlab = "Péntada", ylab = "Estación - satélite (mm)", main = "Diferencias por par")
lines(pentadas_eda, diferencia_media_eda, lwd = 2, col = "darkred")
abline(h = 0, lty = 3)
plot(pentadas_eda, diferencia_media_eda, type = "l", lwd = 2, col = "darkred",
     ylim = range(0, inferior_diferencia_eda, superior_diferencia_eda),
     xlab = "Péntada", ylab = "Estación - satélite (mm)", main = "Diferencia media +- 2 EE")
lines(pentadas_eda, inferior_diferencia_eda, lty = 2, col = "darkred")
lines(pentadas_eda, superior_diferencia_eda, lty = 2, col = "darkred")
abline(h = 0, lty = 3)
guardar_figura("05_exploratorio_diferencia_funcional.png", 12, 5)
par(par_diferencia_eda)

# Tramos del año donde la banda +- 2 EE no contiene el cero
signo_banda <- ifelse(inferior_diferencia_eda > 0, 1, ifelse(superior_diferencia_eda < 0, -1, 0))
tramos <- rle(signo_banda)
fin_tramo <- cumsum(tramos$lengths)
inicio_tramo <- fin_tramo - tramos$lengths + 1
intervalos_diferencia_eda <- data.frame()
for (k in which(tramos$values != 0)) {
  idx <- inicio_tramo[k]:fin_tramo[k]
  intervalos_diferencia_eda <- rbind(intervalos_diferencia_eda, data.frame(
    inicio_pentada = pentadas_eda[inicio_tramo[k]], fin_pentada = pentadas_eda[fin_tramo[k]],
    signo = ifelse(tramos$values[k] > 0, "Estacion mayor", "Satelite mayor"),
    diferencia_min_mm = min(diferencia_media_eda[idx]), diferencia_max_mm = max(diferencia_media_eda[idx])))
}
intervalos_diferencia_eda
guardar_tabla(intervalos_diferencia_eda, "05_exploratorio_intervalos_diferencia.csv")
# La diferencia media es negativa cerca de ambos picos lluviosos
# La banda de +- 2 EE es descriptiva porque hay estaciones repetidas y píxeles compartidos

# Profundidad de Fraiman-Muniz y boxplot funcional con MBD
profundidades_eda <- data.frame()
boxplot_eda <- list()
limites_boxplot_eda <- range(unlist(valores_eda))
par_boxplot_eda <- par(mfrow = c(1, 2))
for (fuente in names(valores_eda)) {
  valores <- valores_eda[[fuente]]
  datos_fdata <- fda.usc::fdata(t(valores), pentadas_eda, range(pentadas_eda))
  profundidad_FM <- fda.usc::depth.FM(datos_fdata)$dep
  
  # Las bandas cerradas incluyen las curvas que coinciden con sus límites
  n_curvas <- ncol(valores)
  n_bandas <- choose(n_curvas, 2)
  bandas_que_contienen <- apply(valores, 1, function(lluvia) {
    debajo <- rank(lluvia, ties.method = "min") - 1
    encima <- n_curvas - rank(lluvia, ties.method = "max")
    n_bandas - choose(debajo, 2) - choose(encima, 2)
  })
  profundidad_MBD <- rowMeans(bandas_que_contienen) / n_bandas
  
  boxplot_eda[[fuente]] <- fda::fbplot(valores, pentadas_eda, method = "MBD",
                                       depth = profundidad_MBD, prob = 0.5, factor = 1.5,
                                       xlab = "Péntada", ylab = "Precipitación (mm)",
                                       xlim = range(pentadas_eda), ylim = limites_boxplot_eda)
  title(c(terrestre = "Estaciones", satelital = "Satélite")[fuente])
  
  indices <- match(colnames(valores), metadatos_eda$curva)
  profundidades_eda <- rbind(profundidades_eda,
                             data.frame(fuente = fuente, metadatos_eda[indices, c("STATION_ID", "YEAR", "pixel", "curva")],
                                        FM = unname(profundidad_FM), MBD = profundidad_MBD,
                                        atipica = seq_len(n_curvas) %in% boxplot_eda[[fuente]]$outpoint))
}
guardar_figura("05_exploratorio_boxplot_funcional.png", 12, 5)
par(par_boxplot_eda)
# La caja contiene el 50 % central de las curvas; las líneas rojas discontinuas son las señaladas como atípicas

# Identificación de curvas atípicas por estación y año
profundidades_eda <- as_tibble(profundidades_eda)
profundidades_eda
guardar_tabla(profundidades_eda, "05_exploratorio_profundidades.csv")
# En cada fuente, una profundidad mayor indica una curva más central respecto al conjunto de estaciones y años
# FM y MBD resumen la posición de toda la curva; un pico en noviembre puede ser parte del patrón habitual
# La profundidad no identifica por sí misma una fase ENSO ni una condición de relieve
atipicas_eda <- profundidades_eda %>%
  filter(atipica) %>%
  arrange(fuente, STATION_ID, YEAR) %>%
  select(-atipica) %>%
  left_join(metadatos_eda[, c("curva", "altitud_completa", "fuente_altitud")], by = "curva") %>%
  rename(altitud_m = altitud_completa) %>%
  left_join(enso_anual[, c("YEAR", "secuencia_fases")], by = "YEAR") %>%
  as.data.frame()
atipicas_eda
guardar_tabla(atipicas_eda, "05_exploratorio_curvas_atipicas.csv")
# Hay casos señalados en 1987 (El Niño), 1990 (neutral) y 1999 (La Niña), según la clasificación ONI previa
# También aparecen en años de transición, como 1988; asignar una sola fase a todo ese año ocultaría el cambio
# Estos casos muestran que las curvas atípicas no aparecen exclusivamente en una fase ENSO
# La secuencia anual no identifica en qué fase ocurrió el tramo de la curva que salió de la banda
resumen_atipicas_eda <- profundidades_eda %>%
  group_by(fuente) %>%
  summarize(curvas = n(), atipicas = sum(atipica), porcentaje = 100 * mean(atipica), .groups = "drop")
resumen_atipicas_eda
guardar_tabla(resumen_atipicas_eda, "05_exploratorio_resumen_atipicas.csv")
# El boxplot señala 64 pares estación-año en estaciones y 28 en satélite
# Representan el 2,04 % y el 0,89 % de las curvas de cada fuente, respectivamente

atipicas_por_estacion_eda <- profundidades_eda %>%
  group_by(fuente, STATION_ID) %>%
  summarize(curvas = n(), atipicas = sum(atipica), .groups = "drop") %>%
  filter(atipicas > 0) %>%
  left_join(unique(metadatos_eda[, c("STATION_ID", "altitud_completa", "fuente_altitud")]), by = "STATION_ID") %>%
  rename(altitud_m = altitud_completa) %>%
  arrange(fuente, desc(atipicas))
atipicas_por_estacion_eda
guardar_tabla(atipicas_por_estacion_eda, "05_exploratorio_atipicas_por_estacion.csv")
# En tierra, la estación 5301500106, a 232 m, aparece señalada en 17 de sus 19 años disponibles
# La estación 5301000104, a 20 m, aparece en sus 11 años en tierra y en 10 de los 11 en satélite
# La repetición muestra que el boxplot detecta diferencias que persisten de algunos sitios frente a la red,
# además de posibles años excepcionales; la marca no significa necesariamente un año anómalo para ese sitio
# Hay casos terrestres desde 20 m hasta 2100 m; algunas estaciones tienen la altitud tomada del modelo digital de elevación


#################### Comparación funcional por fuente ####################
##### Colapso a píxel-año único y suavizado simétrico #####
# Dispersión entre estaciones que comparten un mismo píxel CHIRPS y año
clave_pixel_anio <- paste(metadatos_eda$pixel, metadatos_eda$YEAR, sep = "_")
estaciones_unicas_eda <- unique(metadatos_eda[, c("STATION_ID", "pixel", "altitud_completa", "altitud_DEM")])
conteo_estaciones_pixel <- table(estaciones_unicas_eda$pixel)
pixeles_compartidos_eda <- sort(as.numeric(names(conteo_estaciones_pixel[conteo_estaciones_pixel > 1])))

dispersion_intrapixel <- do.call(rbind, lapply(pixeles_compartidos_eda, function(px) {
  est_px <- estaciones_unicas_eda[estaciones_unicas_eda$pixel == px, ]
  claves_px <- unique(clave_pixel_anio[metadatos_eda$pixel == px])
  claves_activas <- claves_px[sapply(claves_px, function(cl) sum(clave_pixel_anio == cl) >= 2)]
  de_pentadas <- numeric()
  rangos_anuales <- numeric()
  for (cl in claves_activas) {
    filas <- which(clave_pixel_anio == cl)
    matriz_px <- fuentes_eda$terrestre[filas, , drop = FALSE]
    de_pentadas <- c(de_pentadas, apply(matriz_px, 2, sd, na.rm = TRUE))
    totales_anuales <- rowSums(matriz_px, na.rm = TRUE)
    rangos_anuales <- c(rangos_anuales, max(totales_anuales) - min(totales_anuales))
  }
  data.frame(
    pixel = px,
    estaciones = nrow(est_px),
    pares_pixel_anio = length(claves_activas),
    altitud_min_m = min(est_px$altitud_completa, na.rm = TRUE),
    altitud_max_m = max(est_px$altitud_completa, na.rm = TRUE),
    desnivel_m = max(est_px$altitud_completa, na.rm = TRUE) - min(est_px$altitud_completa, na.rm = TRUE),
    DE_intrapixel_media_mm = mean(de_pentadas, na.rm = TRUE),
    DE_intrapixel_max_mm = max(de_pentadas, na.rm = TRUE),
    rango_anual_medio_mm = mean(rangos_anuales, na.rm = TRUE)
  )
}))
dispersion_intrapixel
guardar_tabla(dispersion_intrapixel, "06_comparacion_dispersion_intrapixel.csv")
# En los 24 píxeles compartidos por 53 estaciones (362 pares píxel-año con dos o más estaciones activas),
# el desnivel interno varía de 5 a 827 m (mediana 103,5 m) y la DE intra-píxel media por péntada va de
# 2,96 a 35,44 mm (mediana 8,44 mm; media 9,54 mm), cuantificando la diferencia de soporte entre
# pluviómetros puntuales y el promedio areal del píxel de 0,05 grados

# Matrices colapsadas a nivel de píxel-año único con las mismas péntadas observadas en ambas fuentes
claves_pixel_anio <- unique(clave_pixel_anio)
indices_por_clave <- split(seq_len(nrow(metadatos_eda)), factor(clave_pixel_anio, levels = claves_pixel_anio))
metadatos_comparacion <- do.call(rbind, lapply(seq_along(claves_pixel_anio), function(k) {
  idx <- indices_por_clave[[k]]
  sub <- metadatos_eda[idx, , drop = FALSE]
  part_sub <- unique(particion_curva(sub$STATION_ID, sub$YEAR))
  stopifnot(length(part_sub) == 1)
  data.frame(
    clave_pixel_anio = claves_pixel_anio[k],
    pixel = sub$pixel[1],
    YEAR = sub$YEAR[1],
    estaciones_en_pixel = nrow(sub),
    longitud = mean(sub$longitud),
    latitud = mean(sub$latitud),
    altitud_completa = if (all(is.na(sub$altitud_completa))) NA_real_ else mean(sub$altitud_completa, na.rm = TRUE),
    altitud_DEM = mean(sub$altitud_DEM, na.rm = TRUE),
    particion = part_sub
  )
}))
rownames(metadatos_comparacion) <- NULL

Y_pixel_anio <- t(sapply(indices_por_clave, function(idx) {
  if (length(idx) == 1) fuentes_eda$terrestre[idx, ] else colMeans(fuentes_eda$terrestre[idx, , drop = FALSE], na.rm = TRUE)
}))
Y_pixel_anio[is.nan(Y_pixel_anio)] <- NA_real_
X_pixel_anio <- t(sapply(indices_por_clave, function(idx) {
  if (length(idx) == 1) fuentes_eda$satelital[idx, ] else colMeans(fuentes_eda$satelital[idx, , drop = FALSE], na.rm = TRUE)
}))
X_pixel_anio[is.na(Y_pixel_anio)] <- NA_real_
rownames(Y_pixel_anio) <- metadatos_comparacion$clave_pixel_anio
rownames(X_pixel_anio) <- metadatos_comparacion$clave_pixel_anio
colnames(Y_pixel_anio) <- colnames(fuentes_eda$terrestre)
colnames(X_pixel_anio) <- colnames(fuentes_eda$satelital)
D_pixel_anio <- Y_pixel_anio - X_pixel_anio
c(pares_pixel_anio = nrow(metadatos_comparacion), pixeles_unicos = n_distinct(metadatos_comparacion$pixel))
table(metadatos_comparacion$particion)
# El colapso reduce los 3139 pares estación-año (176 estaciones) a 2721 pares píxel-año en 147 píxeles únicos
# (1690 en entrenamiento, 110 en prueba y 921 fuera de la partición del modelo)

# Suavizado simétrico con la misma base B-spline cúbica y lambda común elegido por VCG en entrenamiento
suavizar_comparacion <- function(datos, base, lambda) {
  fd_par <- fdPar(base, int2Lfd(2), lambda)
  completas <- which(rowSums(is.na(datos)) == 0)
  incompletas <- which(rowSums(is.na(datos)) > 0)
  coeficientes <- matrix(NA_real_, base$nbasis, nrow(datos),
                         dimnames = list(NULL, rownames(datos)))
  vcg <- numeric(nrow(datos))
  gl <- numeric(nrow(datos))
  if (length(completas) > 0) {
    ajuste <- smooth.basis(tiempos, t(datos[completas, , drop = FALSE]), fd_par)
    coeficientes[, completas] <- ajuste$fd$coefs
    vcg[completas] <- ajuste$gcv
    gl[completas] <- ajuste$df
  }
  for (i in incompletas) {
    obs <- which(!is.na(datos[i, ]))
    ajuste <- smooth.basis(tiempos[obs], datos[i, obs], fd_par)
    coeficientes[, i] <- ajuste$fd$coefs[, 1]
    vcg[i] <- ajuste$gcv
    gl[i] <- ajuste$df
  }
  list(fd = fd(coeficientes, base, list("Péntada", rownames(datos), "Precipitación (mm)")),
       gcv = sum(vcg), df = mean(gl))
}

en_entrenamiento_comp <- metadatos_comparacion$particion == "entrenamiento"
log_lambda_comp <- seq(-1, 5, by = 0.25)
tabla_vcg_comparacion <- do.call(rbind, lapply(log_lambda_comp, function(ll) {
  lam <- 10^ll
  aj_est <- suavizar_comparacion(Y_pixel_anio[en_entrenamiento_comp, , drop = FALSE], base_bspline, lam)
  aj_sat <- suavizar_comparacion(X_pixel_anio[en_entrenamiento_comp, , drop = FALSE], base_bspline, lam)
  aj_dif <- suavizar_comparacion(D_pixel_anio[en_entrenamiento_comp, , drop = FALSE], base_bspline, lam)
  data.frame(
    log10_lambda = ll,
    lambda = lam,
    gl_efectivos = aj_est$df,
    VCG_estaciones = aj_est$gcv,
    VCG_satelite = aj_sat$gcv,
    VCG_total = aj_est$gcv + aj_sat$gcv,
    VCG_diferencia = aj_dif$gcv
  )
}))
tabla_vcg_comparacion
guardar_tabla(tabla_vcg_comparacion, "06_comparacion_VCG_suavizado_simetrico.csv")

idx_opt_comp <- which.min(tabla_vcg_comparacion$VCG_total)
log_lambda_opt <- tabla_vcg_comparacion$log10_lambda[idx_opt_comp]
lambda_opt <- tabla_vcg_comparacion$lambda[idx_opt_comp]
gl_opt <- tabla_vcg_comparacion$gl_efectivos[idx_opt_comp]
tabla_vcg_comparacion[idx_opt_comp, ]

png("Figuras/06_comparacion_VCG_lambda_comun.png", width = 6.5, height = 4.5, units = "in", res = 150)
plot(tabla_vcg_comparacion$log10_lambda, tabla_vcg_comparacion$VCG_total, type = "b", pch = 19,
     xlab = "log10(lambda)", ylab = "VCG total (estaciones + satélite)",
     main = "Selección de lambda común para comparación funcional")
abline(v = log_lambda_opt, col = "red", lty = 2, lwd = 2)
dev.off()
# La VCG conjunta sobre los 1690 pares píxel-año de entrenamiento selecciona log10(lambda) = 2,25
# (lambda = 177,83; 7,96 grados de libertad efectivos) para suavizar ambas fuentes con el mismo operador

ajuste_comp_est <- suavizar_comparacion(Y_pixel_anio, base_bspline, lambda_opt)
ajuste_comp_sat <- suavizar_comparacion(X_pixel_anio, base_bspline, lambda_opt)
ajuste_comp_dif <- suavizar_comparacion(D_pixel_anio, base_bspline, lambda_opt)
valores_comparacion <- list(
  terrestre = eval.fd(pentadas_eda, ajuste_comp_est$fd),
  satelital = eval.fd(pentadas_eda, ajuste_comp_sat$fd)
)
diferencia_comparacion <- eval.fd(pentadas_eda, ajuste_comp_dif$fd)
error_linealidad_comp <- max(abs((valores_comparacion$terrestre - valores_comparacion$satelital) - diferencia_comparacion))
error_linealidad_comp
stopifnot(error_linealidad_comp < 1e-11)
# La diferencia máxima entre restar las curvas suavizadas con el mismo operador y suavizar directamente
# las diferencias discretas es 2,56e-13 mm en los 143 puntos y 2721 pares píxel-año

# Sensibilidad a lambda (lambda/10, lambda óptimo y 10*lambda) frente a datos discretos y suavizado asimétrico
if (!exists("diferencia_curvas_eda")) {
  diferencia_curvas_eda <- valores_eda$terrestre - valores_eda$satelital
}
media_discreta_comp <- colMeans(D_pixel_anio, na.rm = TRUE)
idx_72_comp <- match(tiempos, pentadas_eda)
diferencia_asimetrica_pa <- t(sapply(indices_por_clave, function(idx) {
  if (length(idx) == 1) diferencia_curvas_eda[idx_72_comp, idx] else rowMeans(diferencia_curvas_eda[idx_72_comp, idx, drop = FALSE])
}))
media_asimetrica_pa <- colMeans(diferencia_asimetrica_pa)
media_asimetrica_pares <- rowMeans(diferencia_curvas_eda[idx_72_comp, ])
aj_dif_bajo <- suavizar_comparacion(D_pixel_anio, base_bspline, lambda_opt / 10)
aj_dif_alto <- suavizar_comparacion(D_pixel_anio, base_bspline, lambda_opt * 10)
media_sim_bajo <- rowMeans(eval.fd(tiempos, aj_dif_bajo$fd))
media_sim_opt <- rowMeans(diferencia_comparacion[idx_72_comp, ])
media_sim_alto <- rowMeans(eval.fd(tiempos, aj_dif_alto$fd))

evaluar_sensibilidad_comp <- function(nombre, log_lam, gl_ef, media_curva) {
  cambio <- media_curva - media_discreta_comp
  mismo_signo <- sign(media_curva) == sign(media_discreta_comp)
  data.frame(
    esquema = nombre,
    log10_lambda = log_lam,
    gl_efectivos = gl_ef,
    MAE_total_mm = mean(abs(cambio)),
    MAE_inicio_1_6_mm = mean(abs(cambio[1:6])),
    MAE_interior_7_66_mm = mean(abs(cambio[7:66])),
    MAE_final_67_72_mm = mean(abs(cambio[67:72])),
    concordancia_signo_pct = mean(mismo_signo) * 100,
    media_P27_mm = media_curva[27],
    media_P68_mm = media_curva[68],
    media_P72_mm = media_curva[72]
  )
}

sensibilidad_suavizado_comp <- rbind(
  evaluar_sensibilidad_comp("Discreta pixel-anio", NA_real_, 72, media_discreta_comp),
  evaluar_sensibilidad_comp("Asimetrico K6 vs K11 (estacion-anio)", NA_real_, NA_real_, media_asimetrica_pares),
  evaluar_sensibilidad_comp("Asimetrico K6 vs K11 (pixel-anio)", NA_real_, NA_real_, media_asimetrica_pa),
  evaluar_sensibilidad_comp("Simetrico lambda/10", log_lambda_opt - 1, aj_dif_bajo$df, media_sim_bajo),
  evaluar_sensibilidad_comp("Simetrico lambda opt", log_lambda_opt, gl_opt, media_sim_opt),
  evaluar_sensibilidad_comp("Simetrico 10*lambda", log_lambda_opt + 1, aj_dif_alto$df, media_sim_alto)
)
rownames(sensibilidad_suavizado_comp) <- NULL
sensibilidad_suavizado_comp
guardar_tabla(sensibilidad_suavizado_comp, "06_comparacion_sensibilidad_suavizado_simetrico.csv")

png("Figuras/06_comparacion_sensibilidad_suavizado_simetrico.png", width = 9.5, height = 5, units = "in", res = 150)
plot(tiempos, media_discreta_comp, type = "l", lwd = 1.5, col = "grey45",
     ylim = range(media_discreta_comp, media_asimetrica_pa, media_sim_bajo, media_sim_opt, media_sim_alto),
     xlab = "Péntada", ylab = "Estación - satélite (mm)",
     main = "Diferencia media: suavizado asimétrico frente a suavizado simétrico")
lines(tiempos, media_asimetrica_pa, lwd = 2, lty = 2, col = "firebrick")
lines(tiempos, media_sim_bajo, lwd = 1.8, lty = 3, col = "steelblue")
lines(tiempos, media_sim_opt, lwd = 2.4, lty = 1, col = "navy")
lines(tiempos, media_sim_alto, lwd = 1.8, lty = 4, col = "darkcyan")
abline(h = 0, lty = 3)
legend("bottomleft",
       c("Discreta píxel-año", "Asimétrico K=6 vs K=11", "Simétrico lambda/10 (gl=13,35)",
         "Simétrico lambda opt (gl=7,96)", "Simétrico 10*lambda (gl=4,92)"),
       col = c("grey45", "firebrick", "steelblue", "navy", "darkcyan"),
       lty = c(1, 2, 3, 1, 4), lwd = c(1.5, 2, 1.8, 2.4, 1.8), bty = "n", cex = 0.78)
dev.off()
# Con suavizado simétrico desaparecen la inversión de signo en la péntada 68 (-4,56 mm discreto ->
# -1,35 mm simétrico frente a +6,48 mm asimétrico) y el artefacto de borde en la péntada 72
# (-5,34 mm discreto -> -3,05 mm simétrico frente a -20,64 mm asimétrico), reduciendo el error medio
# en las péntadas 67-72 de 6,33 mm a 3,02 mm y subiendo la concordancia de signo de 65,3 % a 73,6 %

##### Comparación funcional entre estación y satélite #####
# Funciones auxiliares para intervalos y pruebas por permutación
B_comparacion <- 5000
alfa_comparacion <- 0.05
p_maximo <- function(maximos_nulos, observado) {
  (1 + sum(maximos_nulos >= observado - 1e-12)) / (length(maximos_nulos) + 1)
}
intervalos_pentadas <- function(indicador, malla = pentadas_eda) {
  r <- rle(indicador)
  if (!any(r$values)) {
    return(data.frame(inicio_pentada = numeric(0), fin_pentada = numeric(0), puntos = integer(0)))
  }
  fin <- cumsum(r$lengths)
  ini <- fin - r$lengths + 1
  data.frame(inicio_pentada = malla[ini[r$values]],
             fin_pentada = malla[fin[r$values]],
             puntos = r$lengths[r$values])
}
stopifnot(identical(colnames(diferencia_comparacion), as.character(metadatos_comparacion$clave_pixel_anio)),
          !anyDuplicated(metadatos_comparacion[, c("pixel", "YEAR")]))
diferencia_por_anio <- function(incluir, anios) {
  sapply(anios, function(anio) {
    pares <- which(incluir & metadatos_comparacion$YEAR == anio)
    rowMeans(diferencia_comparacion[, pares, drop = FALSE])
  })
}
anios_red <- sort(unique(metadatos_comparacion$YEAR))
D_red <- diferencia_por_anio(rep(TRUE, nrow(metadatos_comparacion)), anios_red)
n_red <- length(anios_red)

# Curvas climatológicas medias por píxel único (P = 147) de cada fuente por separado
pixeles_unicos <- sort(unique(metadatos_comparacion$pixel))
indices_por_pixel <- split(seq_len(nrow(metadatos_comparacion)),
                           factor(metadatos_comparacion$pixel, levels = pixeles_unicos))
metadatos_pixeles <- do.call(rbind, lapply(seq_along(pixeles_unicos), function(k) {
  idx <- indices_por_pixel[[k]]
  sub <- metadatos_comparacion[idx, , drop = FALSE]
  alt_comp <- if (all(is.na(sub$altitud_completa))) NA_real_ else mean(sub$altitud_completa, na.rm = TRUE)
  alt_dem <- mean(sub$altitud_DEM, na.rm = TRUE)
  alt_ref <- ifelse(!is.na(alt_comp), alt_comp, alt_dem)
  piso <- case_when(
    alt_ref < 700 ~ "Pacífico bajo (<700 m)",
    alt_ref < 1240 ~ "Valle (700-1239 m)",
    alt_ref <= 2000 ~ "Montaña media (1240-2000 m)",
    TRUE ~ "Alta montaña (>2000 m)"
  )
  data.frame(
    pixel = pixeles_unicos[k],
    anios_activos = nrow(sub),
    estaciones_max = max(sub$estaciones_en_pixel),
    longitud = mean(sub$longitud),
    latitud = mean(sub$latitud),
    altitud_m = alt_ref,
    piso_altitudinal = piso
  )
}))
rownames(metadatos_pixeles) <- NULL

X1_pixeles <- sapply(indices_por_pixel, function(idx) rowMeans(valores_comparacion$terrestre[, idx, drop = FALSE]))
X2_pixeles <- sapply(indices_por_pixel, function(idx) rowMeans(valores_comparacion$satelital[, idx, drop = FALSE]))
coefs1_pixeles <- sapply(indices_por_pixel, function(idx) rowMeans(ajuste_comp_est$fd$coefs[, idx, drop = FALSE]))
coefs2_pixeles <- sapply(indices_por_pixel, function(idx) rowMeans(ajuste_comp_sat$fd$coefs[, idx, drop = FALSE]))
colnames(X1_pixeles) <- colnames(X2_pixeles) <- as.character(pixeles_unicos)
colnames(coefs1_pixeles) <- colnames(coefs2_pixeles) <- as.character(pixeles_unicos)
table(metadatos_pixeles$piso_altitudinal)

# Estadístico t con varianza combinada y estadístico de tipo Welch T(t) para dos muestras
Estadisticos_t <- function(A, B) {
  nA <- ncol(A)
  nB <- ncol(B)
  mA <- rowMeans(A)
  mB <- rowMeans(B)
  vA <- rowSums((A - mA)^2) / (nA - 1)
  vB <- rowSums((B - mB)^2) / (nB - 1)
  sp2 <- ((nA - 1) * vA + (nB - 1) * vB) / (nA + nB - 2)
  list(t_comb = (mA - mB) / sqrt(sp2 * (1 / nA + 1 / nB)),
       T_welch = abs(mA - mB) / sqrt(vA / nA + vB / nB),
       dif = mA - mB,
       ee_welch = sqrt(vA / nA + vB / nB))
}

# Prueba t funcional de dos muestras por permutación de etiquetas con criterio del máximo
prueba_t_dos_muestras <- function(A, B_mat, nombre_grupo) {
  n1 <- ncol(A)
  n2 <- ncol(B_mat)
  n <- n1 + n2
  Est_obs <- Estadisticos_t(A, B_mat)
  X_todas <- cbind(A, B_mat)
  set.seed(2026)
  Grupos1 <- replicate(B_comparacion, sort(sample(n, n1)))
  T_nulos <- matrix(NA_real_, nrow(A), B_comparacion)
  for (b in seq_len(B_comparacion)) {
    g1 <- Grupos1[, b]
    T_nulos[, b] <- Estadisticos_t(X_todas[, g1, drop = FALSE], X_todas[, -g1, drop = FALSE])$T_welch
  }
  T_max_nulos <- apply(T_nulos, 2, max)
  T_obs <- Est_obs$T_welch
  c_puntual <- apply(T_nulos, 1, quantile, probs = 1 - alfa_comparacion)
  c_maximo <- unname(quantile(T_max_nulos, 1 - alfa_comparacion))
  p_puntual <- (1 + rowSums(T_nulos >= T_obs - 1e-12)) / (B_comparacion + 1)
  pico <- which.max(T_obs)
  p_global <- p_maximo(T_max_nulos, T_obs[pico])

  por_pentada <- data.frame(
    grupo = nombre_grupo,
    pentada = pentadas_eda,
    media_estaciones_mm = rowMeans(A),
    media_satelite_mm = rowMeans(B_mat),
    DE_estaciones_mm = apply(A, 1, sd),
    DE_satelite_mm = apply(B_mat, 1, sd),
    diferencia_mm = Est_obs$dif,
    EE_welch_mm = Est_obs$ee_welch,
    t_comb = Est_obs$t_comb,
    T_welch = T_obs,
    c_puntual = c_puntual,
    c_maximo = c_maximo,
    p_puntual = p_puntual,
    significativa = T_obs > c_maximo
  )

  iv <- intervalos_pentadas(T_obs > c_maximo)
  if (nrow(iv) > 0) {
    iv$grupo <- nombre_grupo
    for (k in seq_len(nrow(iv))) {
      idx <- which(pentadas_eda >= iv$inicio_pentada[k] & pentadas_eda <= iv$fin_pentada[k])
      j <- idx[which.max(T_obs[idx])]
      iv$dif_min_mm[k] <- min(Est_obs$dif[idx])
      iv$dif_max_mm[k] <- max(Est_obs$dif[idx])
      iv$pico_pentada[k] <- pentadas_eda[j]
      iv$pico_T_welch[k] <- T_obs[j]
    }
    iv <- iv[, c("grupo", "inicio_pentada", "fin_pentada", "puntos",
                 "dif_min_mm", "dif_max_mm", "pico_pentada", "pico_T_welch")]
  }

  resumen <- data.frame(
    grupo = nombre_grupo,
    n_estaciones = n1,
    n_satelite = n2,
    media_estaciones_mm = mean(A),
    media_satelite_mm = mean(B_mat),
    diferencia_media_mm = mean(Est_obs$dif),
    DE_media_estaciones_mm = mean(apply(A, 1, sd)),
    DE_media_satelite_mm = mean(apply(B_mat, 1, sd)),
    pentada_pico = pentadas_eda[pico],
    diferencia_pico_mm = Est_obs$dif[pico],
    EE_welch_pico_mm = Est_obs$ee_welch[pico],
    T_welch_max = T_obs[pico],
    c_maximo = c_maximo,
    p_global = p_global,
    puntos_sig_maximo = sum(T_obs > c_maximo),
    puntos_sig_puntual = sum(p_puntual < alfa_comparacion)
  )
  list(resumen = resumen, por_pentada = por_pentada, intervalos = iv)
}

# Evaluación en la red completa, cobertura >= 18 años, dominio andino >= 700 m y pisos altitudinales
subconjuntos_t <- list(
  "Red completa (P=147)" = rep(TRUE, nrow(metadatos_pixeles)),
  "Cobertura >=18 años (P=105)" = metadatos_pixeles$anios_activos >= 18,
  "Dominio andino >=700 m (P=140)" = metadatos_pixeles$altitud_m >= 700,
  "Pacífico bajo (<700 m)" = metadatos_pixeles$piso_altitudinal == "Pacífico bajo (<700 m)",
  "Valle (700-1239 m)" = metadatos_pixeles$piso_altitudinal == "Valle (700-1239 m)",
  "Montaña media (1240-2000 m)" = metadatos_pixeles$piso_altitudinal == "Montaña media (1240-2000 m)",
  "Alta montaña (>2000 m)" = metadatos_pixeles$piso_altitudinal == "Alta montaña (>2000 m)"
)
resultados_t <- lapply(names(subconjuntos_t), function(nom) {
  sel <- subconjuntos_t[[nom]]
  prueba_t_dos_muestras(X1_pixeles[, sel, drop = FALSE], X2_pixeles[, sel, drop = FALSE], nom)
})
names(resultados_t) <- names(subconjuntos_t)

t_resumen <- do.call(rbind, lapply(resultados_t, `[[`, "resumen"))
rownames(t_resumen) <- NULL
t_resumen
guardar_tabla(t_resumen, "06_comparacion_t_funcional_resumen.csv")

t_por_pentada <- do.call(rbind, lapply(resultados_t, `[[`, "por_pentada"))
rownames(t_por_pentada) <- NULL
guardar_tabla(t_por_pentada, "06_comparacion_t_por_pentada.csv")

intervalos_t <- do.call(rbind, Filter(function(x) nrow(x) > 0, lapply(resultados_t, `[[`, "intervalos")))
rownames(intervalos_t) <- NULL
intervalos_t
guardar_tabla(intervalos_t, "06_comparacion_t_intervalos.csv")

# Comprobación con tperm.fd de la librería fda en la red completa
set.seed(2026)
fd_est_pixeles <- fd(coefs1_pixeles, base_bspline)
fd_sat_pixeles <- fd(coefs2_pixeles, base_bspline)
tperm_red <- tperm.fd(fd_est_pixeles, fd_sat_pixeles, nperm = 1000, q = alfa_comparacion,
                      argvals = pentadas_eda, plotres = FALSE)
verificacion_tperm <- data.frame(
  diferencia_max_T_welch = max(abs(resultados_t[["Red completa (P=147)"]]$por_pentada$T_welch - tperm_red$Tvals)),
  c_maximo_manual = resultados_t[["Red completa (P=147)"]]$resumen$c_maximo,
  c_maximo_tperm = unname(tperm_red$qval),
  p_global_manual = resultados_t[["Red completa (P=147)"]]$resumen$p_global,
  p_global_tperm = tperm_red$pval
)
verificacion_tperm
stopifnot(verificacion_tperm$diferencia_max_T_welch < 1e-12)
# El estadístico T(t) de Welch manual reproduce tperm.fd() con diferencia máxima de 4,33e-15

# Figura 1: Media, DE y prueba t funcional por permutación en la red completa y dominio andino
p_red <- resultados_t[["Red completa (P=147)"]]$por_pentada
p_and <- resultados_t[["Dominio andino >=700 m (P=140)"]]$por_pentada
png("Figuras/06_comparacion_t_funcional_red.png", width = 11, height = 8, units = "in", res = 150)
par(mfrow = c(2, 2))
matplot(pentadas_eda, cbind(p_red$media_estaciones_mm, p_red$media_satelite_mm), type = "l", lty = c(1, 2),
        lwd = 2, col = c("steelblue", "darkorange"), xlab = "Péntada", ylab = "Precipitación (mm)",
        main = "Media funcional por píxel - Red completa (P=147)")
legend("topleft", c("Estaciones", "Satélite"), col = c("steelblue", "darkorange"), lty = c(1, 2), lwd = 2, bty = "n", cex = 0.8)
matplot(pentadas_eda, cbind(p_red$DE_estaciones_mm, p_red$DE_satelite_mm), type = "l", lty = c(1, 2),
        lwd = 2, col = c("steelblue", "darkorange"), xlab = "Péntada", ylab = "Desviación estándar (mm)",
        main = "DE funcional por píxel - Red completa (P=147)")
legend("topleft", c("Estaciones", "Satélite"), col = c("steelblue", "darkorange"), lty = c(1, 2), lwd = 2, bty = "n", cex = 0.8)
plot(pentadas_eda, p_red$T_welch, type = "l", lwd = 2, col = "red", xlab = "Péntada", ylab = "T(t) de Welch",
     main = "Prueba t funcional - Red completa (P=147)",
     ylim = c(0, max(p_red$T_welch, p_red$c_maximo) * 1.35))
lines(pentadas_eda, p_red$c_puntual, col = "blue", lty = 3, lwd = 2)
abline(h = p_red$c_maximo[1], col = "black", lty = 2)
legend("topright", c("T(t) observado", "Crítico puntual (0.05)", "Crítico máximo"),
       col = c("red", "blue", "black"), lty = c(1, 3, 2), lwd = c(2, 2, 1), bty = "n", cex = 0.75)
plot(pentadas_eda, p_and$T_welch, type = "l", lwd = 2, col = "red", xlab = "Péntada", ylab = "T(t) de Welch",
     main = "Prueba t funcional - Dominio andino >=700 m (P=140)",
     ylim = c(0, max(p_and$T_welch, p_and$c_maximo) * 1.35))
lines(pentadas_eda, p_and$c_puntual, col = "blue", lty = 3, lwd = 2)
abline(h = p_and$c_maximo[1], col = "black", lty = 2)
legend("topleft", c("T(t) observado", "Crítico puntual (0.05)", "Crítico máximo"),
       col = c("red", "blue", "black"), lty = c(1, 3, 2), lwd = c(2, 2, 1), bty = "n", cex = 0.75)
dev.off()

# Figura 2: Prueba t funcional de dos muestras por permutación en cada piso altitudinal
nombres_pisos_t <- c("Pacífico bajo (<700 m)", "Valle (700-1239 m)",
                     "Montaña media (1240-2000 m)", "Alta montaña (>2000 m)")
png("Figuras/06_comparacion_t_funcional_pisos.png", width = 11, height = 8, units = "in", res = 150)
par(mfrow = c(2, 2))
for (nom_p in nombres_pisos_t) {
  pp <- resultados_t[[nom_p]]$por_pentada
  plot(pentadas_eda, pp$T_welch, type = "l", lwd = 2, col = "red", xlab = "Péntada", ylab = "T(t) de Welch",
       main = paste0(nom_p, " (P=", resultados_t[[nom_p]]$resumen$n_estaciones, ")"),
       ylim = c(0, max(pp$T_welch, pp$c_maximo) * 1.35))
  lines(pentadas_eda, pp$c_puntual, col = "blue", lty = 3, lwd = 2)
  abline(h = pp$c_maximo[1], col = "black", lty = 2)
  legend("topright", c("T(t) observado", "Crítico puntual (0.05)", "Crítico máximo"),
         col = c("red", "blue", "black"), lty = c(1, 3, 2), lwd = c(2, 2, 1), bty = "n", cex = 0.72)
}
dev.off()
# En la red completa (P=147), la DE media de estaciones (17,98 mm) casi duplica la del satélite (9,44 mm) por el cambio de soporte
# En promedio anual las dos fuentes casi se igualan (+0,11 mm) porque la subestimación en Pacífico bajo (+33,91 mm, p = 0,0148)
# y Alta montaña (+2,24 mm, p = 0,0104) se compensa con la sobreestimación en Valle (-3,85 mm, p = 0,0002 en 107 de 143 puntos)
# En el dominio andino >=700 m (P=140), T(t) supera c_maximo = 2,54 en cuatro tramos (53 puntos; p_global = 0,0002)

# Diagnóstico de autocorrelación espacial entre píxeles vecinos (Índice de Moran)
lat0_pixeles <- mean(metadatos_pixeles$latitud)
x_pixeles_km <- metadatos_pixeles$longitud * 111.32 * cos(lat0_pixeles * pi / 180)
y_pixeles_km <- metadatos_pixeles$latitud * 110.57
dist_pixeles_km <- as.matrix(dist(cbind(x_pixeles_km, y_pixeles_km)))
umbral_moran_km <- unname(quantile(dist(cbind(x_pixeles_km, y_pixeles_km)), 0.05))
W_moran <- ifelse(dist_pixeles_km > 0 & dist_pixeles_km < umbral_moran_km, 1, 0)
suma_filas_W <- rowSums(W_moran)
W_moran[suma_filas_W > 0, ] <- W_moran[suma_filas_W > 0, ] / suma_filas_W[suma_filas_W > 0]

calcular_indice_moran <- function(z, W, nombre_variable) {
  n <- length(z)
  z_centrado <- z - mean(z)
  denom <- sum(z_centrado^2)
  I_obs <- (n / sum(W)) * sum(W * (z_centrado %o% z_centrado)) / denom
  set.seed(2026)
  I_nulos <- replicate(B_comparacion, {
    zp <- sample(z_centrado)
    (n / sum(W)) * sum(W * (zp %o% zp)) / denom
  })
  data.frame(
    variable = nombre_variable,
    pixeles = n,
    umbral_vecindad_km = umbral_moran_km,
    I_Moran = I_obs,
    esperanza_H0 = -1 / (n - 1),
    p_permutacion = p_maximo(I_nulos, I_obs)
  )
}

media_est_pixel <- colMeans(X1_pixeles)
media_sat_pixel <- colMeans(X2_pixeles)
dif_media_pixel <- media_est_pixel - media_sat_pixel
residuo_piso_pixel <- dif_media_pixel - ave(dif_media_pixel, metadatos_pixeles$piso_altitudinal)
moran_pixeles <- rbind(
  calcular_indice_moran(media_est_pixel, W_moran, "Media por píxel - Estaciones"),
  calcular_indice_moran(media_sat_pixel, W_moran, "Media por píxel - Satélite"),
  calcular_indice_moran(dif_media_pixel, W_moran, "Diferencia media por píxel"),
  calcular_indice_moran(residuo_piso_pixel, W_moran, "Residuo centrado por piso altitudinal")
)
moran_pixeles
guardar_tabla(moran_pixeles, "06_comparacion_t_moran_espacial.csv")
# Con umbral del 5 % de distancias más cortas (16,01 km), existe autocorrelación espacial positiva tanto en
# estaciones (I = 0,4218) y satélite (I = 0,5795) como en los residuos dentro de cada piso altitudinal (I = 0,4454; p = 0,0002),
# lo que indica que el tamaño muestral espacial efectivo es menor que P = 147

##### ANOVA funcional por pisos altitudinales #####
# Comparación entre la zonificación cartográfica IGAC/IDEAM (que pierde 66 estaciones y deja Costa con S = 1)
# y la zonificación por pisos altitudinales (Ospina Galindez, 2016; Hurtado-Montoya y Mesa-Sánchez, 2014) sin pérdida de estaciones
archivo_cartografia <- "data/Paisajes_IGAC_hidrografia_IDEAM.gpkg"
paisajes <- vect(archivo_cartografia, layer = "paisajes_IGAC")
hidrografia <- vect(archivo_cartografia, layer = "hidrografia_IDEAM")
estaciones_zonas <- unique(metadatos_eda[, c("STATION_ID", "pixel", "longitud", "latitud", "altitud_completa", "altitud_DEM")])
stopifnot(!anyDuplicated(estaciones_zonas$STATION_ID))
puntos_zonas <- vect(estaciones_zonas, geom = c("longitud", "latitud"), crs = "EPSG:4326")
paisaje_punto <- terra::extract(paisajes, project(puntos_zonas, crs(paisajes)))
cuenca_punto <- terra::extract(hidrografia, project(puntos_zonas, crs(hidrografia)))
stopifnot(identical(as.integer(paisaje_punto$id.y), seq_len(nrow(estaciones_zonas))),
          identical(as.integer(cuenca_punto$id.y), seq_len(nrow(estaciones_zonas))))
estaciones_zonas$en_departamento <- is.related(project(puntos_zonas, crs(valle)), valle, "intersects")
estaciones_zonas$paisaje <- paisaje_punto$PAISAJE
estaciones_zonas$relieve <- paisaje_punto$TIPO_RELIEVE
estaciones_zonas$area_hidrografica <- cuenca_punto$nom_ah
estaciones_zonas$zona_hidrografica <- cuenca_punto$nombre_zh

zona_relieve <- function(paisaje, area_hidrografica, zona_hidrografica) {
  case_when(
    grepl("^Montaña", paisaje) ~ "Montaña",
    grepl("^Planicie", paisaje) & area_hidrografica == "Pacifico" ~ "Costa",
    grepl("^Planicie", paisaje) & zona_hidrografica == "Cauca" ~ "Valle",
    TRUE ~ NA_character_
  )
}
estaciones_zonas$zona_IGAC <- with(estaciones_zonas, zona_relieve(paisaje, area_hidrografica, zona_hidrografica))
estaciones_zonas$zona_IGAC[!estaciones_zonas$en_departamento] <- NA_character_

# Altitud de referencia (con respaldo en altitud_DEM para las 4 estaciones >= 2000 m excluidas en altitud_completa)
niveles_pisos_4 <- c("Pacífico bajo (<700 m)", "Valle (700-1239 m)",
                     "Montaña media (1240-2000 m)", "Alta montaña (>2000 m)")
niveles_pisos_3 <- niveles_pisos_4[-1]
estaciones_zonas$altitud_m <- ifelse(!is.na(estaciones_zonas$altitud_completa),
                                     estaciones_zonas$altitud_completa,
                                     estaciones_zonas$altitud_DEM)
estaciones_zonas$piso_altitudinal <- case_when(
  estaciones_zonas$altitud_m < 700 ~ "Pacífico bajo (<700 m)",
  estaciones_zonas$altitud_m < 1240 ~ "Valle (700-1239 m)",
  estaciones_zonas$altitud_m <= 2000 ~ "Montaña media (1240-2000 m)",
  TRUE ~ "Alta montaña (>2000 m)"
)
stopifnot(!anyNA(estaciones_zonas$piso_altitudinal))
table(estaciones_zonas$zona_IGAC, estaciones_zonas$piso_altitudinal, useNA = "ifany")
guardar_tabla(estaciones_zonas, "06_comparacion_estaciones_subregion.csv")

# Mapa comparativo: zonificación cartográfica IGAC/IDEAM frente a pisos altitudinales (176 estaciones, 147 píxeles)
montanas <- paisajes[grepl("^Montaña", paisajes$PAISAJE), ]
planicies <- paisajes[grepl("^Planicie", paisajes$PAISAJE), "PAISAJE"]
planicies <- terra::intersect(makeValid(planicies), makeValid(hidrografia[, c("nom_ah", "nombre_zh")]))
planicies$zona <- zona_relieve(planicies$PAISAJE, planicies$nom_ah, planicies$nombre_zh)
colores_zonas_igac <- c(Costa = "steelblue", Montaña = "sienna", Valle = "darkgreen")
colores_puntos_igac <- unname(colores_zonas_igac[estaciones_zonas$zona_IGAC])
colores_puntos_igac[is.na(colores_puntos_igac)] <- "grey65"
colores_pisos <- c("Pacífico bajo (<700 m)" = "steelblue",
                   "Valle (700-1239 m)" = "darkgreen",
                   "Montaña media (1240-2000 m)" = "darkorange2",
                   "Alta montaña (>2000 m)" = "firebrick")
colores_puntos_pisos <- unname(colores_pisos[estaciones_zonas$piso_altitudinal])

png("Figuras/06_estaciones_subregion.png", width = 12, height = 6, units = "in", res = 150)
par(mfrow = c(1, 2))
plot(valle, col = "grey97", border = "grey35", main = "Zonificación IGAC/IDEAM (66 sin asignar)")
plot(montanas, col = "wheat2", border = NA, add = TRUE)
plot(planicies[planicies$zona == "Costa", ], col = "lightblue", border = NA, add = TRUE)
plot(planicies[planicies$zona == "Valle", ], col = "darkseagreen2", border = NA, add = TRUE)
lines(valle, col = "grey35")
points(project(puntos_zonas, crs(valle)), pch = 20, col = colores_puntos_igac, cex = 0.9)
legend("topleft", legend = c(names(colores_zonas_igac), "Sin asignar (66)"),
       col = c(colores_zonas_igac, "grey65"), pch = 20, bty = "n", cex = 0.78)
plot(valle, col = "grey97", border = "grey35", main = "Pisos altitudinales (176 estaciones, 147 píxeles)")
plot(montanas, col = "grey92", border = NA, add = TRUE)
lines(valle, col = "grey35")
points(project(puntos_zonas, crs(valle)), pch = 20, col = colores_puntos_pisos, cex = 0.9)
legend("topleft", legend = niveles_pisos_4, col = unname(colores_pisos[niveles_pisos_4]),
       pch = 20, bty = "n", cex = 0.75)
dev.off()

# Curvas climatológicas de diferencia por píxel único D_p(t) = X1_p(t) - X2_p(t) y cobertura anual por piso altitudinal
D_pixeles <- X1_pixeles - X2_pixeles
colnames(D_pixeles) <- as.character(pixeles_unicos)
piso_par <- metadatos_pixeles$piso_altitudinal[match(metadatos_comparacion$pixel, metadatos_pixeles$pixel)]
cobertura_pisos <- table(anio = metadatos_comparacion$YEAR,
                         piso = factor(piso_par, levels = niveles_pisos_4))
cobertura_anova <- data.frame(anio = as.integer(rownames(cobertura_pisos)),
                              as.data.frame.matrix(cobertura_pisos),
                              todos_activos = rowSums(cobertura_pisos > 0) == length(niveles_pisos_4),
                              check.names = FALSE)
cobertura_anova
guardar_tabla(cobertura_anova, "06_comparacion_cobertura_ANOVA.csv")

# Revisión de medias y dispersión espacial entre píxeles dentro de cada piso altitudinal
medias_pisos <- sapply(niveles_pisos_4, function(p) {
  rowMeans(D_pixeles[, metadatos_pixeles$piso_altitudinal == p, drop = FALSE])
})
de_pisos <- sapply(niveles_pisos_4, function(p) {
  apply(D_pixeles[, metadatos_pixeles$piso_altitudinal == p, drop = FALSE], 1, sd)
})
media_global_pisos <- rowMeans(D_pixeles)
dispersion_zonas <- data.frame(
  piso_altitudinal = niveles_pisos_4,
  estaciones = sapply(niveles_pisos_4, function(p) sum(estaciones_zonas$piso_altitudinal == p)),
  pixeles = sapply(niveles_pisos_4, function(p) sum(metadatos_pixeles$piso_altitudinal == p)),
  pares_pixel_anio = sapply(niveles_pisos_4, function(p) sum(piso_par == p)),
  sesgo_medio_mm = colMeans(medias_pisos),
  sesgo_anual_mm = 72 * colMeans(medias_pisos),
  DE_min_mm = apply(de_pisos, 2, min),
  DE_max_mm = apply(de_pisos, 2, max),
  pentada_max_DE = pentadas_eda[apply(de_pisos, 2, which.max)],
  DE_media_mm = colMeans(de_pisos),
  row.names = NULL
)
dispersion_zonas
guardar_tabla(dispersion_zonas, "06_comparacion_dispersion_zonas.csv")

png("Figuras/06_comparacion_ANOVA_medias_DE.png", width = 11, height = 4.5, units = "in", res = 150)
par(mfrow = c(1, 2))
matplot(pentadas_eda, medias_pisos, type = "l", lty = 1, lwd = 2, col = colores_pisos[niveles_pisos_4],
        ylim = range(medias_pisos, media_global_pisos) + c(-2, 12),
        xlab = "Péntada", ylab = "Estación - satélite (mm)", main = "Medias funcionales por piso altitudinal")
lines(pentadas_eda, media_global_pisos, lwd = 2, lty = 2, col = "black")
legend("topleft", c(niveles_pisos_4, "Global (P=147)"),
       col = c(unname(colores_pisos[niveles_pisos_4]), "black"),
       lty = c(rep(1, length(niveles_pisos_4)), 2), lwd = 2, bty = "n", cex = 0.68)
matplot(pentadas_eda, de_pisos, type = "l", lty = 1, lwd = 2, col = colores_pisos[niveles_pisos_4],
        ylim = c(0, max(de_pisos) * 1.25),
        xlab = "Péntada", ylab = "Desviación estándar entre píxeles (mm)",
        main = "DE funcional por piso altitudinal")
legend("topleft", niveles_pisos_4, col = unname(colores_pisos[niveles_pisos_4]),
       lwd = 2, bty = "n", cex = 0.68)
dev.off()
# La zonificación por pisos altitudinales conserva las 176 estaciones (147 píxeles) y separa un gradiente monótono:
# Valle (-3,85 mm/péntada; DE media 5,51 mm), Montaña media (-0,46 mm; DE 7,91 mm), Alta montaña (+2,24 mm; DE 11,18 mm)
# y Pacífico bajo (+33,91 mm; DE 35,15 mm con P = 7 píxeles independientes en vez de S = 1 estación)

# Función ANOVA funcional (SSR_t, SSE_t, SST_t, F_t, Zhang L2 y F_n,
# T_GPF de Zhang y Liang, permutación de etiquetas con c_maximo y prueba heterocedástica V_n de Cuevas et al., 2004)
w_eda <- rep(1, length(pentadas_eda))
w_eda[c(1, length(pentadas_eda))] <- 0.5

ejecutar_anova_funcional <- function(X, Grupo, nombre_esquema, w = w_eda, malla = pentadas_eda) {
  Grupo <- factor(Grupo, levels = unique(levels(Grupo)[levels(Grupo) %in% Grupo]))
  niveles <- levels(Grupo)
  k <- length(niveles)
  n_i <- table(Grupo)
  n <- sum(n_i)
  M <- nrow(X)

  Z_grupo <- model.matrix(~ Grupo - 1)
  Medias <- X %*% Z_grupo %*% diag(1 / as.numeric(n_i))
  colnames(Medias) <- niveles
  Media_global <- rowMeans(X)
  Varianzas <- sapply(niveles, function(g) apply(X[, Grupo == g, drop = FALSE], 1, var))
  DE_grupos <- sqrt(Varianzas)

  # Sumas de cuadrados puntuales y comprobación de la identidad SST(t) = SSR(t) + SSE(t)
  SSR_t <- as.numeric((Medias - Media_global)^2 %*% as.numeric(n_i))
  SSE_t <- rowSums((X - Medias[, as.integer(Grupo)])^2)
  SST_t <- rowSums((X - Media_global)^2)
  stopifnot(max(abs(SST_t - (SSR_t + SSE_t))) < 1e-9)
  F_t <- (SSR_t / (k - 1)) / (SSE_t / (n - k))

  # Pruebas globales L2 y tipo F con corrección por sesgo finito de trazas de Zhang (2013)
  SSR_n <- sum(w * SSR_t)
  SSE_n <- sum(w * SSE_t)
  SST_n <- sum(w * SST_t)
  N_gl <- n - k
  Resid <- X - Medias[, as.integer(Grupo)]
  G_comb <- Resid %*% t(Resid) / N_gl
  DG <- G_comb * w
  tr_gamma <- sum(diag(DG))
  tr_gamma2 <- sum(DG * t(DG))
  tr_gamma2_c <- N_gl^2 / ((N_gl - 1) * (N_gl + 2)) * (tr_gamma2 - tr_gamma^2 / N_gl)
  tr2_gamma_c <- N_gl * (N_gl + 1) / ((N_gl - 1) * (N_gl + 2)) * (tr_gamma^2 - 2 * tr_gamma2 / (N_gl + 1))
  beta_hat <- tr_gamma2_c / tr_gamma
  d_hat <- tr2_gamma_c / tr_gamma2_c
  p_L2 <- pchisq(SSR_n / beta_hat, df = (k - 1) * d_hat, lower.tail = FALSE)
  F_n <- (SSR_n / (k - 1)) / (SSE_n / (n - k))
  p_F_param <- pf(F_n, (k - 1) * d_hat, (n - k) * d_hat, lower.tail = FALSE)

  # Prueba F puntual globalizada T_GPF (Zhang y Liang, 2014)
  T_GPF <- sum(w * F_t)
  s_t <- sqrt(diag(G_comb))
  Rho <- G_comb / (s_t %o% s_t)
  DR <- Rho * w
  tr_rho <- sum(diag(DR))
  tr_rho2 <- sum(DR * t(DR))
  beta_rho <- tr_rho2 / ((k - 1) * tr_rho)
  d_rho <- (k - 1) * tr_rho^2 / tr_rho2
  p_GPF_param <- pchisq(T_GPF / beta_rho, df = d_rho, lower.tail = FALSE)

  # Calibración por permutación de etiquetas de piso altitudinal para F_n, T_GPF y F_max
  set.seed(2026)
  F_t_nulos <- matrix(NA_real_, M, B_comparacion)
  SSR_n_nulos <- numeric(B_comparacion)
  n_vec <- as.numeric(n_i)
  diag_inv <- diag(1 / n_vec)
  for (b in seq_len(B_comparacion)) {
    Zb <- Z_grupo[sample(n), ]
    Med_b <- X %*% Zb %*% diag_inv
    SSR_b <- as.numeric((Med_b - Media_global)^2 %*% n_vec)
    F_t_nulos[, b] <- (SSR_b / (k - 1)) / ((SST_t - SSR_b) / (n - k))
    SSR_n_nulos[b] <- sum(w * SSR_b)
  }
  F_n_nulos <- (SSR_n_nulos / (k - 1)) / ((SST_n - SSR_n_nulos) / (n - k))
  GPF_nulos <- colSums(w * F_t_nulos)
  Fmax_nulos <- apply(F_t_nulos, 2, max)

  F_max <- max(F_t)
  pico_F <- which.max(F_t)
  p_F_perm <- p_maximo(F_n_nulos, F_n)
  p_GPF_perm <- p_maximo(GPF_nulos, T_GPF)
  p_Fmax_perm <- p_maximo(Fmax_nulos, F_max)
  c_puntual <- apply(F_t_nulos, 1, quantile, probs = 1 - alfa_comparacion)
  c_maximo <- unname(quantile(Fmax_nulos, 1 - alfa_comparacion))
  p_puntual <- (1 + rowSums(F_t_nulos >= F_t - 1e-12)) / (B_comparacion + 1)

  # Prueba heterocedástica de Cuevas, Febrero y Fraiman (2004): estadístico V_n con covarianzas por grupo
  pares_idx <- combn(k, 2)
  V_n <- sum(apply(pares_idx, 2, function(p) {
    n_i[p[1]] * sum(w * (Medias[, p[1]] - Medias[, p[2]])^2)
  }))
  Raiz_cov <- lapply(niveles, function(g) {
    R_g <- X[, Grupo == g, drop = FALSE] - Medias[, g]
    E_g <- eigen(R_g %*% t(R_g) / (n_i[g] - 1), symmetric = TRUE)
    r <- min(n_i[g] - 1, M)
    E_g$vectors[, 1:r, drop = FALSE] %*% diag(sqrt(pmax(E_g$values[1:r], 0)), r)
  })
  set.seed(2026)
  V_nulos <- numeric(B_comparacion)
  for (b in seq_len(B_comparacion)) {
    Zs <- sapply(Raiz_cov, function(A_cov) A_cov %*% rnorm(ncol(A_cov)))
    V_nulos[b] <- sum(apply(pares_idx, 2, function(p) {
      sum(w * (Zs[, p[1]] - sqrt(n_i[p[1]] / n_i[p[2]]) * Zs[, p[2]])^2)
    }))
  }
  p_CFF <- p_maximo(V_nulos, V_n)

  tabla_anova <- data.frame(
    esquema = nombre_esquema,
    Fuente = c("Entre grupos", "Dentro de grupos", "Total"),
    gl = c(k - 1, n - k, n - 1),
    SC_integrada = c(SSR_n, SSE_n, SST_n),
    CM_integrado = c(SSR_n / (k - 1), SSE_n / (n - k), NA_real_),
    F_n = c(F_n, NA_real_, NA_real_),
    gl1_adj = c((k - 1) * d_hat, NA_real_, NA_real_),
    gl2_adj = c((n - k) * d_hat, NA_real_, NA_real_),
    Valor_p_param = c(p_F_param, NA_real_, NA_real_),
    Valor_p_perm = c(p_F_perm, NA_real_, NA_real_)
  )

  resumen <- data.frame(
    esquema = nombre_esquema,
    grupos = k,
    pixeles = n,
    gl_entre = k - 1,
    gl_dentro = N_gl,
    d_Zhang = d_hat,
    SSR_n = SSR_n,
    SSE_n = SSE_n,
    p_L2 = p_L2,
    F_n = F_n,
    p_F_param = p_F_param,
    p_F_perm = p_F_perm,
    T_GPF = T_GPF,
    d_rho = d_rho,
    p_GPF_param = p_GPF_param,
    p_GPF_perm = p_GPF_perm,
    pentada_pico = malla[pico_F],
    F_max = F_max,
    c_maximo = c_maximo,
    p_Fmax_perm = p_Fmax_perm,
    puntos_sig_Fmax = sum(F_t > c_maximo),
    V_n_CFF = V_n,
    p_CFF = p_CFF
  )

  list(resumen = resumen, tabla_anova = tabla_anova, Medias = Medias, Media_global = Media_global,
       DE_grupos = DE_grupos, SSR_t = SSR_t, SSE_t = SSE_t, SST_t = SST_t, F_t = F_t,
       c_puntual = c_puntual, c_maximo = c_maximo, p_puntual = p_puntual,
       F_n = F_n, F_n_nulos = F_n_nulos, d_hat = d_hat, tr_gamma = tr_gamma,
       N_gl = N_gl, n_i = n_i)
}

# Esquema principal E5: Dominio andino bimodal de 3 pisos (>= 700 m, P = 140) siguiendo a Ospina Galindez (2016)
sel_andino <- metadatos_pixeles$altitud_m >= 700
anova_E5 <- ejecutar_anova_funcional(
  D_pixeles[, sel_andino, drop = FALSE],
  factor(metadatos_pixeles$piso_altitudinal[sel_andino], levels = niveles_pisos_3),
  "E5: Dominio andino 3 pisos (>=700 m, P=140)"
)

# Esquema completo E3b: Red completa de 4 pisos altitudinales incluyendo Pacífico bajo (<700 m, P = 147)
anova_E3b <- ejecutar_anova_funcional(
  D_pixeles,
  factor(metadatos_pixeles$piso_altitudinal, levels = niveles_pisos_4),
  "E3b: Red completa 4 pisos (P=147)"
)

anova_zonas <- rbind(anova_E5$resumen, anova_E3b$resumen)
anova_zonas
guardar_tabla(anova_zonas, "06_comparacion_ANOVA_zonas.csv")

tabla_anova_zonas <- rbind(anova_E5$tabla_anova, anova_E3b$tabla_anova)
tabla_anova_zonas
guardar_tabla(tabla_anova_zonas, "06_comparacion_ANOVA_tabla.csv")

anova_por_pentada <- rbind(
  data.frame(
    esquema = "E5: Dominio andino 3 pisos (>=700 m, P=140)",
    pentada = pentadas_eda,
    Pacífico_bajo_mm = NA_real_,
    Valle_mm = anova_E5$Medias[, "Valle (700-1239 m)"],
    Montaña_media_mm = anova_E5$Medias[, "Montaña media (1240-2000 m)"],
    Alta_montaña_mm = anova_E5$Medias[, "Alta montaña (>2000 m)"],
    SC_entre = anova_E5$SSR_t,
    SC_dentro = anova_E5$SSE_t,
    F = anova_E5$F_t,
    c_puntual = anova_E5$c_puntual,
    c_maximo = anova_E5$c_maximo,
    p_puntual = anova_E5$p_puntual,
    significativa = anova_E5$F_t > anova_E5$c_maximo
  ),
  data.frame(
    esquema = "E3b: Red completa 4 pisos (P=147)",
    pentada = pentadas_eda,
    Pacífico_bajo_mm = anova_E3b$Medias[, "Pacífico bajo (<700 m)"],
    Valle_mm = anova_E3b$Medias[, "Valle (700-1239 m)"],
    Montaña_media_mm = anova_E3b$Medias[, "Montaña media (1240-2000 m)"],
    Alta_montaña_mm = anova_E3b$Medias[, "Alta montaña (>2000 m)"],
    SC_entre = anova_E3b$SSR_t,
    SC_dentro = anova_E3b$SSE_t,
    F = anova_E3b$F_t,
    c_puntual = anova_E3b$c_puntual,
    c_maximo = anova_E3b$c_maximo,
    p_puntual = anova_E3b$p_puntual,
    significativa = anova_E3b$F_t > anova_E3b$c_maximo
  )
)
guardar_tabla(anova_por_pentada, "06_comparacion_ANOVA_por_pentada.csv")

# Tramos del año donde F(t) supera el valor crítico máximo simultáneo c_maximo en E5 y E3b
extraer_intervalos_anova <- function(obj_anova, nombre_esquema) {
  iv <- intervalos_pentadas(obj_anova$F_t > obj_anova$c_maximo)
  if (nrow(iv) == 0) return(data.frame())
  iv$esquema <- nombre_esquema
  for (k in seq_len(nrow(iv))) {
    idx <- which(pentadas_eda >= iv$inicio_pentada[k] & pentadas_eda <= iv$fin_pentada[k])
    j <- idx[which.max(obj_anova$F_t[idx])]
    iv$pico_pentada[k] <- pentadas_eda[j]
    iv$pico_F[k] <- obj_anova$F_t[j]
    iv$Pacífico_bajo_mm[k] <- if ("Pacífico bajo (<700 m)" %in% colnames(obj_anova$Medias)) {
      obj_anova$Medias[j, "Pacífico bajo (<700 m)"]
    } else NA_real_
    iv$Valle_mm[k] <- obj_anova$Medias[j, "Valle (700-1239 m)"]
    iv$Montaña_media_mm[k] <- obj_anova$Medias[j, "Montaña media (1240-2000 m)"]
    iv$Alta_montaña_mm[k] <- obj_anova$Medias[j, "Alta montaña (>2000 m)"]
  }
  iv[, c("esquema", "inicio_pentada", "fin_pentada", "puntos", "pico_pentada", "pico_F",
         "Pacífico_bajo_mm", "Valle_mm", "Montaña_media_mm", "Alta_montaña_mm")]
}
intervalos_anova <- rbind(
  extraer_intervalos_anova(anova_E5, "E5: Dominio andino 3 pisos (>=700 m, P=140)"),
  extraer_intervalos_anova(anova_E3b, "E3b: Red completa 4 pisos (P=147)")
)
intervalos_anova
guardar_tabla(intervalos_anova, "06_comparacion_ANOVA_intervalos.csv")

png("Figuras/06_comparacion_ANOVA_permutacion.png", width = 11, height = 8, units = "in", res = 150)
par(mfrow = c(2, 2))
plot(pentadas_eda, anova_E5$F_t, type = "l", lwd = 2, col = "red",
     xlab = "Péntada", ylab = "F(t)",
     main = "F(t) puntual - E5 Dominio andino 3 pisos (P=140)",
     ylim = c(0, max(anova_E5$F_t, anova_E5$c_maximo) * 1.35))
lines(pentadas_eda, anova_E5$c_puntual, col = "blue", lty = 3, lwd = 2)
abline(h = anova_E5$c_maximo, col = "black", lty = 2)
legend("topleft", c("F(t) observado", "Crítico puntual (0.05)", "Crítico máximo"),
       col = c("red", "blue", "black"), lty = c(1, 3, 2), lwd = c(2, 2, 1), bty = "n", cex = 0.75)
hist(anova_E5$F_n_nulos, breaks = 50, freq = FALSE, col = "grey85", border = "white",
     xlim = c(0, max(anova_E5$F_n_nulos, anova_E5$F_n) * 1.25),
     xlab = expression(F[n]), main = "Permutación de F integrado - E5 (3 pisos)")
abline(v = anova_E5$F_n, col = "red", lwd = 2)
legend("topright", c("F observado"), col = "red", lwd = 2, bty = "n", cex = 0.8)
plot(pentadas_eda, anova_E3b$F_t, type = "l", lwd = 2, col = "red",
     xlab = "Péntada", ylab = "F(t)",
     main = "F(t) puntual - E3b Red completa 4 pisos (P=147)",
     ylim = c(0, max(anova_E3b$F_t, anova_E3b$c_maximo) * 1.35))
lines(pentadas_eda, anova_E3b$c_puntual, col = "blue", lty = 3, lwd = 2)
abline(h = anova_E3b$c_maximo, col = "black", lty = 2)
legend("topleft", c("F(t) observado", "Crítico puntual (0.05)", "Crítico máximo"),
       col = c("red", "blue", "black"), lty = c(1, 3, 2), lwd = c(2, 2, 1), bty = "n", cex = 0.75)
hist(anova_E3b$F_n_nulos, breaks = 50, freq = FALSE, col = "grey85", border = "white",
     xlim = c(0, max(anova_E3b$F_n_nulos, anova_E3b$F_n) * 1.15),
     xlab = expression(F[n]), main = "Permutación de F integrado - E3b (4 pisos)")
abline(v = anova_E3b$F_n, col = "red", lwd = 2)
legend("topright", c("F observado"), col = "red", lwd = 2, bty = "n", cex = 0.8)
dev.off()
# En el dominio andino bimodal E5 (P=140, gl = 2 y 137), todas las pruebas rechazan H0:
# F_n = 5,95 (p_param = 0,00095; p_perm = 0,0022), T_GPF = 969,69 (p_perm = 0,0004), F_max = 14,02 en la péntada 39,5
# (c_maximo = 4,67; p_Fmax = 0,0002 en 97 de 143 puntos) y la prueba heterocedástica de Cuevas et al. da V_n = 562.485,6 (p_CFF = 0,0040)
# Al incluir Pacífico bajo en E3b (P=147, gl = 3 y 143), F_n sube a 33,48 (p_perm = 0,0002), F_max = 60,11 (143/143 puntos) y p_CFF = 0,0006

# Comparaciones por pares entre pisos altitudinales independientes (T(t) de Welch por permutación y F de dos muestras con ajuste de Holm)
comparar_pares_pisos <- function(X, Grupo, obj_anova, nombre_esquema) {
  Grupo <- factor(Grupo, levels = colnames(obj_anova$Medias))
  niveles <- levels(Grupo)
  pares_lista <- combn(niveles, 2, simplify = FALSE)
  n_pares <- length(pares_lista)
  res_pares <- data.frame()
  iv_pares <- data.frame()

  for (i in seq_len(n_pares)) {
    par_g <- pares_lista[[i]]
    A <- X[, Grupo == par_g[1], drop = FALSE]
    B_mat <- X[, Grupo == par_g[2], drop = FALSE]
    nA <- ncol(A)
    nB <- ncol(B_mat)
    n_ab <- nA + nB
    est_ab <- Estadisticos_t(A, B_mat)

    # Estadístico F de dos muestras con covarianza combinada global
    F_par_glob <- (nA * nB / (nA + nB)) * sum(w_eda * est_ab$dif^2) / obj_anova$tr_gamma
    p_F_glob <- pf(F_par_glob, obj_anova$d_hat, obj_anova$N_gl * obj_anova$d_hat, lower.tail = FALSE)

    # Estadístico T(t) de Welch de dos muestras independientes por permutación de etiquetas entre los dos pisos
    X_ab <- cbind(A, B_mat)
    set.seed(2026)
    Grupos_A <- replicate(B_comparacion, sort(sample(n_ab, nA)))
    T_nulos_ab <- matrix(NA_real_, nrow(X), B_comparacion)
    for (b in seq_len(B_comparacion)) {
      g1 <- Grupos_A[, b]
      T_nulos_ab[, b] <- Estadisticos_t(X_ab[, g1, drop = FALSE], X_ab[, -g1, drop = FALSE])$T_welch
    }
    T_max_nulos_ab <- apply(T_nulos_ab, 2, max)
    pico <- which.max(est_ab$T_welch)
    c_max_05 <- unname(quantile(T_max_nulos_ab, 1 - alfa_comparacion))
    c_max_adj <- unname(quantile(T_max_nulos_ab, 1 - alfa_comparacion / n_pares))
    p_T_welch <- p_maximo(T_max_nulos_ab, est_ab$T_welch[pico])

    res_pares <- rbind(res_pares, data.frame(
      esquema = nombre_esquema,
      piso_A = par_g[1],
      piso_B = par_g[2],
      contraste = paste(par_g[1], "-", par_g[2]),
      n_A = nA,
      n_B = nB,
      dif_media_mm = mean(est_ab$dif),
      pico_pentada = pentadas_eda[pico],
      dif_pico_mm = est_ab$dif[pico],
      EE_welch_pico_mm = est_ab$ee_welch[pico],
      T_welch_max = est_ab$T_welch[pico],
      c_max_05 = c_max_05,
      c_max_adj = c_max_adj,
      p_T_welch = p_T_welch,
      F_dos_muestras = F_par_glob,
      p_F = p_F_glob
    ))

    iv <- intervalos_pentadas(est_ab$T_welch > c_max_adj)
    if (nrow(iv) > 0) {
      iv$esquema <- nombre_esquema
      iv$contraste <- paste(par_g[1], "-", par_g[2])
      for (k in seq_len(nrow(iv))) {
        idx <- which(pentadas_eda >= iv$inicio_pentada[k] & pentadas_eda <= iv$fin_pentada[k])
        j <- idx[which.max(est_ab$T_welch[idx])]
        iv$dif_min_mm[k] <- min(est_ab$dif[idx])
        iv$dif_max_mm[k] <- max(est_ab$dif[idx])
        iv$pico_pentada[k] <- pentadas_eda[j]
        iv$pico_T_welch[k] <- est_ab$T_welch[j]
      }
      iv_pares <- rbind(iv_pares, iv[, c("esquema", "contraste", "inicio_pentada", "fin_pentada",
                                         "puntos", "dif_min_mm", "dif_max_mm",
                                         "pico_pentada", "pico_T_welch")])
    }
  }
  res_pares$p_Holm_T_welch <- p.adjust(res_pares$p_T_welch, method = "holm")
  res_pares$p_Holm_F <- p.adjust(res_pares$p_F, method = "holm")
  list(resumen = res_pares, intervalos = iv_pares)
}

pares_E5 <- comparar_pares_pisos(
  D_pixeles[, sel_andino, drop = FALSE],
  factor(metadatos_pixeles$piso_altitudinal[sel_andino], levels = niveles_pisos_3),
  anova_E5, "E5: Dominio andino 3 pisos (>=700 m, P=140)"
)
pares_E3b <- comparar_pares_pisos(
  D_pixeles,
  factor(metadatos_pixeles$piso_altitudinal, levels = niveles_pisos_4),
  anova_E3b, "E3b: Red completa 4 pisos (P=147)"
)

anova_pares <- rbind(pares_E5$resumen, pares_E3b$resumen)
anova_pares
guardar_tabla(anova_pares, "06_comparacion_ANOVA_pares.csv")

intervalos_pares <- rbind(pares_E5$intervalos, pares_E3b$intervalos)
intervalos_pares
guardar_tabla(intervalos_pares, "06_comparacion_ANOVA_pares_intervalos.csv")

# Mapa de calor de valores p ajustados por Holm (T(t) de Welch por permutación) en los 4 pisos altitudinales
matriz_p_pisos <- matrix(NA_real_, length(niveles_pisos_4), length(niveles_pisos_4),
                         dimnames = list(niveles_pisos_4, niveles_pisos_4))
for (r in seq_len(nrow(pares_E3b$resumen))) {
  pA <- pares_E3b$resumen$piso_A[r]
  pB <- pares_E3b$resumen$piso_B[r]
  matriz_p_pisos[pA, pB] <- pares_E3b$resumen$p_Holm_T_welch[r]
  matriz_p_pisos[pB, pA] <- pares_E3b$resumen$p_Holm_T_welch[r]
}
datos_p_pisos <- as.data.frame(as.table(matriz_p_pisos)) %>%
  filter(!is.na(Freq)) %>%
  mutate(Decision = ifelse(Freq < alfa_comparacion, "p < 0.05", "p >= 0.05"))
png("Figuras/06_comparacion_ANOVA_pares_Holm.png", width = 7.5, height = 6, units = "in", res = 150)
print(
  ggplot(datos_p_pisos, aes(Var1, Var2, fill = Decision)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.4f", Freq)), size = 3.5) +
    scale_fill_manual(values = c("p < 0.05" = "cyan", "p >= 0.05" = "grey85"), name = "Holm") +
    labs(x = "", y = "", title = "Valores p ajustados por Holm (T(t) de Welch entre pisos)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
)
dev.off()
# Con T(t) de Welch de dos muestras independientes y ajuste de Holm, todos los pares de pisos difieren al 5 %:
# En E5, Valle frente a Montaña media (p_Holm = 0,0054; pico en péntada 10,5 con -4,56 mm), Valle frente a Alta montaña
# (p_Holm = 0,0056; pico en péntada 3,5 con -7,69 mm) y Montaña media frente a Alta montaña (p_Holm = 0,0292; pico en péntada 42 con -6,44 mm)
# En E3b, los tres contrastes con Pacífico bajo alcanzan p_Holm <= 0,0152 (T_welch entre 4,25 y 4,89 en las péntadas 56,5 a 58)


##### ANOVA funcional por fase ENSO #####
# Fase ENSO mensual (criterio NOAA en enso_evaluacion) asignada a las 72 péntadas (6 por mes)
fases_enso <- c("El Niño", "Neutral", "La Niña")
colores_enso <- c("El Niño" = "firebrick", Neutral = "grey40", "La Niña" = "steelblue")
idx_72 <- match(tiempos, pentadas_eda)
D_red_enso <- D_red[idx_72, , drop = FALSE]
colnames(D_red_enso) <- as.character(anios_red)
mes_pentada <- ceiling(tiempos / 6)

fase_mes_enso <- matrix(NA_character_, nrow = 12, ncol = n_red,
                        dimnames = list(1:12, as.character(anios_red)))
for (r in seq_len(nrow(enso_evaluacion))) {
  anio_col <- as.character(enso_evaluacion$YEAR[r])
  if (anio_col %in% colnames(fase_mes_enso)) {
    fase_mes_enso[enso_evaluacion$mes_central[r], anio_col] <- enso_evaluacion$fase_enso[r]
  }
}
fase_pentada_enso <- fase_mes_enso[mes_pentada, , drop = FALSE]
stopifnot(!anyNA(fase_pentada_enso))

# Cobertura de años por fase en cada mes
cobertura_enso_mes <- data.frame(
  mes = 1:12,
  sapply(fases_enso, function(f) rowSums(fase_mes_enso == f)),
  check.names = FALSE
)
stopifnot(all(as.matrix(cobertura_enso_mes[, fases_enso]) >= 2))
cobertura_enso_mes
guardar_tabla(cobertura_enso_mes, "06_comparacion_ENSO_cobertura.csv")
# Cada mes reúne entre 4 y 7 años en El Niño y en La Niña, y entre 9 y 15 en Neutral
# La asignación mensual conserva las transiciones dentro del año en 18 de los 23 años

# Revisión de medias y dispersión por fase ENSO en las 72 péntadas
medias_enso <- sapply(fases_enso, function(f) {
  rowSums(D_red_enso * (fase_pentada_enso == f)) / rowSums(fase_pentada_enso == f)
})
media_global_enso <- rowMeans(D_red_enso)
de_enso <- sapply(fases_enso, function(f) {
  sapply(seq_along(tiempos), function(i) sd(D_red_enso[i, fase_pentada_enso[i, ] == f]))
})
dispersion_enso <- data.frame(
  fase = fases_enso,
  anios_mes_min = sapply(fases_enso, function(f) min(cobertura_enso_mes[[f]])),
  anios_mes_max = sapply(fases_enso, function(f) max(cobertura_enso_mes[[f]])),
  media_anual_mm = colMeans(medias_enso),
  DE_min_mm = apply(de_enso, 2, min),
  DE_max_mm = apply(de_enso, 2, max),
  pentada_max_DE = tiempos[apply(de_enso, 2, which.max)],
  DE_media_mm = colMeans(de_enso)
)
dispersion_enso
guardar_tabla(dispersion_enso, "06_comparacion_ENSO_dispersion.csv")

png("Figuras/06_comparacion_ENSO_medias_DE.png", width = 11, height = 4.5, units = "in", res = 150)
par(mfrow = c(1, 2))
matplot(tiempos, medias_enso, type = "l", lty = 1, lwd = 2, col = colores_enso,
        ylim = range(medias_enso, media_global_enso) + c(0, 6),
        xlab = "Péntada", ylab = "Estación - satélite (mm)", main = "Medias funcionales por fase ENSO")
lines(tiempos, media_global_enso, lwd = 2, lty = 2)
legend("top", c(fases_enso, "Global"), col = c(colores_enso, "black"), lty = c(rep(1, length(fases_enso)), 2),
       lwd = 2, bty = "n", ncol = 4, cex = 0.75)
matplot(tiempos, de_enso, type = "l", lty = 1, lwd = 2, col = colores_enso,
        ylim = c(0, max(de_enso) * 1.25),
        xlab = "Péntada", ylab = "Desviación estándar (mm)", main = "Desviación estándar funcional por fase ENSO")
legend("top", fases_enso, col = colores_enso, lwd = 2, bty = "n", ncol = 3, cex = 0.75)
dev.off()
# El sesgo medio pasa de -3,64 mm en El Niño a -0,29 mm en Neutral y +2,88 mm en La Niña
# Con suavizado simétrico las DE medias quedan entre 3,81 y 4,99 mm y desaparece el pico de borde en La Niña (DE máxima 8,78 mm)

# ANOVA funcional por fase ENSO con permutación de años completos
anova_fase_enso <- function(Y, Fmat) {
  med_glob <- rowMeans(Y)
  SST <- rowSums((Y - med_glob)^2)
  SSR <- numeric(nrow(Y))
  for (f in fases_enso) {
    ind <- Fmat == f
    nf <- rowSums(ind)
    mf <- rowSums(Y * ind) / nf
    SSR <- SSR + nf * (mf - med_glob)^2
  }
  SSE <- SST - SSR
  F_t <- (SSR / (length(fases_enso) - 1)) / (SSE / (ncol(Y) - length(fases_enso)))
  list(SSR = SSR, SSE = SSE, SST = SST, F_t = F_t)
}

obs_enso <- anova_fase_enso(D_red_enso, fase_pentada_enso)
gl_enso <- length(fases_enso) - 1
gl_error_enso <- n_red - length(fases_enso)
F_aov_enso <- sapply(seq_along(tiempos), function(i) {
  summary(aov(D_red_enso[i, ] ~ factor(fase_pentada_enso[i, ])))[[1]][1, "F value"]
})
stopifnot(max(abs(obs_enso$SST - (obs_enso$SSR + obs_enso$SSE))) < 1e-10,
          max(abs(obs_enso$F_t - F_aov_enso)) < 1e-10)

w_72 <- rep(1, length(tiempos))
w_72[c(1, length(tiempos))] <- 0.5
F_obs_enso <- obs_enso$F_t
F_n_enso <- (sum(w_72 * obs_enso$SSR) / gl_enso) / (sum(w_72 * obs_enso$SSE) / gl_error_enso)
GPF_enso <- sum(w_72 * F_obs_enso)

# Permutar años completos conserva la autocorrelación anual de cada curva y la cobertura mensual de fases
set.seed(2026)
perms_enso <- replicate(B_comparacion, sample(n_red))
F_t_nulos_enso <- matrix(NA_real_, length(tiempos), B_comparacion)
F_n_nulos_enso <- numeric(B_comparacion)
GPF_nulos_enso <- numeric(B_comparacion)
for (b in seq_len(B_comparacion)) {
  rb <- anova_fase_enso(D_red_enso[, perms_enso[, b]], fase_pentada_enso)
  F_t_nulos_enso[, b] <- rb$F_t
  F_n_nulos_enso[b] <- (sum(w_72 * rb$SSR) / gl_enso) / (sum(w_72 * rb$SSE) / gl_error_enso)
  GPF_nulos_enso[b] <- sum(w_72 * rb$F_t)
}
F_max_nulos_enso <- apply(F_t_nulos_enso, 2, max)
c_puntual_enso <- apply(F_t_nulos_enso, 1, quantile, probs = 1 - alfa_comparacion)
c_maximo_enso <- unname(quantile(F_max_nulos_enso, 1 - alfa_comparacion))
p_puntual_enso <- (1 + rowSums(F_t_nulos_enso >= F_obs_enso - 1e-12)) / (B_comparacion + 1)

pico_F_enso <- which.max(F_obs_enso)
anova_enso <- data.frame(
  anios = n_red, gl_fases = gl_enso, gl_error = gl_error_enso,
  pentada_max = tiempos[pico_F_enso],
  F_max = F_obs_enso[pico_F_enso], c_maximo = c_maximo_enso,
  p_F_max = p_maximo(F_max_nulos_enso, F_obs_enso[pico_F_enso]),
  F_integrado = F_n_enso, p_F_integrado = p_maximo(F_n_nulos_enso, F_n_enso),
  GPF = GPF_enso, p_GPF = p_maximo(GPF_nulos_enso, GPF_enso)
)
anova_enso
guardar_tabla(anova_enso, "06_comparacion_ENSO_ANOVA.csv")

anova_enso_por_pentada <- data.frame(
  pentada = tiempos, mes = mes_pentada, medias_enso,
  SC_fases = obs_enso$SSR, SC_error = obs_enso$SSE,
  F = F_obs_enso, c_puntual = c_puntual_enso, c_maximo = c_maximo_enso,
  p_puntual = p_puntual_enso, significativa = F_obs_enso > c_maximo_enso,
  check.names = FALSE
)
guardar_tabla(anova_enso_por_pentada, "06_comparacion_ENSO_por_pentada.csv")

intervalos_enso <- intervalos_pentadas(F_obs_enso > c_maximo_enso, malla = tiempos)
for (k in seq_len(nrow(intervalos_enso))) {
  idx <- which(tiempos >= intervalos_enso$inicio_pentada[k] & tiempos <= intervalos_enso$fin_pentada[k])
  j <- idx[which.max(F_obs_enso[idx])]
  intervalos_enso$pico_pentada[k] <- tiempos[j]
  intervalos_enso$pico_F[k] <- F_obs_enso[j]
  intervalos_enso$El_Nino_mm[k] <- medias_enso[j, "El Niño"]
  intervalos_enso$Neutral_mm[k] <- medias_enso[j, "Neutral"]
  intervalos_enso$La_Nina_mm[k] <- medias_enso[j, "La Niña"]
}
intervalos_enso
guardar_tabla(intervalos_enso, "06_comparacion_ENSO_intervalos.csv")

png("Figuras/06_comparacion_ENSO_permutacion.png", width = 11, height = 4.5, units = "in", res = 150)
par(mfrow = c(1, 2))
plot(tiempos, F_obs_enso, type = "l", lwd = 2, col = "red",
     xlab = "Péntada", ylab = "F(t)",
     main = "Prueba F puntual por permutación (ENSO)",
     ylim = c(0, max(F_obs_enso, c_maximo_enso) * 1.45))
lines(tiempos, c_puntual_enso, col = "blue", lty = 3, lwd = 2)
abline(h = c_maximo_enso, col = "black", lty = 2)
legend("topleft", c("Estadístico observado", "Valor crítico puntual (0.05)", "Valor crítico máximo"),
       col = c("red", "blue", "black"), lty = c(1, 3, 2), lwd = c(2, 2, 1), bty = "n", cex = 0.75)
hist(F_n_nulos_enso, breaks = 50, freq = FALSE, col = "grey85", border = "white",
     xlim = c(0, max(F_n_nulos_enso, F_n_enso) * 1.35),
     xlab = expression(F[n]), main = "Distribución por permutación de F (ENSO)")
abline(v = F_n_enso, col = "red", lwd = 2)
legend("topright", c("F observado"), col = "red", lwd = 2, bty = "n", cex = 0.8)
dev.off()
# Las pruebas globales integradas rechazan H0: F_n = 3,36 (p = 0,0010) y GPF = 253,55 (p = 0,0018), mientras F_max = 7,67 en la péntada 38 da p = 0,0584 (c_maximo = 7,93)
# F_n y GPF muestran señal más fuerte que F_max porque el corrimiento entre fases se sostiene durante casi todo el año sin el pico artificial de borde en la péntada 68

# Comparaciones por pares entre fases ENSO tras rechazar H0, con ajuste de Holm
pares_enso <- combn(fases_enso, 2, simplify = FALSE)
n_pares_enso <- length(pares_enso)
anova_pares_enso <- data.frame()
intervalos_pares_enso <- data.frame()
curvas_pares_enso <- vector("list", n_pares_enso)

par_fase_enso <- function(Y, Fmat, f1, f2, SSE_glob) {
  ind1 <- Fmat == f1
  ind2 <- Fmat == f2
  n1 <- rowSums(ind1)
  n2 <- rowSums(ind2)
  m1 <- rowSums(Y * ind1) / n1
  m2 <- rowSums(Y * ind2) / n2
  v1 <- rowSums(((Y - m1) * ind1)^2) / (n1 - 1)
  v2 <- rowSums(((Y - m2) * ind2)^2) / (n2 - 1)
  sp2 <- ((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2)
  t_comb <- (m1 - m2) / sqrt(sp2 * (1 / n1 + 1 / n2))
  T_welch <- abs(m1 - m2) / sqrt(v1 / n1 + v2 / n2)
  SSR_par <- (n1 * n2 / (n1 + n2)) * (m1 - m2)^2
  F_par <- sum(w_72 * SSR_par) / (sum(w_72 * SSE_glob) / gl_error_enso)
  list(dif = m1 - m2, t_comb = t_comb, T_welch = T_welch, F_par = F_par)
}

SSE_nulos_enso <- matrix(NA_real_, length(tiempos), B_comparacion)
for (b in seq_len(B_comparacion)) {
  SSE_nulos_enso[, b] <- anova_fase_enso(D_red_enso[, perms_enso[, b]], fase_pentada_enso)$SSE
}

for (i in seq_len(n_pares_enso)) {
  pf <- pares_enso[[i]]
  obs_par <- par_fase_enso(D_red_enso, fase_pentada_enso, pf[1], pf[2], obs_enso$SSE)
  F_par_nulos <- numeric(B_comparacion)
  t_max_par_nulos <- numeric(B_comparacion)
  T_welch_max_nulos <- numeric(B_comparacion)
  for (b in seq_len(B_comparacion)) {
    p_b <- par_fase_enso(D_red_enso[, perms_enso[, b]], fase_pentada_enso, pf[1], pf[2], SSE_nulos_enso[, b])
    F_par_nulos[b] <- p_b$F_par
    t_max_par_nulos[b] <- max(abs(p_b$t_comb))
    T_welch_max_nulos[b] <- max(p_b$T_welch)
  }
  pico_par <- which.max(abs(obs_par$t_comb))
  pico_welch <- which.max(obs_par$T_welch)
  c_adj_par <- unname(quantile(t_max_par_nulos, 1 - alfa_comparacion / n_pares_enso))
  c_adj_welch <- unname(quantile(T_welch_max_nulos, 1 - alfa_comparacion / n_pares_enso))
  anova_pares_enso <- rbind(anova_pares_enso, data.frame(
    fase_A = pf[1], fase_B = pf[2],
    contraste = paste(pf[1], "-", pf[2]),
    dif_media_mm = mean(obs_par$dif),
    pico_pentada = tiempos[pico_par],
    dif_pico_mm = obs_par$dif[pico_par],
    t_max = abs(obs_par$t_comb[pico_par]),
    c_max_adj = c_adj_par,
    p_t_max = p_maximo(t_max_par_nulos, abs(obs_par$t_comb[pico_par])),
    T_welch_max = obs_par$T_welch[pico_welch],
    c_max_adj_welch = c_adj_welch,
    p_T_welch = p_maximo(T_welch_max_nulos, obs_par$T_welch[pico_welch]),
    F_integrado = obs_par$F_par,
    p_F_integrado = p_maximo(F_par_nulos, obs_par$F_par)
  ))
  curvas_pares_enso[[i]] <- list(contraste = paste(pf[1], "-", pf[2]),
                                 dif = obs_par$dif, t_comb = obs_par$t_comb, c_adj = c_adj_par)
}
anova_pares_enso$p_Holm_F <- p.adjust(anova_pares_enso$p_F_integrado, method = "holm")
anova_pares_enso$p_Holm_t <- p.adjust(anova_pares_enso$p_t_max, method = "holm")
anova_pares_enso$p_Holm_T_welch <- p.adjust(anova_pares_enso$p_T_welch, method = "holm")
anova_pares_enso <- anova_pares_enso[order(anova_pares_enso$p_F_integrado), ]
anova_pares_enso
guardar_tabla(anova_pares_enso, "06_comparacion_ENSO_pares.csv")

for (i in seq_len(n_pares_enso)) {
  cp <- curvas_pares_enso[[i]]
  iv <- intervalos_pentadas(abs(cp$t_comb) > cp$c_adj, malla = tiempos)
  if (nrow(iv) == 0) next
  iv$contraste <- cp$contraste
  for (k in seq_len(nrow(iv))) {
    idx <- which(tiempos >= iv$inicio_pentada[k] & tiempos <= iv$fin_pentada[k])
    j <- idx[which.max(abs(cp$t_comb[idx]))]
    iv$dif_min_mm[k] <- min(cp$dif[idx])
    iv$dif_max_mm[k] <- max(cp$dif[idx])
    iv$pico_pentada[k] <- tiempos[j]
    iv$pico_t[k] <- cp$t_comb[j]
  }
  intervalos_pares_enso <- rbind(intervalos_pares_enso, iv[, c("contraste", "inicio_pentada", "fin_pentada",
                                                               "puntos", "dif_min_mm", "dif_max_mm",
                                                               "pico_pentada", "pico_t")])
}
intervalos_pares_enso
guardar_tabla(intervalos_pares_enso, "06_comparacion_ENSO_pares_intervalos.csv")

matriz_p_enso <- matrix(NA, length(fases_enso), length(fases_enso), dimnames = list(fases_enso, fases_enso))
for (r in seq_len(nrow(anova_pares_enso))) {
  matriz_p_enso[anova_pares_enso$fase_A[r], anova_pares_enso$fase_B[r]] <- anova_pares_enso$p_Holm_F[r]
  matriz_p_enso[anova_pares_enso$fase_B[r], anova_pares_enso$fase_A[r]] <- anova_pares_enso$p_Holm_F[r]
}
datos_p_enso <- as.data.frame(as.table(matriz_p_enso)) %>%
  filter(!is.na(Freq)) %>%
  mutate(Decision = ifelse(Freq < alfa_comparacion, "p < 0.05", "p >= 0.05"))
png("Figuras/06_comparacion_ENSO_pares_Holm.png", width = 6.5, height = 5.5, units = "in", res = 150)
print(
  ggplot(datos_p_enso, aes(Var1, Var2, fill = Decision)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.4f", Freq)), size = 3.5) +
    scale_fill_manual(values = c("p < 0.05" = "cyan", "p >= 0.05" = "grey85"), name = "Holm") +
    labs(x = "", y = "", title = "Valores p ajustados por Holm (prueba F por permutación)") +
    theme_minimal()
)
dev.off()
# Tras ajustar por Holm en F integrado, difieren El Niño - La Niña (p = 0,0006; p_T_welch sin ajuste = 0,0216) y El Niño - Neutral (p = 0,0352)
# Neutral - La Niña no alcanza significancia al 5 % (p = 0,0914)

# Comprobación por fuente separada: respuesta de estaciones y satélite a las fases ENSO
E_red_enso <- sapply(anios_red, function(anio) {
  rowMeans(valores_comparacion$terrestre[idx_72, metadatos_comparacion$YEAR == anio, drop = FALSE])
})
S_red_enso <- sapply(anios_red, function(anio) {
  rowMeans(valores_comparacion$satelital[idx_72, metadatos_comparacion$YEAR == anio, drop = FALSE])
})
fuentes_enso_lista <- list(Estaciones = E_red_enso, Satélite = S_red_enso)
medias_fuentes_enso <- list()
anova_fuentes_enso <- data.frame()

for (nom_f in names(fuentes_enso_lista)) {
  Y_f <- fuentes_enso_lista[[nom_f]]
  obs_f <- anova_fase_enso(Y_f, fase_pentada_enso)
  F_n_f <- (sum(w_72 * obs_f$SSR) / gl_enso) / (sum(w_72 * obs_f$SSE) / gl_error_enso)
  F_max_f <- max(obs_f$F_t)
  F_n_nul_f <- numeric(B_comparacion)
  F_max_nul_f <- numeric(B_comparacion)
  for (b in seq_len(B_comparacion)) {
    rb <- anova_fase_enso(Y_f[, perms_enso[, b]], fase_pentada_enso)
    F_n_nul_f[b] <- (sum(w_72 * rb$SSR) / gl_enso) / (sum(w_72 * rb$SSE) / gl_error_enso)
    F_max_nul_f[b] <- max(rb$F_t)
  }
  meds_f <- sapply(fases_enso, function(f) {
    rowSums(Y_f * (fase_pentada_enso == f)) / rowSums(fase_pentada_enso == f)
  })
  medias_fuentes_enso[[nom_f]] <- meds_f
  anova_fuentes_enso <- rbind(anova_fuentes_enso, data.frame(
    fuente = nom_f,
    media_El_Nino_mm = mean(meds_f[, "El Niño"]),
    media_Neutral_mm = mean(meds_f[, "Neutral"]),
    media_La_Nina_mm = mean(meds_f[, "La Niña"]),
    F_max = F_max_f, p_F_max = p_maximo(F_max_nul_f, F_max_f),
    F_integrado = F_n_f, p_F_integrado = p_maximo(F_n_nul_f, F_n_f)
  ))
}
anova_fuentes_enso
guardar_tabla(anova_fuentes_enso, "06_comparacion_ENSO_fuentes.csv")

limites_fuentes_enso <- range(unlist(medias_fuentes_enso)) + c(0, 8)
png("Figuras/06_comparacion_ENSO_fuentes.png", width = 11, height = 4.5, units = "in", res = 150)
par(mfrow = c(1, 2))
for (nom_f in names(fuentes_enso_lista)) {
  matplot(tiempos, medias_fuentes_enso[[nom_f]], type = "l", lty = 1, lwd = 2, col = colores_enso,
          ylim = limites_fuentes_enso, xlab = "Péntada", ylab = "Precipitación media (mm)",
          main = paste("Medias por fase ENSO -", nom_f))
  lines(tiempos, rowMeans(fuentes_enso_lista[[nom_f]]), lwd = 2, lty = 2)
  legend("top", c(fases_enso, "Global"), col = c(colores_enso, "black"),
         lty = c(rep(1, length(fases_enso)), 2), lwd = 2, bty = "n", ncol = 4, cex = 0.75)
}
dev.off()
# En estaciones la lluvia media varía con ENSO (21,62 mm en El Niño; 25,22 en Neutral; 27,57 en La Niña; F_n = 3,39, p = 0,0006)
# En satélite las tres curvas permanecen alrededor de 25 mm (25,26; 25,51; 24,69 mm; F_n = 1,12, p = 0,3429), por lo que sobrestima en El Niño y subestima en La Niña

# Sensibilidad por cohorte de píxeles (dominio andino >= 700 m y cohorte fija de 23 años) y diagnóstico de autocorrelación serial lag-1
pixeles_23_anios <- metadatos_pixeles$pixel[metadatos_pixeles$anios_activos == 23]
pixeles_andinos_enso <- metadatos_pixeles$pixel[metadatos_pixeles$altitud_m >= 700]
cohortes_enso_lista <- list(
  "Red completa (147 píxeles)" = diferencia_por_anio(rep(TRUE, nrow(metadatos_comparacion)), anios_red)[idx_72, , drop = FALSE],
  "Dominio andino >=700 m (140 píxeles)" = diferencia_por_anio(metadatos_comparacion$pixel %in% pixeles_andinos_enso, anios_red)[idx_72, , drop = FALSE],
  "Cohorte fija 23 años (59 píxeles)" = diferencia_por_anio(metadatos_comparacion$pixel %in% pixeles_23_anios, anios_red)[idx_72, , drop = FALSE]
)
sensibilidad_enso_cohortes <- do.call(rbind, lapply(names(cohortes_enso_lista), function(nom_c) {
  Y_c <- cohortes_enso_lista[[nom_c]]
  obs_c <- anova_fase_enso(Y_c, fase_pentada_enso)
  F_n_c <- (sum(w_72 * obs_c$SSR) / gl_enso) / (sum(w_72 * obs_c$SSE) / gl_error_enso)
  GPF_c <- sum(w_72 * obs_c$F_t)
  F_max_c <- max(obs_c$F_t)
  F_n_nul_c <- numeric(B_comparacion)
  GPF_nul_c <- numeric(B_comparacion)
  F_max_nul_c <- numeric(B_comparacion)
  for (b in seq_len(B_comparacion)) {
    rb <- anova_fase_enso(Y_c[, perms_enso[, b]], fase_pentada_enso)
    F_n_nul_c[b] <- (sum(w_72 * rb$SSR) / gl_enso) / (sum(w_72 * rb$SSE) / gl_error_enso)
    GPF_nul_c[b] <- sum(w_72 * rb$F_t)
    F_max_nul_c[b] <- max(rb$F_t)
  }
  meds_c <- sapply(fases_enso, function(f) {
    rowSums(Y_c * (fase_pentada_enso == f)) / rowSums(fase_pentada_enso == f)
  })
  serie_anual_c <- colMeans(Y_c)
  acf1_anual <- acf(serie_anual_c, lag.max = 1, plot = FALSE)$acf[2]
  acf1_pentadas <- apply(Y_c, 1, function(fila) acf(fila, lag.max = 1, plot = FALSE)$acf[2])
  data.frame(
    cohorte = nom_c,
    media_El_Nino_mm = mean(meds_c[, "El Niño"]),
    media_Neutral_mm = mean(meds_c[, "Neutral"]),
    media_La_Nina_mm = mean(meds_c[, "La Niña"]),
    F_max = F_max_c,
    p_F_max = p_maximo(F_max_nul_c, F_max_c),
    F_integrado = F_n_c,
    p_F_integrado = p_maximo(F_n_nul_c, F_n_c),
    GPF = GPF_c,
    p_GPF = p_maximo(GPF_nul_c, GPF_c),
    acf_lag1_media_anual = acf1_anual,
    acf_lag1_mediana_pentadas = median(acf1_pentadas)
  )
}))
sensibilidad_enso_cohortes
guardar_tabla(sensibilidad_enso_cohortes, "06_comparacion_ENSO_cohortes_autocorrelacion.csv")
# El corrimiento del sesgo entre El Niño y La Niña se conserva al excluir el Pacífico (<700 m) y en la cohorte balanceada de 59 píxeles
# con los 23 años completos (F_n entre 2,98 y 3,36; p <= 0,0030), y la autocorrelación interanual lag-1 de la media anual es baja (0,061 a 0,116)
# COMPARACION FUNCIONAL ANUAL TABLAS Y FIGURAS IMPORTANTES (1981-2003)
anios<-sort(unique(metadatos_comparacion$YEAR));tt<-pentadas_eda;pisos<-niveles_pisos_4
cols_pisos<-c("Pacífico bajo (<700 m)"="steelblue4","Valle (700-1239 m)"="darkgreen","Montaña media (1240-2000 m)"="darkorange2","Alta montaña (>2000 m)"="firebrick");cols_enso<-c("El Niño"="firebrick","Neutral"="grey40","La Niña"="steelblue")
guardar<-function(x,n)write.csv(x,paste0(n,".csv"),row.names=FALSE)
abrir<-function(n,mf=c(5,5),mar=c(3,3.2,2.5,.7),oma=c(2.8,2.2,1.3,.5)){par(mfrow=mf,mar=mar,oma=oma,cex.axis=1.15,cex.main=1.35,cex.lab=1.2)}
cerrar<-function(xlab="Péntada"){mtext(xlab,side=1,outer=TRUE,line=.6,cex=1.35)}
panel<-function(i){if(i>length(anios)){plot.new();return(FALSE)};TRUE}
datos<-lapply(anios,function(a){i<-which(metadatos_comparacion$YEAR==a);list(E=valores_comparacion$terrestre[,i,drop=FALSE],S=valores_comparacion$satelital[,i,drop=FALSE],D=diferencia_comparacion[,i,drop=FALSE],meta=metadatos_comparacion[i,,drop=FALSE])});names(datos)<-as.character(anios)

# EMPAREJAMIENTO ESTACION-PIXEL
emp<-estaciones_zonas[,c("STATION_ID","pixel","longitud","latitud","altitud_m","piso_altitudinal")];emp$estaciones_en_pixel<-as.integer(table(emp$pixel)[as.character(emp$pixel)])
guardar(emp,"01_emparejamiento_estacion_pixel");guardar(emp[emp$estaciones_en_pixel>1,],"01_pixeles_compartidos")
res_emp<-data.frame(estaciones=nrow(emp),pixeles_unicos=length(unique(emp$pixel)),pixeles_compartidos=length(unique(emp$pixel[emp$estaciones_en_pixel>1])),estaciones_en_pixeles_compartidos=sum(emp$estaciones_en_pixel>1));guardar(res_emp,"01_resumen_emparejamiento")
cat("Emparejamiento:",nrow(emp),"estaciones,",length(unique(emp$pixel)),"pixeles\n")

# DIFERENCIA FUNCIONAL ANUAL Y BANDA +/-2EE
bandas<-do.call(rbind,lapply(seq_along(anios),function(i){D<-datos[[i]]$D;n<-ncol(D);m<-rowMeans(D);ee<-apply(D,1,sd)/sqrt(n);data.frame(anio=anios[i],pentada=tt,pixeles=n,d=m,EE=ee,LI=m-2*ee,LS=m+2*ee,excluye_cero=(m-2*ee>0)|(m+2*ee<0))}))
guardar(bandas,"02_diferencia_bandas_anuales");lim<-range(bandas[,c("LI","LS")],0,na.rm=TRUE)
abrir("02_diferencias_23_anios")
for(i in 1:25){if(!panel(i))next;b<-bandas[bandas$anio==anios[i],];plot(tt,b$d,type="n",ylim=lim,main=paste0(anios[i]," (P=",b$pixeles[1],")"),xlab="",ylab="");polygon(c(tt,rev(tt)),c(b$LI,rev(b$LS)),col=adjustcolor("steelblue",alpha.f=.2),border=NA);lines(tt,b$d,col="navy",lwd=1.8);lines(tt,b$LI,col="steelblue",lty=2);lines(tt,b$LS,col="steelblue",lty=2);abline(h=0,col="red",lty=3);if(i==1)legend("topright",c("d(t)=estación-satélite","Banda ±2EE","Cero"),col=c("navy","steelblue","red"),lty=c(1,2,3),lwd=c(1.8,1,1),bty="n",cex=.9)}
cerrar()
intervalos<-do.call(rbind,lapply(split(bandas,bandas$anio),function(b){r<-rle(b$excluye_cero);f<-cumsum(r$lengths);ini<-f-r$lengths+1;if(!any(r$values))return(NULL);data.frame(anio=b$anio[1],inicio=b$pentada[ini[r$values]],fin=b$pentada[f[r$values]],signo=ifelse(b$d[ini[r$values]]>0,"Positivo: estación>satélite","Negativo: satélite>estación"))}))
guardar(intervalos,"02_intervalos_fuera_cero")

# CORRELACION CRUZADA FUNCIONAL Y SU DIAGONAL
corES<-lapply(datos,function(z)cor(t(z$E),t(z$S)));diagES<-sapply(corES,diag);Cmed<-apply(simplify2array(corES),c(1,2),median,na.rm=TRUE);Dmed<-apply(diagES,1,median,na.rm=TRUE)
tabla_cor<-data.frame(anio=anios,correlacion_media=round(colMeans(diagES),3),correlacion_min=round(apply(diagES,2,min),3),correlacion_max=round(apply(diagES,2,max),3));guardar(tabla_cor,"03_correlacion_anual_resumen")
pal<-colorRampPalette(c("blue3","white","red3"))(100);par(mfrow=c(1,2),mar=c(5,5,3.7,1.4),oma=c(0,0,1,0),cex.axis=1.3,cex.lab=1.4,cex.main=1.45)
image(tt,tt,Cmed,zlim=c(-1,1),col=pal,xlab="Péntada satélite",ylab="Péntada estación",main="Correlación cruzada\nMediana de 23 años",useRaster=TRUE);contour(tt,tt,Cmed,levels=seq(-1,1,.2),add=TRUE,drawlabels=TRUE,col="black",lwd=.8,labcex=1);abline(0,1,lty=3,col="black");legend("bottomleft",legend=c("Azul: negativa","Blanco: cercana a 0","Rojo: positiva"),fill=c("blue3","white","red3"),bty="n",cex=1)
matplot(tt,diagES,type="l",lty=1,col=adjustcolor("steelblue",alpha.f=.25),ylim=c(-1,1),xlab="Péntada",ylab="Correlación",main="Diagonales anuales (1981-2003)");lines(tt,Dmed,col="navy",lwd=3);abline(h=0,lty=3);legend("bottomleft",c("23 años","Mediana anual"),col=c("steelblue","navy"),lwd=c(1,3),bty="n",cex=1.1)

# TABLA DE COVARIANZA POR FUENTE (SIN GRAFICAS)
covE<-lapply(datos,function(z)cov(t(z$E)));covS<-lapply(datos,function(z)cov(t(z$S)))
tabla_cov<-data.frame(anio=anios,pixeles=sapply(datos,function(z)ncol(z$D)),var_est=round(sapply(covE,function(M)mean(diag(M))),2),var_sat=round(sapply(covS,function(M)mean(diag(M))),2),cov_est_consec=round(sapply(covE,function(M)mean(M[cbind(1:(nrow(M)-1),2:nrow(M))])),2),cov_sat_consec=round(sapply(covS,function(M)mean(M[cbind(1:(nrow(M)-1),2:nrow(M))])),2))
tabla_cov$razon_var_est_sat<-round(tabla_cov$var_est/tabla_cov$var_sat,2);guardar(tabla_cov,"04_tabla_covarianza_anual")
cat("\nTABLA DE COVARIANZA ANUAL\n");print(tabla_cov,row.names=FALSE);cat("\nAños con varianza de estaciones mayor que satélite:",sum(tabla_cov$var_est>tabla_cov$var_sat),"/",nrow(tabla_cov),"\n")

# COMPARACION POR ALTITUD DENTRO DE CADA ANO
alt<-do.call(rbind,lapply(seq_along(anios),function(i){z<-datos[[i]];g<-factor(metadatos_pixeles$piso_altitudinal[match(z$meta$pixel,metadatos_pixeles$pixel)],levels=pisos);do.call(rbind,lapply(pisos,function(p){j<-which(g==p);data.frame(anio=anios[i],piso=p,pixeles=length(j),pentada=tt,d=if(length(j))rowMeans(z$D[,j,drop=FALSE]) else NA_real_)}))}))
guardar(alt,"05_sesgo_pisos_por_anio");res_alt<-aggregate(d~anio+piso,alt,mean);res_alt$pixeles<-mapply(function(a,p)unique(alt$pixeles[alt$anio==a & alt$piso==p]),res_alt$anio,res_alt$piso);guardar(res_alt,"05_resumen_anual_altitud")
abrir("05_altitud_23_anios")
for(i in 1:25){if(i==24){plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));legend(x=.02,y=.88,legend=pisos,col=cols_pisos[pisos],lty=1,lwd=2,bty="n",cex=1.05,xjust=0,yjust=1,xpd=FALSE,seg.len=1.5);next};if(!panel(i))next;z<-alt[alt$anio==anios[i],];mat<-sapply(pisos,function(p)z$d[z$piso==p]);matplot(tt,mat,type="l",lty=1,col=cols_pisos[pisos],lwd=1.8,main=anios[i],xlab="",ylab="",ylim=range(alt$d,0,na.rm=TRUE));abline(h=0,lty=3)}
cerrar()

# ENSO: CURVA ANUAL COLOREADA POR FASE MENSUAL
enso_tab<-data.frame(anio=rep(anios,each=72),pentada=rep(tiempos,length(anios)),fase=as.vector(fase_pentada_enso),d=as.vector(D_red_enso));guardar(enso_tab,"06_ENSO_anual_pentadas");lim_enso<-range(D_red_enso,0,finite=TRUE)
abrir("06_ENSO_23_anios")
for(i in 1:25){if(i==24){plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));legend(x=.05,y=.85,legend=names(cols_enso),col=cols_enso,lty=1,lwd=2,bty="n",cex=1.35,xjust=0,yjust=1,xpd=FALSE,seg.len=2);next};if(!panel(i))next;z<-enso_tab[enso_tab$anio==anios[i],];plot(z$pentada,z$d,type="n",ylim=lim_enso,main=anios[i],xlab="",ylab="");segments(head(z$pentada,-1),head(z$d,-1),tail(z$pentada,-1),tail(z$d,-1),col=cols_enso[head(z$fase,-1)],lwd=2);points(z$pentada,z$d,col=cols_enso[z$fase],pch=16,cex=.4);abline(h=0,lty=3)}
cerrar()

#################### Regresion funcional ####################
particion_completa <- particion_curva(metadatos_completos$STATION_ID, metadatos_completos$YEAR)
idx_entrenamiento <- which(particion_completa == "entrenamiento")
idx_prueba <- which(particion_completa == "prueba")
grupo_entrenamiento <- grupo_vc(metadatos_completos$STATION_ID[idx_entrenamiento])
stopifnot(!anyNA(grupo_entrenamiento))
table(grupo_entrenamiento)

# Raiz de la curva en la misma base
malla_regresion <- seq(0, 72, by = 0.1)
raiz_suavizada <- function(coeficientes, base) {
  valores <- sqrt(pmax(eval.fd(malla_regresion, fd(t(coeficientes), base)), 0))
  smooth.basis(malla_regresion, valores, base)
}
raiz_estaciones <- raiz_suavizada(construccion_completa$terrestre$coeficientes, base_estaciones)
satelite_fd <- raiz_suavizada(construccion_completa$satelital$coeficientes, base_satelite)$fd
estaciones_fd <- raiz_estaciones$fd
a_mm <- function(y) pmax(y, 0)^2

# Covariables estandarizadas con entrenamiento
oni_curvas <- oni_anual$oni_promedio[match(metadatos_completos$YEAR, oni_anual$YEAR)]
estandarizar <- function(x) (x - mean(x[idx_entrenamiento])) / sd(x[idx_entrenamiento])
altitud_std <- estandarizar(metadatos_completos$altitud_completa)
oni_std <- estandarizar(oni_curvas)

# Modelo 1: funcion a funcion 
beta0Par <- fdPar(base_estaciones, int2Lfd(2), 0)
Ajuste_linmod <- function(lam, idx = idx_entrenamiento) {
  beta1bifd <- bifd(matrix(0, K_satelite, K_estaciones), base_satelite, base_estaciones)
  beta1Par <- bifdPar(beta1bifd, int2Lfd(2), int2Lfd(2), lam, lam)
  linmod(satelite_fd[idx], estaciones_fd[idx], list(beta0Par, beta1Par))
}
# Prediccion en las 72 pentadas
puntos_medios <- seq(0.5, 71.5, 1)
predecir_linmod <- function(ajuste, idx) {
  b0 <- eval.fd(tiempos, ajuste$beta0estfd)
  B1 <- eval.bifd(puntos_medios, tiempos, ajuste$beta1estbifd)
  X <- eval.fd(puntos_medios, satelite_fd[idx])
  as.vector(b0) + t(B1) %*% X
}

# Modelo 2: concurrente con covariables
Ajuste_concurrente <- function(lam, idx = idx_entrenamiento, y2cMap = NULL, SigmaE = NULL) {
  xfdlist <- list(const = rep(1, length(idx)), satelite = satelite_fd[idx],
                  altitud = altitud_std[idx], oni = oni_std[idx])
  betaPar <- fdPar(base_estaciones, int2Lfd(2), lam)
  betalist <- list(const = betaPar, satelite = betaPar, altitud = betaPar, oni = betaPar)
  fRegress(estaciones_fd[idx], xfdlist, betalist, y2cMap = y2cMap, SigmaE = SigmaE)
}
predecir_concurrente <- function(ajuste, idx) {
  b <- sapply(ajuste$betaestlist, function(beta) eval.fd(tiempos, beta$fd))
  b[, 1] + b[, 2] * eval.fd(tiempos, satelite_fd[idx]) +
    outer(b[, 3], altitud_std[idx]) + outer(b[, 4], oni_std[idx])
}

# Lambda elegido con validacion cruzada de 5 grupos de estaciones
# Cada grupo se predice con el modelo ajustado en los otros cuatro
prediccion_vc <- function(ajustar, predecir, lam) {
  prediccion <- matrix(NA, length(tiempos), length(idx_entrenamiento))
  for (g in 1:5) {
    fuera <- grupo_entrenamiento == g
    prediccion[, fuera] <- predecir(ajustar(lam, idx_entrenamiento[!fuera]), idx_entrenamiento[fuera])
  }
  prediccion
}
Y_entrenamiento <- eval.fd(tiempos, estaciones_fd[idx_entrenamiento])
loglam_regresion <- seq(-2, 8, 1)
ECM_linmod <- c(); ECM_concurrente <- c()
for (l in seq_along(loglam_regresion)) {
  ECM_linmod[l] <- mean((Y_entrenamiento - prediccion_vc(Ajuste_linmod, predecir_linmod, 10^loglam_regresion[l]))^2)
  ECM_concurrente[l] <- mean((Y_entrenamiento - prediccion_vc(Ajuste_concurrente, predecir_concurrente,
                                                               10^loglam_regresion[l]))^2)
}
seleccion_regresion <- data.frame(log10_lambda = loglam_regresion, ECM_VC_funcion_a_funcion = ECM_linmod,
                                  ECM_VC_concurrente = ECM_concurrente)
seleccion_regresion
guardar_tabla(seleccion_regresion, "06_regresion_seleccion_lambda.csv")
lam_linmod <- 10^loglam_regresion[which.min(ECM_linmod)]
lam_concurrente <- 10^loglam_regresion[which.min(ECM_concurrente)]

plot(loglam_regresion, ECM_linmod, type = "b", pch = 19, ylim = range(ECM_linmod, ECM_concurrente),
     xlab = expression(log[10](lambda)), ylab = "ECM de validacion cruzada (escala raiz)", main = "Eleccion de lambda")
lines(loglam_regresion, ECM_concurrente, type = "b", pch = 17, col = "red")
legend("topright", c("Funcion a funcion", "Concurrente con covariables"), col = c("black", "red"),
       pch = c(19, 17), bty = "n")
guardar_figura("06_regresion_seleccion_lambda.png", 7, 5)

# Modelos finales con todo el entrenamiento y el lambda elegido
ajuste_linmod <- Ajuste_linmod(lam_linmod)
ajuste_concurrente <- Ajuste_concurrente(lam_concurrente)

# Errores estandar del concurrente
Y_entrenamiento_malla <- eval.fd(malla_regresion, estaciones_fd[idx_entrenamiento])
residuos_malla <- Y_entrenamiento_malla - eval.fd(malla_regresion, ajuste_concurrente$yhatfdobj)
SigmaE_concurrente <- var(t(residuos_malla))
error_estandar_concurrente <- Ajuste_concurrente(lam_concurrente, y2cMap = raiz_estaciones$y2cMap,
                                                 SigmaE = SigmaE_concurrente)

# R2 funcional en entrenamiento
R2_funcional <- function(Y, Yhat) 1 - rowSums((Y - Yhat)^2) / rowSums((Y - rowMeans(Y))^2)
R2_linmod <- R2_funcional(Y_entrenamiento, predecir_linmod(ajuste_linmod, idx_entrenamiento))
R2_concurrente <- R2_funcional(Y_entrenamiento, predecir_concurrente(ajuste_concurrente, idx_entrenamiento))

# Validacion cruzada frente a la curva de la estacion (mm), con el lambda elegido
Y_entrenamiento_mm <- pmax(eval.fd(tiempos, fd(t(construccion_completa$terrestre$coeficientes[idx_entrenamiento, ]),
                                               base_estaciones)), 0)
prediccion_vc_linmod <- prediccion_vc(Ajuste_linmod, predecir_linmod, lam_linmod)
prediccion_vc_concurrente <- prediccion_vc(Ajuste_concurrente, predecir_concurrente, lam_concurrente)
predicciones_vc <- list(
  "Funcion a funcion" = a_mm(prediccion_vc_linmod),
  "Concurrente con covariables" = a_mm(prediccion_vc_concurrente),
  "Satelite sin corregir" = pmax(eval.fd(tiempos, fd(t(construccion_completa$satelital$coeficientes[idx_entrenamiento, ]),
                                                     base_satelite)), 0))
comparacion_regresion <- data.frame(
  modelo = names(predicciones_vc),
  log10_lambda = c(log10(lam_linmod), log10(lam_concurrente), NA),
  RMSE_VC_mm = sapply(predicciones_vc, function(p) sqrt(mean((p - Y_entrenamiento_mm)^2))),
  sesgo_acumulado_VC_pct = sapply(predicciones_vc, function(p) 100 * (sum(p) / sum(Y_entrenamiento_mm) - 1)),
  R2_medio_entrenamiento = c(mean(R2_linmod), mean(R2_concurrente), NA))
comparacion_regresion
guardar_tabla(comparacion_regresion, "06_regresion_comparacion.csv")

# Funcion a funcion: R2, superficie beta y correlacion de los residuos
par(mfrow = c(1, 3))
plot(tiempos, R2_linmod, type = "l", lwd = 2, xlab = "Pentada", ylab = expression(R^2), ylim = c(0, 1),
     main = expression(R^2 * "(t) en entrenamiento"))
B1_regresion <- eval.bifd(puntos_medios, tiempos, ajuste_linmod$beta1estbifd)
contour(puntos_medios, tiempos, B1_regresion, xlab = "s: pentada del satelite", ylab = "t: pentada de la estacion",
        main = expression(hat(beta)[1](s, t)), nlevels = 12)
abline(0, 1, lty = 2, col = "red")
# Correlacion de los residuos entre pentadas (Hooker, 2010)
residuos_linmod <- Y_entrenamiento - predecir_linmod(ajuste_linmod, idx_entrenamiento)
contour(tiempos, tiempos, cor(t(residuos_linmod)), xlab = "Pentada (s)", ylab = "Pentada (t)",
        main = "Correlacion de los residuos", nlevels = 10)
guardar_figura("06_regresion_funcion_a_funcion.png", 15, 5)

# Concurrente: coeficientes con bandas y R2
nombres_beta <- c("Intercepto: beta_0(t)", "Satelite: beta_1(t)", "Altitud: beta_2(t)", "ONI: beta_3(t)")
par(mfrow = c(2, 3))
for (j in 1:4) {
  bj <- eval.fd(malla_regresion, ajuste_concurrente$betaestlist[[j]]$fd)
  sej <- eval.fd(malla_regresion, error_estandar_concurrente$betastderrlist[[j]])
  matplot(malla_regresion, cbind(bj, bj + 2 * sej, bj - 2 * sej), type = "l", lty = c(1, 2, 2), lwd = c(2, 1, 1),
          col = c("black", "red", "red"), xlab = "Pentada", ylab = "", main = nombres_beta[j])
  abline(h = 0, lty = 3)
}
plot(tiempos, R2_concurrente, type = "l", lwd = 2, ylim = c(0, 1), xlab = "Pentada", ylab = expression(R^2 * (t)),
     main = "R2(t) en entrenamiento")
lines(tiempos, R2_linmod, lwd = 2, lty = 2, col = "blue")
legend("bottomright", c("Concurrente", "Funcion a funcion"), lty = c(1, 2), col = c("black", "blue"), bty = "n")
guardar_figura("06_regresion_concurrente.png", 11, 6.5)

# Error por pentada en validacion cruzada
par(mfrow = c(1, 1))
RMSE_pentada <- sapply(predicciones_vc, function(p) sqrt(rowMeans((p - Y_entrenamiento_mm)^2)))
matplot(tiempos, RMSE_pentada, type = "l", lty = 1, lwd = 2, col = c("black", "red", "darkorange"),
        xlab = "Pentada", ylab = "RMSE frente a la curva de la estacion (mm)", main = "Error por pentada en validacion cruzada")
legend("topleft", colnames(RMSE_pentada), col = c("black", "red", "darkorange"), lty = 1, lwd = 2, bty = "n")
guardar_figura("06_regresion_error_por_pentada.png", 9, 5)

# Predicciones en mm de los dos modelos
# Entrenamiento: prediccion de validacion cruzada; prueba: modelo final
curvas_regresion <- rownames(construccion_completa$terrestre$coeficientes)[c(idx_entrenamiento, idx_prueba)]
prediccion_regresion <- rbind(
  data.frame(modelo = "Funcion a funcion", curva = curvas_regresion,
             t(a_mm(cbind(prediccion_vc_linmod, predecir_linmod(ajuste_linmod, idx_prueba))))),
  data.frame(modelo = "Concurrente con covariables", curva = curvas_regresion,
             t(a_mm(cbind(prediccion_vc_concurrente, predecir_concurrente(ajuste_concurrente, idx_prueba))))))
names(prediccion_regresion)[-(1:2)] <- paste0("P", sprintf("%02d", tiempos))
write.csv(prediccion_regresion, "Curvas/regresion_funcional_predicciones.csv", row.names = FALSE)


#################### Curvas para Python ####################
# Curvas MCP, cada fuente con su base
nombres_estaciones <- paste0("c", sprintf("%02d", 1:K_estaciones))
nombres_satelite <- paste0("c", sprintf("%02d", 1:K_satelite))
nombres_pentadas <- paste0("P", sprintf("%02d", tiempos))

# Metadatos de cada curva estacion-año
metadatos_python <- metadatos_eda[, c("curva", "STATION_ID", "YEAR", "longitud", "latitud", "pixel",
                                      "altitud", "altitud_completa", "fuente_altitud")]
metadatos_python$particion <- particion_curva(metadatos_python$STATION_ID, metadatos_python$YEAR)
# Grupo de validacion cruzada de cada curva de entrenamiento (el mismo de la regresion)
metadatos_python$grupo_vc <- ifelse(metadatos_python$particion == "entrenamiento", grupo_vc(metadatos_python$STATION_ID), NA)
stopifnot(!anyNA(metadatos_python$grupo_vc[metadatos_python$particion == "entrenamiento"]))
pares_modelo <- metadatos_python$curva[metadatos_python$particion != "sin_uso"]
metadatos_python <- merge(metadatos_python, oni_anual, by = "YEAR", all.x = TRUE)
metadatos_python <- metadatos_python[match(pares_modelo, metadatos_python$curva), ]
table(metadatos_python$particion)
table(metadatos_python$grupo_vc)

# Coeficientes
coeficientes_estaciones <- data.frame(curva = pares_modelo, t(coeficientes_eda$terrestre[, pares_modelo]))
coeficientes_satelite <- data.frame(curva = pares_modelo, t(coeficientes_eda$satelital[, pares_modelo]))
names(coeficientes_estaciones)[-1] <- nombres_estaciones
names(coeficientes_satelite)[-1] <- nombres_satelite

# Curvas evaluadas en las 72 pentadas
curvas_estaciones <- data.frame(curva = pares_modelo, t(Phi_estaciones %*% coeficientes_eda$terrestre[, pares_modelo]))
curvas_satelite <- data.frame(curva = pares_modelo, t(Phi_satelite %*% coeficientes_eda$satelital[, pares_modelo]))
names(curvas_estaciones)[-1] <- nombres_pentadas
names(curvas_satelite)[-1] <- nombres_pentadas

# Datos observados por pentada
observado_estaciones <- data.frame(curva = pares_modelo, fuentes_eda$terrestre[pares_modelo, ])
observado_satelite <- data.frame(curva = pares_modelo, fuentes_eda$satelital[pares_modelo, ])
names(observado_estaciones)[-1] <- nombres_pentadas
names(observado_satelite)[-1] <- nombres_pentadas

# Bases en malla fina para la matriz de Gram
nodos_estaciones <- data.frame(nodo = c(0, base_estaciones$params, 72))
nodos_satelite <- data.frame(nodo = c(0, base_satelite$params, 72))
malla_fina <- seq(0, 72, by = 0.1)
malla_estaciones <- data.frame(t = malla_fina, eval.basis(malla_fina, base_estaciones))
malla_satelite <- data.frame(t = malla_fina, eval.basis(malla_fina, base_satelite))
names(malla_estaciones)[-1] <- nombres_estaciones
names(malla_satelite)[-1] <- nombres_satelite
write.csv(data.frame(pentada = tiempos, estaciones = pesos_estaciones, satelite = pesos_satelite),
          "Curvas/pesos_MCP.csv", row.names = FALSE)

write.csv(metadatos_python, "Curvas/metadatos_curvas.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(coeficientes_estaciones, "Curvas/coeficientes_estaciones.csv", row.names = FALSE)
write.csv(coeficientes_satelite, "Curvas/coeficientes_satelite.csv", row.names = FALSE)
write.csv(curvas_estaciones, "Curvas/curvas_estaciones_72_pentadas.csv", row.names = FALSE)
write.csv(curvas_satelite, "Curvas/curvas_satelite_72_pentadas.csv", row.names = FALSE)
write.csv(observado_estaciones, "Curvas/observado_estaciones.csv", row.names = FALSE)
write.csv(observado_satelite, "Curvas/observado_satelite.csv", row.names = FALSE)
write.csv(nodos_estaciones, "Curvas/base_estaciones_nodos.csv", row.names = FALSE)
write.csv(nodos_satelite, "Curvas/base_satelite_nodos.csv", row.names = FALSE)
write.csv(malla_estaciones, "Curvas/base_estaciones_malla_fina.csv", row.names = FALSE)
write.csv(malla_satelite, "Curvas/base_satelite_malla_fina.csv", row.names = FALSE)
write.csv(parametros_construccion, "Curvas/parametros_construccion.csv", row.names = FALSE)

# Comprobacion de las curvas exportadas
max(abs(as.matrix(curvas_estaciones[, -1]) -
          as.matrix(coeficientes_estaciones[, -1]) %*% t(Phi_estaciones)))















