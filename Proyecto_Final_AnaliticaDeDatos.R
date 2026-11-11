# PROYECTO 3
# Por: Natalia Paredes, Mariana Cifuentes, John Franco
# y Estefanía Hernández
# LIBRERIAS
#install.packages("fda")
#install.packages("fda.usc")
#install.packages("readxl")
#install.packages("tidyr")
library(fda)
library(fda.usc)
library(readxl)
library(terra)
library(dplyr)
library(ggplot2)
library(tidyterra)
library(geodata)
library(tidyr)

set.seed(2026)

# DIRECTORIO DE TRABAJO
rm(list=ls())
setwd("C:/Users/maria/Downloads/proyecto_final_analiticaDeDatos")
dir.create("Figuras", showWarnings = FALSE)

# CARGA DE DATOS: PRECIPITACION SATELITE Y TIERRA
p_tierra_1930_1979 <- read_excel("data/P_1930_1979.xlsx")
p_tierra_1980_1990 <- read.csv("data/P_1980_1990.csv")
p_tierra_1991_2003 <- read.csv("data/P_1991_2003.csv")

# de 1981 a 2013
p_satelite <- rast(list.files("data/ImagenesSatelitePentadal", 
                              pattern = "ppt.*\\.tif$", full.names = TRUE))

# Periodo principal del proyecto: 1981-2003
p_tierra <- bind_rows( 
  p_tierra_1980_1990, 
  p_tierra_1991_2003) %>%
  filter(between(YEAR, 1981, 2003))

# Conservar las columnas necesarias para construir la fecha.
# No eliminar columnas por posicion: antes se borraba YEAR/MONTH/DAY
# cuando cambiaba el orden de los archivos.
columnas_requeridas <- c("STATION_ID", "YEAR", "MONTH", "DAY", "VALUE")
faltantes <- setdiff(columnas_requeridas, names(p_tierra))
if (length(faltantes) > 0) {
  stop("Faltan columnas en p_tierra: ", paste(faltantes, collapse = ", "),
       ". Columnas disponibles: ", paste(names(p_tierra), collapse = ", "))
}
p_tierra$fecha <- as.Date(sprintf("%04d-%02d-%02d", p_tierra$YEAR, p_tierra$MONTH, p_tierra$DAY))
p_tierra <- p_tierra[order(p_tierra$STATION_ID, p_tierra$fecha), ]

# Reservar prueba sin someterla al preprocesamiento de desarrollo
p_prueba <- p_tierra %>% filter(between(YEAR, 2001, 2003))

# Aplicar los bloques existentes solo a entrenamiento y validacion
p_tierra <- p_tierra %>% filter(between(YEAR, 1981, 2000))

# Cargar coordenadas y altitudes
id_txt <- function(x) trimws(format(x, scientific = FALSE, trim = TRUE))

coord <- read.csv("data/CVC_Wheather_Stations_Long_Lat.csv",
                  fileEncoding = "latin1") %>%
  transmute(STATION_ID = id_txt(CODIGO),
            longitud = as.numeric(Lon),
            latitud = as.numeric(lat))

alt <- read_excel("data/CatalogoEstacionesElevacion.xlsx") %>%
  transmute(STATION_ID = id_txt(codigo),
            altitud = as.numeric(alt))

stopifnot(!anyDuplicated(coord$STATION_ID),
          !anyDuplicated(alt$STATION_ID))

info_estaciones <- full_join(coord, alt, by = "STATION_ID")

# Borde del Valle
departamentos_col <- gadm(country = "COL", level = 1, path = tempdir())  
valle <- departamentos_col[departamentos_col$NAME_1 == "Valle del Cauca", ]

# Valores imposibles
summary(p_tierra)

quantile(p_tierra$VALUE, probs = c(0.5, 0.9, 0.95, 0.99, 0.999, 0.99999), na.rm = TRUE)
tapply(p_tierra$VALUE, p_tierra$YEAR, summary)

# 1984: VALOR 859375
sum(p_tierra$VALUE > 300, na.rm = TRUE)

mayores_300 <- p_tierra[!is.na(p_tierra$VALUE) & p_tierra$VALUE > 300,
                        c("STATION_ID", "YEAR", "MONTH", "DAY", "VALUE")]
mayores_300 <- mayores_300[order(-mayores_300$YEAR), ]
print(head(mayores_300, 300))

# VALORES MENORES AL PERCENTIL 99.9 (120 MM)
hist(p_tierra$VALUE[p_tierra$VALUE < 120], 
     main = "Precipitación diaria",
     xlab = "Precipitación (mm)", ylab = "Frecuencia")

# OUTLIERS EN CADA ESTACION Y AÑO
casos <- unique(mayores_300[, c("STATION_ID", "YEAR")])
casos <- casos[order(casos$YEAR, casos$STATION_ID), ]
nrow(casos)   

par(mfrow = c(2, 3))

for (i in seq_len(nrow(casos))) {
  d <- p_tierra[p_tierra$STATION_ID == casos$STATION_ID[i] & p_tierra$YEAR == casos$YEAR[i], ]
  d <- d[order(d$fecha), ]
  
  plot(d$fecha, d$VALUE, type = "b", pch = 16, cex = 0.4, lwd = 0.5, col = "steelblue",
       ylim = c(0, min(max(d$VALUE, na.rm = TRUE), 1000)),
       xlab = "Fecha", ylab = "Precipitación (mm)",
       main = paste("Estación", casos$STATION_ID[i], "-", casos$YEAR[i]))
  
  fuertes <- d[!is.na(d$VALUE) & d$VALUE > 200, ]
  points(fuertes$fecha, fuertes$VALUE, col = "red", pch = 19)
}

par(mfrow = c(1, 1))

# ELIMINAR LOS VALORES AISLADOS > 400
fin_mes <- format(p_tierra$fecha + 1, "%d") == "01"

n <- nrow(p_tierra)
anterior  <- c(NA, p_tierra$VALUE[-n])
siguiente <- c(p_tierra$VALUE[-1], NA)
vecino_max <- pmax(anterior, siguiente, na.rm = TRUE)

error <- which(p_tierra$VALUE > 400 & p_tierra$VALUE > 5 * vecino_max & !fin_mes)

data.frame(p_tierra[error, c("STATION_ID", "YEAR", "MONTH", "DAY")], anterior = anterior[error],
           VALUE = p_tierra$VALUE[error], siguiente = siguiente[error])

length(error)
length(unique(p_tierra$STATION_ID[error]))

p_tierra$VALUE[error] <- NA

sum(p_tierra$VALUE > 300, na.rm = TRUE)
# ACUMULADOS MENSUALES Y CEROS FALSOS
# DATOS PARA DESARROLLO DEL PREPROCESAMIENTO (1981-2000)
p_m <- p_tierra %>%
  filter(YEAR >= 1981, YEAR <= 2000) %>%
  arrange(STATION_ID, fecha) %>%
  mutate(fin_mes = format(fecha + 1, "%d") == "01")

cat("Registros del periodo:", nrow(p_m), "\n")
cat("Estaciones:", n_distinct(p_m$STATION_ID), "\n")
cat("Pares estacion-anio:", nrow(distinct(p_m, STATION_ID, YEAR)), "\n")


# RESUMEN DE REGISTROS POR ESTACION Y MES
resumen_mensual <- p_m %>%
  group_by(STATION_ID, YEAR, MONTH) %>%
  summarise(n_registros = n_distinct(fecha),
            n_validos = sum(!is.na(VALUE)),
            n_ceros = sum(VALUE == 0, na.rm = TRUE),
            n_positivos = sum(VALUE > 0, na.rm = TRUE),
            n_na = sum(is.na(VALUE)),
            total_observado = sum(VALUE, na.rm = TRUE),
            valor_fin_mes = if (any(fin_mes & !is.na(VALUE))) {
              VALUE[which(fin_mes & !is.na(VALUE))[1]]
            } else {
              NA_real_
            },
            .groups = "drop")

# INDICADORES DE POSIBLES ACUMULADOS
resumen_mensual <- resumen_mensual %>%
  mutate(proporcion_ceros = ifelse(n_validos > 0, n_ceros / n_validos, NA_real_),
         proporcion_fin_mes = ifelse(total_observado > 0,
                                     valor_fin_mes / total_observado, NA_real_),
         # Regla exploratoria, no clasificacion definitiva
         candidato = !is.na(valor_fin_mes) &
           valor_fin_mes >= 100 &
           proporcion_fin_mes >= 0.80)


# MESES CANDIDATOS
candidatos_M <- resumen_mensual %>%
  filter(candidato) %>%
  arrange(STATION_ID, YEAR, MONTH)

cat("Meses candidatos:", nrow(candidatos_M), "\n")

print(tibble::as_tibble(candidatos_M %>%
                          select(STATION_ID, YEAR, MONTH, n_registros, n_ceros, n_positivos,
                                 valor_fin_mes, total_observado, proporcion_fin_mes)),
      n = 30)

# PATRONES QUE SE REPITEN EN VARIOS MESES DEL MISMO ANIO
patrones_M <- candidatos_M %>%
  group_by(STATION_ID, YEAR) %>%
  summarise(meses_sospechosos = n(), .groups = "drop") %>%
  arrange(desc(meses_sospechosos))

print(tibble::as_tibble(patrones_M), n = 30)

# CASOS CONOCIDOS EN DETALLE
caso_1981 <- resumen_mensual %>% filter(STATION_ID == 5421200201, YEAR == 1981)
caso_1986 <- resumen_mensual %>% filter(STATION_ID == 5300000201, YEAR == 1986)

print(tibble::as_tibble(caso_1981), n = 12)
print(tibble::as_tibble(caso_1986), n = 12)

# CEROS EN LOS MESES CANDIDATOS
# Son candidatos de revision, NO ceros falsos confirmados
ceros_M <- p_m %>%
  semi_join(candidatos_M, by = c("STATION_ID", "YEAR", "MONTH")) %>%
  filter(VALUE == 0)

cat("Ceros en meses candidatos:", nrow(ceros_M), "\n")

# CONTEOS DEL DIAGNOSTICO
diagnostico_M <- data.frame(indicador = c("Registros del periodo",
                                          "Estaciones",
                                          "Pares estacion-anio",
                                          "Meses evaluados",
                                          "Meses candidatos",
                                          "Estaciones con meses candidatos",
                                          "Pares con al menos 3 meses candidatos"),
                            cantidad = c(nrow(p_m),
                                         n_distinct(p_m$STATION_ID),
                                         nrow(distinct(p_m, STATION_ID, YEAR)),
                                         nrow(resumen_mensual),
                                         nrow(candidatos_M),
                                         n_distinct(candidatos_M$STATION_ID),
                                         sum(patrones_M$meses_sospechosos >= 3)))

print(diagnostico_M)


# TOTALES MENSUALES SEGUN LA ESTRUCTURA DE REGISTROS
# Meses con una sola observacion valida ubicada en el ultimo dia del mes
mensuales_M <- resumen_mensual %>%
  filter(n_registros == 1, n_validos == 1, !is.na(valor_fin_mes)) %>%
  arrange(STATION_ID, YEAR, MONTH)

cat("Meses con un unico registro al final:", nrow(mensuales_M), "\n")

print(tibble::as_tibble(mensuales_M %>%
                          select(STATION_ID, YEAR, MONTH, n_registros, valor_fin_mes)),
      n = 40)

# PATRONES REPETIDOS POR ESTACION Y ANIO
patrones_mensuales <- mensuales_M %>%
  group_by(STATION_ID, YEAR) %>%
  summarise(meses_con_un_registro = n(), .groups = "drop") %>%
  arrange(desc(meses_con_un_registro))

print(tibble::as_tibble(patrones_mensuales), n = 30)


# Considerar todos los meses con un unico registro al final
alta_sospecha_M <- patrones_mensuales %>%
  filter(meses_con_un_registro >= 1)

cat("Pares estacion-anio con registro unico al final:", nrow(alta_sospecha_M), "\n")

print(tibble::as_tibble(alta_sospecha_M), n = 30)

# ESTACIONES PROBLEMATICAS A TRAVES DE LOS ANIOS
problemas_estacion <- patrones_mensuales %>%
  group_by(STATION_ID) %>%
  summarise(anios_afectados = n_distinct(YEAR),
            meses_sospechosos = sum(meses_con_un_registro),
            anios = paste(sort(unique(YEAR)), collapse = ", "),
            .groups = "drop") %>%
  arrange(desc(meses_sospechosos))

print(tibble::as_tibble(problemas_estacion), n = Inf)


# AUDITORIA DE ACUMULADOS: FECHAS EXACTAS DE LOS CASOS PRIORITARIOS
auditoria_acumulados_M <- p_m %>%
  inner_join(mensuales_M %>%
               semi_join(alta_sospecha_M, by = c("STATION_ID", "YEAR")) %>%
               select(STATION_ID, YEAR, MONTH),
             by = c("STATION_ID", "YEAR", "MONTH")) %>%
  filter(fin_mes, !is.na(VALUE)) %>%
  transmute(STATION_ID, YEAR, MONTH, fecha,
            valor_original = VALUE,
            motivo = "Unico registro del mes ubicado en el ultimo dia",
            # Se excluye de la serie DIARIA solo si se confirma que es un total mensual
            decision = "Propuesta: excluir como dato diario") %>%
  arrange(STATION_ID, fecha)

cat("Registros prioritarios:", nrow(auditoria_acumulados_M), "\n")

print(tibble::as_tibble(auditoria_acumulados_M), n = Inf)


# DIAGNOSTICO DE CEROS FALSOS
# Meses con muchos ceros y un pico final
diagnostico_ceros_M <- resumen_mensual %>%
  filter(n_registros >= 25,
         n_ceros >= 20,
         n_positivos <= 3,
         !is.na(valor_fin_mes),
         valor_fin_mes >= 30,
         proporcion_fin_mes >= 0.80) %>%
  select(STATION_ID, YEAR, MONTH, n_registros, n_ceros, n_positivos,
         valor_fin_mes, proporcion_fin_mes) %>%
  arrange(STATION_ID, YEAR, MONTH)

cat("Meses con ceros potencialmente falsos:", nrow(diagnostico_ceros_M), "\n")

print(tibble::as_tibble(diagnostico_ceros_M), n = Inf)

# REGISTROS CERO PARA REVISION (NO SE MODIFICAN)
auditoria_ceros_M <- p_m %>%
  semi_join(diagnostico_ceros_M, by = c("STATION_ID", "YEAR", "MONTH")) %>%
  filter(VALUE == 0) %>%
  transmute(STATION_ID, YEAR, MONTH, fecha,
            valor_original = VALUE,
            motivo = "Cero en mes con patron sospechoso",
            decision = "Pendiente de revision")

cat("Ceros candidatos a revision:", nrow(auditoria_ceros_M), "\n")

print(tibble::as_tibble(auditoria_ceros_M), n = Inf)

# DIAS CON LLUVIA POSITIVA DENTRO DE ESOS MESES
revision_ceros <- p_m %>%
  semi_join(diagnostico_ceros_M, by = c("STATION_ID", "YEAR", "MONTH")) %>%
  filter(!is.na(VALUE), VALUE > 0) %>%
  arrange(STATION_ID, YEAR, MONTH, fecha) %>%
  group_by(STATION_ID, YEAR, MONTH) %>%
  summarise(dias_con_lluvia = n(),
            fechas_lluvia = paste(format(fecha, "%d"), collapse = ", "),
            valores_lluvia = paste(VALUE, collapse = ", "),
            .groups = "drop")

print(tibble::as_tibble(revision_ceros), n = Inf)

# CONTEXTO TEMPORAL: MES ANTERIOR, SOSPECHOSO Y POSTERIOR
meses_ceros <- diagnostico_ceros_M %>%
  distinct(STATION_ID, YEAR, MONTH) %>%
  mutate(fecha_mes = as.Date(sprintf("%04d-%02d-01", YEAR, MONTH)))

contexto_ceros_M <- bind_rows(meses_ceros %>% select(STATION_ID, fecha_mes),
                              meses_ceros %>% transmute(STATION_ID, fecha_mes = fecha_mes - 1),
                              meses_ceros %>% transmute(STATION_ID, fecha_mes = fecha_mes + 32)) %>%
  mutate(fecha_mes = as.Date(format(fecha_mes, "%Y-%m-01"))) %>%
  distinct() %>%
  mutate(YEAR = as.integer(format(fecha_mes, "%Y")),
         MONTH = as.integer(format(fecha_mes, "%m"))) %>%
  inner_join(resumen_mensual, by = c("STATION_ID", "YEAR", "MONTH")) %>%
  arrange(STATION_ID, fecha_mes) %>%
  select(STATION_ID, YEAR, MONTH, n_registros, n_ceros, n_positivos, valor_fin_mes)

print(tibble::as_tibble(contexto_ceros_M), n = Inf)

# Evaluacion historica:frecuencia de dias lluviosos por mes
historico_M <- resumen_mensual %>%
  filter(YEAR >= 1981, YEAR <= 1995, n_registros >= 25, n_validos >= 25) %>%
  mutate(frecuencia_lluvia = n_positivos / n_validos) %>%
  select(STATION_ID, YEAR, MONTH, frecuencia_lluvia)

# Comparar cada mes candidato con otros anios de la misma estacion y el mismo mes
comparacion_historica_M <- diagnostico_ceros_M %>%
  mutate(frecuencia_observada = n_positivos / n_registros) %>%
  left_join(historico_M,
            by = c("STATION_ID", "MONTH"),
            relationship = "many-to-many") %>%
  filter(is.na(YEAR.y) | YEAR.x != YEAR.y) %>%
  group_by(STATION_ID, YEAR.x, MONTH, n_ceros, n_positivos, frecuencia_observada) %>%
  summarise(n_anios_referencia = sum(!is.na(frecuencia_lluvia)),
            frecuencia_historica = if (all(is.na(frecuencia_lluvia))) {
              NA_real_
            } else {
              mean(frecuencia_lluvia, na.rm = TRUE)
            },
            .groups = "drop") %>%
  rename(YEAR = YEAR.x) %>%
  mutate(temporada = ifelse(MONTH %in% c(3, 4, 5, 9, 10, 11),
                            "Lluviosa", "Menos lluviosa"),
         diferencia = frecuencia_observada - frecuencia_historica,
         # Clasificacion para revisar, no falsedad confirmada
         prioridad = case_when(n_anios_referencia < 5 ~ "Referencia insuficiente",
                               temporada == "Lluviosa" & diferencia <= -0.25 ~ "Prioridad alta",
                               diferencia <= -0.25 ~ "Prioridad media",
                               TRUE ~ "Sin evidencia fuerte")) %>%
  arrange(prioridad, diferencia)

print(tibble::as_tibble(comparacion_historica_M), n = Inf)


# Comparacion historica de los 11 meses con ceros sospechosos
comparacion_ceros_historica_M <- resumen_mensual %>%
  semi_join(diagnostico_ceros_M, by = c("STATION_ID", "MONTH")) %>%
  filter(n_registros >= 25, n_validos >= 25, between(YEAR, 1981, 1995)) %>%
  mutate(pct_ceros = 100 * n_ceros / n_validos,
         pct_lluvia = 100 * n_positivos / n_validos) %>%
  inner_join(diagnostico_ceros_M %>%
               transmute(STATION_ID, MONTH, anio_sospechoso = YEAR,
                         ceros_sospechosos = n_ceros,
                         lluvia_fin_mes = valor_fin_mes),
             by = c("STATION_ID", "MONTH")) %>%
  filter(YEAR != anio_sospechoso) %>%
  group_by(STATION_ID, MONTH, anio_sospechoso, ceros_sospechosos, lluvia_fin_mes) %>%
  summarise(anios_historicos = n(),
            mediana_ceros = median(n_ceros),
            max_ceros_historicos = max(n_ceros),
            mediana_pct_ceros = round(median(pct_ceros), 1),
            mediana_dias_lluvia = median(n_positivos),
            .groups = "drop") %>%
  arrange(STATION_ID, anio_sospechoso)

print(tibble::as_tibble(comparacion_ceros_historica_M), n = Inf)


# Comparar los 11 meses sospechosos con otras estaciones
# y con el historico de esas estaciones (1981-1995)

meses_sospechosos_M <- diagnostico_ceros_M %>%
  transmute(STATION_ID, YEAR, MONTH,
            ceros_sospechosos = n_ceros,
            lluvia_fin_mes = valor_fin_mes,
            pct_ceros_sospechosos = round(100 * n_ceros / n_registros, 1))

referencia_estaciones_M <- resumen_mensual %>%
  filter(n_registros >= 25, n_validos >= 25) %>%
  anti_join(diagnostico_ceros_M, by = c("STATION_ID", "YEAR", "MONTH")) %>%
  mutate(pct_ceros = 100 * n_ceros / n_validos)

comparacion_otras_estaciones_M <- meses_sospechosos_M %>%
  rowwise() %>%
  mutate(
    otras_estaciones_mismo_anio = n_distinct(
      referencia_estaciones_M$STATION_ID[
        referencia_estaciones_M$YEAR == YEAR &
          referencia_estaciones_M$MONTH == MONTH &
          referencia_estaciones_M$STATION_ID != STATION_ID
      ]),
    mediana_ceros_mismo_anio = median(
      referencia_estaciones_M$pct_ceros[
        referencia_estaciones_M$YEAR == YEAR &
          referencia_estaciones_M$MONTH == MONTH &
          referencia_estaciones_M$STATION_ID != STATION_ID
      ], na.rm = TRUE),
    mediana_ceros_historica_otras = median(
      referencia_estaciones_M$pct_ceros[
        between(referencia_estaciones_M$YEAR, 1981, 1995) &
          referencia_estaciones_M$YEAR != YEAR &
          referencia_estaciones_M$MONTH == MONTH &
          referencia_estaciones_M$STATION_ID != STATION_ID
      ], na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(diferencia_mismo_anio = round(
    pct_ceros_sospechosos - mediana_ceros_mismo_anio, 1),
    diferencia_historica = round(
      pct_ceros_sospechosos - mediana_ceros_historica_otras, 1)) %>%
  arrange(desc(diferencia_mismo_anio))

print(tibble::as_tibble(comparacion_otras_estaciones_M), n = Inf, width = Inf)

# Historico de lluvia concentrada al final del mes
historico_picos_M <- diagnostico_ceros_M %>%
  transmute(STATION_ID, MONTH, anio_sospechoso = YEAR) %>%
  left_join(resumen_mensual %>%
              filter(between(YEAR, 1981, 1995),
                     n_registros >= 25, n_validos >= 25) %>%
              transmute(STATION_ID, MONTH, anio_historico = YEAR,
                        pico_final = !is.na(proporcion_fin_mes) &
                          proporcion_fin_mes >= 0.80,
                        unica_lluvia_final = n_positivos == 1 &
                          !is.na(valor_fin_mes) & valor_fin_mes > 0),
            by = c("STATION_ID", "MONTH"),
            relationship = "many-to-many") %>%
  filter(is.na(anio_historico) | anio_historico != anio_sospechoso) %>%
  group_by(STATION_ID, MONTH, anio_sospechoso) %>%
  summarise(anios_historicos = sum(!is.na(anio_historico)),
            veces_pico_final = sum(pico_final, na.rm = TRUE),
            veces_unica_lluvia_final = sum(unica_lluvia_final, na.rm = TRUE),
            .groups = "drop")

print(tibble::as_tibble(historico_picos_M), n = Inf, width = Inf)

# Detectar meses completos con precipitacion igual a cero
resumen_ceros_completos_M <- resumen_mensual %>%
  mutate(inicio_mes = as.Date(sprintf("%04d-%02d-01", YEAR, MONTH)),
         siguiente_mes = as.Date(format(inicio_mes + 32, "%Y-%m-01")),
         dias_mes = as.integer(siguiente_mes - inicio_mes),
         pct_ceros = 100 * n_ceros / n_validos)

meses_todo_cero_M <- resumen_ceros_completos_M %>%
  filter(n_registros == dias_mes,
         n_validos == dias_mes,
         n_ceros == dias_mes) %>%
  select(STATION_ID, YEAR, MONTH, dias_mes, n_ceros)

# Comparar con otros años de la misma estacion
historico_todo_cero_M <- meses_todo_cero_M %>%
  rename(anio_sospechoso = YEAR) %>%
  left_join(resumen_ceros_completos_M %>%
              filter(between(YEAR, 1981, 1995), n_validos >= 25) %>%
              transmute(STATION_ID, MONTH, anio_historico = YEAR,
                        pct_ceros_historico = pct_ceros),
            by = c("STATION_ID", "MONTH"),
            relationship = "many-to-many") %>%
  filter(is.na(anio_historico) | anio_historico != anio_sospechoso) %>%
  group_by(STATION_ID, anio_sospechoso, MONTH, dias_mes) %>%
  summarise(anios_historicos = sum(!is.na(anio_historico)),
            mediana_historica = if (all(is.na(pct_ceros_historico)))
              NA_real_ else median(pct_ceros_historico, na.rm = TRUE),
            .groups = "drop")

# Comparar con otras estaciones en el mismo mes y año
referencia_regional_M <- resumen_ceros_completos_M %>%
  filter(n_validos >= 25) %>%
  group_by(YEAR, MONTH) %>%
  summarise(mediana_regional = median(pct_ceros, na.rm = TRUE),
            estaciones_referencia = n_distinct(STATION_ID),
            .groups = "drop")

comparacion_todo_cero_M <- historico_todo_cero_M %>%
  left_join(referencia_regional_M,
            by = c("anio_sospechoso" = "YEAR", "MONTH")) %>%
  mutate(diferencia_historica = round(100 - mediana_historica, 1),
         diferencia_regional = round(100 - mediana_regional, 1)) %>%
  arrange(desc(diferencia_regional))

cat("Meses completos con todos los valores en cero:",
    nrow(meses_todo_cero_M), "\n")

print(tibble::as_tibble(head(comparacion_todo_cero_M, 30)), width = Inf)


# LIMPIEZA FINAL Y TRAZABILIDAD
p_m_depurado <- p_m

# Seleccionar exclusivamente los 33 meses candidatos a acumulados
# Criterio: lluvia final >= 100 mm y >= 80% del total mensual
acumulados_tratar_M <- p_m %>%
  semi_join(
    candidatos_M %>%
      select(STATION_ID, YEAR, MONTH),
    by = c("STATION_ID", "YEAR", "MONTH")
  ) %>%
  filter(fin_mes, !is.na(VALUE)) %>%
  transmute(
    STATION_ID, YEAR, MONTH, fecha,
    valor_original = VALUE,
    motivo = "Precipitacion final >= 100 mm y >= 80% del total mensual"
  ) %>%
  distinct(STATION_ID, fecha, .keep_all = TRUE) %>%
  arrange(STATION_ID, fecha)

# Verificar los 33 candidatos
stopifnot(
  nrow(candidatos_M) == 33,
  nrow(acumulados_tratar_M) == nrow(candidatos_M)
)

# Verificar coincidencia con la base original
verificacion_M <- acumulados_tratar_M %>%
  left_join(
    p_m %>% select(STATION_ID, fecha, VALUE),
    by = c("STATION_ID", "fecha")
  )

stopifnot(
  nrow(verificacion_M) == nrow(acumulados_tratar_M),
  all(!is.na(verificacion_M$VALUE)),
  all(verificacion_M$valor_original == verificacion_M$VALUE)
)

# Convertir en NA unicamente los 33 registros de fin de mes
p_m_depurado <- p_m_depurado %>%
  left_join(
    acumulados_tratar_M %>%
      transmute(STATION_ID, fecha, tratar_M = TRUE),
    by = c("STATION_ID", "fecha")
  ) %>%
  mutate(
    VALUE = if_else(
      !is.na(tratar_M),
      NA_real_,
      as.numeric(VALUE)
    )
  ) %>%
  select(-tratar_M)

# Registrar modificaciones de acumulados
trazabilidad_M <- acumulados_tratar_M %>%
  mutate(
    valor_nuevo = NA_real_,
    tratamiento = "Excluido de serie diaria",
    justificacion = motivo
  )

# Verificar que no queden candidatos sin tratar
acumulados_pendientes_M <- candidatos_M %>%
  anti_join(
    acumulados_tratar_M %>%
      select(STATION_ID, YEAR, MONTH),
    by = c("STATION_ID", "YEAR", "MONTH")
  )

stopifnot(nrow(acumulados_pendientes_M) == 0)

# CEROS SOSPECHOSOS
# Cinco meses con ceros sospechosos y lluvia al final
meses_ceros_tratar_M <- tibble(
  STATION_ID = c(
    2623300101, 2611400103, 2612710106,
    2622900301, 2610000105
  ),
  YEAR = c(1987, 1987, 1987, 1992, 1986),
  MONTH = c(11, 9, 7, 8, 12),
  motivo = "Ceros sospechosos y precipitacion concentrada al final"
)

# Meses completamente en cero previamente seleccionados
# Los casos de 2003 se mantienen fuera del desarrollo
meses_todo_cero_tratar_M <- tibble(
  STATION_ID = c(
    2612000113, 2612000110, 2612000112,
    2612000202, 2612700203, 2612720101,
    2622600201
  ),
  YEAR = c(1997, 2003, 2003, 1990, 1990, 1992, 1994),
  MONTH = c(1, 5, 5, 2, 9, 10, 9),
  motivo = "Mes completamente en cero con anomalia historica y regional"
)

# Unificar solo los meses de 1981-2000
meses_adicionales_M <- bind_rows(
  meses_ceros_tratar_M,
  meses_todo_cero_tratar_M
) %>%
  filter(between(YEAR, 1981, 2000))

stopifnot(
  nrow(meses_adicionales_M) == 10,
  nrow(distinct(
    meses_adicionales_M, STATION_ID, YEAR, MONTH
  )) == 10
)

# Guardar los valores originales antes de modificarlos
trazabilidad_adicional_M <- p_m_depurado %>%
  inner_join(
    meses_adicionales_M,
    by = c("STATION_ID", "YEAR", "MONTH")
  ) %>%
  filter(!is.na(VALUE)) %>%
  transmute(
    STATION_ID, YEAR, MONTH, fecha,
    valor_original = VALUE,
    valor_nuevo = NA_real_,
    tratamiento = "Excluido de serie diaria",
    justificacion = motivo
  )

# Verificar que los 10 meses tengan registros validos
stopifnot(
  nrow(distinct(
    trazabilidad_adicional_M,
    STATION_ID, YEAR, MONTH
  )) == nrow(meses_adicionales_M)
)

# Convertir todos los valores de esos 10 meses a NA
p_m_depurado <- p_m_depurado %>%
  left_join(
    meses_adicionales_M %>%
      transmute(
        STATION_ID, YEAR, MONTH,
        tratar_adicional_M = TRUE
      ),
    by = c("STATION_ID", "YEAR", "MONTH")
  ) %>%
  mutate(
    VALUE = if_else(
      !is.na(tratar_adicional_M),
      NA_real_,
      as.numeric(VALUE)
    )
  ) %>%
  select(-tratar_adicional_M)

# TRAZABILIDAD TOTAL
trazabilidad_total_M <- bind_rows(
  trazabilidad_M %>%
    select(
      STATION_ID, YEAR, MONTH, fecha,
      valor_original, valor_nuevo,
      tratamiento, justificacion
    ),
  trazabilidad_adicional_M %>%
    select(
      STATION_ID, YEAR, MONTH, fecha,
      valor_original, valor_nuevo,
      tratamiento, justificacion
    )
) %>%
  arrange(STATION_ID, fecha)

# Verificar que cada registro modificado sea unico
stopifnot(
  !anyDuplicated(
    trazabilidad_total_M[c("STATION_ID", "fecha")]
  )
)

# Verificar coincidencia entre limpieza y trazabilidad
stopifnot(
  nrow(trazabilidad_total_M) ==
    sum(!is.na(p_m$VALUE)) -
    sum(!is.na(p_m_depurado$VALUE))
)

# COMPARACION ANTES Y DESPUES
contar_pares <- function(datos) {
  datos %>%
    filter(!is.na(VALUE)) %>%
    distinct(STATION_ID, YEAR) %>%
    nrow()
}

comparacion_M <- data.frame(
  etapa = c("Antes de M", "Despues de M"),
  registros = c(
    nrow(p_m),
    nrow(p_m_depurado)
  ),
  validos = c(
    sum(!is.na(p_m$VALUE)),
    sum(!is.na(p_m_depurado$VALUE))
  ),
  estaciones = c(
    n_distinct(p_m$STATION_ID[!is.na(p_m$VALUE)]),
    n_distinct(
      p_m_depurado$STATION_ID[!is.na(p_m_depurado$VALUE)]
    )
  ),
  pares_estacion_anio = c(
    contar_pares(p_m),
    contar_pares(p_m_depurado)
  )
)

print(comparacion_M)

# RESULTADOS DE LIMPIEZA
cat("Meses candidatos a acumulados:",
    nrow(candidatos_M), "\n")

cat("Acumulados tratados:",
    nrow(trazabilidad_M), "\n")

cat("Meses adicionales tratados:",
    nrow(meses_adicionales_M), "\n")

cat("Registros adicionales convertidos en NA:",
    nrow(trazabilidad_adicional_M), "\n")

cat("Total de registros modificados en M:",
    nrow(trazabilidad_total_M), "\n")

cat("Ceros modificados:",
    sum(p_m$VALUE == 0, na.rm = TRUE) -
      sum(p_m_depurado$VALUE == 0, na.rm = TRUE),
    "\n")

# TABLAS DE TRAZABILIDAD
print(tibble::as_tibble(trazabilidad_total_M), n = Inf)

print(tibble::as_tibble(acumulados_pendientes_M), n = Inf)


# DIAS FALTANTES Y AÑOS INCOMPLETOS (1981-2000)
# Construir calendario completo por estacion-anio
datos_N <- p_m_depurado %>%
  select(STATION_ID, YEAR, fecha, VALUE) %>%
  arrange(STATION_ID, fecha)

stopifnot(nrow(datos_N) == nrow(distinct(datos_N, STATION_ID, fecha)),
          !any(is.na(datos_N$fecha)))

fechas_N <- tibble(fecha = seq(as.Date("1981-01-01"), as.Date("2000-12-31"), by = "day"),
                   YEAR = as.integer(format(fecha, "%Y")))

calendario_N <- datos_N %>%
  distinct(STATION_ID, YEAR) %>%
  inner_join(fechas_N, by = "YEAR", relationship = "many-to-many") %>%
  left_join(datos_N %>% mutate(registrado = TRUE) %>%
              select(STATION_ID, fecha, VALUE, registrado),
            by = c("STATION_ID", "fecha")) %>%
  arrange(STATION_ID, YEAR, fecha) %>%
  mutate(fecha_ausente = is.na(registrado),
         valor_na = !fecha_ausente & is.na(VALUE),
         valido = !is.na(VALUE)) %>%
  select(STATION_ID, YEAR, fecha, VALUE, fecha_ausente, valor_na, valido)

# Calcular cobertura, rachas de faltantes y ausencias en los extremos
racha_max <- function(v) {
  r <- rle(!v)
  if (any(r$values)) max(r$lengths[r$values]) else 0
}

borde_max <- function(v) {
  f <- !v
  max(sum(cumprod(f)), sum(cumprod(rev(f))))
}

cobertura_N <- calendario_N %>%
  group_by(STATION_ID, YEAR) %>%
  summarise(dias_esperados = n(),
            dias_validos = sum(valido),
            fechas_ausentes = sum(fecha_ausente),
            valores_na = sum(valor_na),
            max_hueco = racha_max(valido),
            max_borde = borde_max(valido),
            .groups = "drop") %>%
  mutate(dias_faltantes = fechas_ausentes + valores_na,
         cobertura = dias_validos / dias_esperados,
         cobertura_pct = round(100 * cobertura, 2)) %>%
  arrange(STATION_ID, YEAR)

stopifnot(all(cobertura_N$dias_validos + cobertura_N$dias_faltantes ==
                cobertura_N$dias_esperados))

# Resumen general de completitud
resumen_N <- cobertura_N %>%
  summarise(estaciones = n_distinct(STATION_ID),
            pares_estacion_anio = n(),
            dias_esperados = sum(dias_esperados),
            dias_validos = sum(dias_validos),
            fechas_ausentes = sum(fechas_ausentes),
            valores_na = sum(valores_na),
            pares_completos = sum(cobertura == 1))

print(resumen_N)
summary(cobertura_N[, c("cobertura_pct", "max_hueco", "max_borde")])

# Diagnosticar patrones en anios con cobertura inferior al 90 %
diagnostico_N <- datos_N %>%
  filter(!is.na(VALUE)) %>%
  semi_join(cobertura_N %>% filter(cobertura < 0.90),
            by = c("STATION_ID", "YEAR")) %>%
  group_by(STATION_ID, YEAR) %>%
  arrange(fecha, .by_group = TRUE) %>%
  summarise(primera_fecha = min(fecha),
            ultima_fecha = max(fecha),
            meses_observados = n_distinct(format(fecha, "%m")),
            intervalo_tipico = if (n() > 1) median(as.numeric(diff(fecha))) else NA_real_,
            .groups = "drop") %>%
  left_join(cobertura_N %>%
              select(STATION_ID, YEAR, dias_validos, cobertura_pct, max_hueco, max_borde),
            by = c("STATION_ID", "YEAR")) %>%
  mutate(patron = case_when(
    meses_observados == 1 ~ "Registros en un solo mes",
    between(intervalo_tipico, 25, 35) ~ "Posible frecuencia mensual",
    between(intervalo_tipico, 5, 9) ~ "Posible frecuencia semanal",
    max_borde >= 30 ~ "Periodo prolongado sin registros",
    TRUE ~ "Faltantes dispersos u otro patron"
  )) %>%
  arrange(cobertura_pct)

print(diagnostico_N %>% count(patron, name = "pares"))
print(head(diagnostico_N, 20))

# Comparar umbrales en el periodo provisional de entrenamiento (1981-1995)
datos_train_N <- cobertura_N %>% filter(between(YEAR, 1981, 1995))

comparacion_umbrales_N <- tibble(umbral = c(0.85, 0.90, 0.95)) %>%
  mutate(pares_solo_cobertura = sapply(umbral, function(u)
    sum(datos_train_N$cobertura >= u)),
    pares_aptos = sapply(umbral, function(u)
      sum(datos_train_N$cobertura >= u & datos_train_N$max_hueco <= 7 &
            datos_train_N$max_borde <= 7)),
    estaciones_aptas = sapply(umbral, function(u)
      n_distinct(datos_train_N$STATION_ID[datos_train_N$cobertura >= u &
                                            datos_train_N$max_hueco <= 7 & datos_train_N$max_borde <= 7])))

print(comparacion_umbrales_N)

# Seleccionar pares aptos para construir curvas funcionales
cobertura_N <- cobertura_N %>%
  mutate(apto = cobertura >= 0.95 & max_hueco <= 7 & max_borde <= 7)

comparacion_N <- cobertura_N %>%
  summarise(pares_antes = n(),
            pares_aptos = sum(apto),
            pares_excluidos = sum(!apto),
            estaciones_antes = n_distinct(STATION_ID),
            estaciones_despues = n_distinct(STATION_ID[apto]))

print(comparacion_N)

pares_validos_N <- cobertura_N %>% filter(apto) %>% select(STATION_ID, YEAR)
p_n_depurado <- p_m_depurado %>%
  semi_join(pares_validos_N, by = c("STATION_ID", "YEAR"))


# Verificar los datos que dejó el preprocesamiento
names(p_n_depurado)

n_distinct(p_n_depurado$STATION_ID)

nrow(distinct(p_n_depurado, STATION_ID, YEAR))

head(p_n_depurado)

# Calendario real de las 72 pentadas
calendario_pentadas <- expand_grid(
  YEAR = 1981:2000,
  pentada = 1:72
) %>%
  mutate(
    MONTH = (pentada - 1L) %/% 6L + 1L,
    pentada_mes = (pentada - 1L) %% 6L + 1L,
    inicio_mes = as.Date(sprintf("%04d-%02d-01", YEAR, MONTH)),
    siguiente_mes = as.Date(format(inicio_mes + 32, "%Y-%m-01")),
    dias_mes = as.integer(siguiente_mes - inicio_mes),
    dia_inicio = (pentada_mes - 1L) * 5L + 1L,
    dia_fin = ifelse(pentada_mes == 6L, dias_mes, pentada_mes * 5L),
    fecha_inicio = inicio_mes + dia_inicio - 1L,
    fecha_fin = inicio_mes + dia_fin - 1L,
    dias_esperados = as.integer(fecha_fin - fecha_inicio) + 1L
  ) %>%
  select(YEAR, pentada, fecha_inicio, fecha_fin, dias_esperados)

stopifnot(all(
  calendario_pentadas %>%
    group_by(YEAR) %>%
    summarise(dias = sum(dias_esperados), .groups = "drop") %>%
    mutate(esperados = ifelse(
      (YEAR %% 4 == 0 & YEAR %% 100 != 0) | YEAR %% 400 == 0,
      366, 365
    )) %>%
    transmute(ok = dias == esperados) %>%
    pull(ok)
))

# EDA terrestre
datos_tierra <- calendario_N %>%
  semi_join(pares_validos_N, by = c("STATION_ID", "YEAR")) %>%
  mutate(
    STATION_ID = id_txt(STATION_ID),
    MONTH = as.integer(format(fecha, "%m")),
    DAY = as.integer(format(fecha, "%d")),
    pentada_mes = pmin(6L, (DAY - 1L) %/% 5L + 1L),
    pentada = (MONTH - 1L) * 6L + pentada_mes
  )

stopifnot(!anyDuplicated(datos_tierra[c("STATION_ID", "fecha")]))

eda_terrestre <- datos_tierra %>%
  group_by(YEAR) %>%
  summarise(
    estaciones = n_distinct(STATION_ID),
    dias_validos = sum(!is.na(VALUE)),
    dias_faltantes = sum(is.na(VALUE)),
    media_diaria = mean(VALUE, na.rm = TRUE),
    lluvia_maxima = max(VALUE, na.rm = TRUE),
    .groups = "drop"
  )

print(tibble::as_tibble(eda_terrestre), n = Inf)
summary(datos_tierra$VALUE)

print(
  ggplot(datos_tierra %>% filter(!is.na(VALUE)), aes(VALUE)) +
    geom_histogram(bins = 60, fill = "steelblue", color = "white") +
    coord_cartesian(xlim = c(0, 120)) +
    labs(title = "Distribucion de precipitacion terrestre",
         x = "Precipitacion diaria (mm)", y = "Frecuencia") +
    theme_bw()
)

print(
  ggplot(eda_terrestre, aes(YEAR, estaciones)) +
    geom_line(color = "steelblue") +
    geom_point() +
    labs(title = "Estaciones disponibles por anio",
         x = "Anio", y = "Estaciones") +
    theme_bw()
)

# Agregar lluvia terrestre por intervalos
tierra_pentadal <- datos_tierra %>%
  group_by(STATION_ID, YEAR, pentada) %>%
  summarise(
    dias_validos = sum(!is.na(VALUE)),
    lluvia_observada = sum(VALUE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(calendario_pentadas,
            by = c("YEAR", "pentada")) %>%
  mutate(
    dias_faltantes = dias_esperados - dias_validos,
    lluvia_mm = ifelse(dias_faltantes == 0,
                       lluvia_observada, NA_real_),
    curva = paste(STATION_ID, YEAR, sep = "_")
  )

stopifnot(
  nrow(tierra_pentadal) == 72 * nrow(pares_validos_N),
  all(tierra_pentadal$dias_faltantes >= 0)
)

eda_pentadas <- tierra_pentadal %>%
  group_by(pentada) %>%
  summarise(
    media_mm = mean(lluvia_mm, na.rm = TRUE),
    mediana_mm = median(lluvia_mm, na.rm = TRUE),
    desviacion = sd(lluvia_mm, na.rm = TRUE),
    cobertura = 100 * mean(!is.na(lluvia_mm)),
    .groups = "drop"
  )


# Matriz terrestre
tabla_Y <- tierra_pentadal %>%
  select(pentada, curva, lluvia_mm) %>%
  pivot_wider(names_from = curva, values_from = lluvia_mm) %>%
  arrange(pentada)

Y_terrestre_total <- as.matrix(tabla_Y[, -1])
rownames(Y_terrestre_total) <- tabla_Y$pentada

# Datos satelitales
ruta_tif <- "data/ImagenesSatelitePentadal"

if (!dir.exists(ruta_tif)) {
  zip_sat <- list.files("data",
                        pattern = "^ImagenesSatelitePentadal.*\\.zip$",
                        full.names = TRUE)
  if (length(zip_sat) != 1)
    stop("No se encontro la carpeta ni un ZIP satelital unico")
  ruta_tif <- tempfile("satelite_")
  dir.create(ruta_tif)
  unzip(zip_sat, exdir = ruta_tif)
}

archivos_tif <- list.files(
  ruta_tif, pattern = "^ppt[0-9]{6}\\.tif$",
  full.names = TRUE, recursive = TRUE
)

if (length(archivos_tif) == 0)
  stop("No se encontraron imagenes pptYYYYPP.tif")

capas <- data.frame(
  archivo = archivos_tif,
  nombre = basename(archivos_tif)
) %>%
  mutate(
    YEAR = as.integer(substr(nombre, 4, 7)),
    pentada = as.integer(substr(nombre, 8, 9))
  ) %>%
  filter(between(YEAR, 1981, 2000)) %>%
  arrange(YEAR, pentada) %>%
  left_join(calendario_pentadas, by = c("YEAR", "pentada"))

stopifnot(
  !anyDuplicated(capas[c("YEAR", "pentada")]),
  all(capas$pentada %in% 1:72),
  !anyNA(capas$fecha_inicio)
)

disponibilidad_sat <- capas %>%
  count(YEAR, name = "imagenes") %>%
  mutate(faltantes = 72 - imagenes)

print(as_tibble(disponibilidad_sat), n = Inf)


# Cargar capas del periodo de desarrollo 1981-2000
capas_modelo <- capas %>%
  filter(between(YEAR, 1981, 2000))

p_sat <- rast(capas_modelo$archivo)

# EDA espacial satelital
print(p_sat[[1]])


# Metadatos y ubicaciones
metadatos_total <- tierra_pentadal %>%
  distinct(STATION_ID, YEAR, curva) %>%
  left_join(info_estaciones, by = "STATION_ID") %>%
  arrange(STATION_ID, YEAR)

estaciones_geo <- metadatos_total %>%
  filter(!is.na(longitud), !is.na(latitud)) %>%
  distinct(STATION_ID, longitud, latitud)

if (nrow(estaciones_geo) == 0)
  stop("No hay estaciones con coordenadas")

puntos <- vect(estaciones_geo,
               geom = c("longitud", "latitud"),
               crs = "EPSG:4326")

puntos <- project(puntos, crs(p_sat))

estaciones_geo$pixel_id <- cellFromXY(
  p_sat[[1]], crds(puntos)
)

# Extraer precipitacion por pixel
extraccion <- terra::extract(p_sat, puntos, ID = FALSE)

stopifnot(
  nrow(extraccion) == nrow(estaciones_geo),
  ncol(extraccion) == nrow(capas_modelo)
)

sat_long <- tibble(
  STATION_ID = rep(estaciones_geo$STATION_ID,
                   each = nrow(capas_modelo)),
  YEAR = rep(capas_modelo$YEAR,
             times = nrow(estaciones_geo)),
  pentada = rep(capas_modelo$pentada,
                times = nrow(estaciones_geo)),
  lluvia_sat_mm = as.vector(t(as.matrix(extraccion)))
) %>%
  mutate(curva = paste(STATION_ID, YEAR, sep = "_"))

# EDA satelital sobre pixeles de estaciones
eda_satelite <- sat_long %>%
  group_by(YEAR, pentada) %>%
  summarise(
    media_mm = mean(lluvia_sat_mm, na.rm = TRUE),
    mediana_mm = median(lluvia_sat_mm, na.rm = TRUE),
    cobertura = 100 * mean(!is.na(lluvia_sat_mm)),
    .groups = "drop"
  )

print(summary(eda_satelite$media_mm))


# Emparejar estaciones y pixeles
metadatos_total <- metadatos_total %>%
  left_join(
    estaciones_geo %>% select(STATION_ID, pixel_id),
    by = "STATION_ID"
  )

pixeles_compartidos <- estaciones_geo %>%
  filter(!is.na(pixel_id)) %>%
  group_by(pixel_id) %>%
  summarise(
    estaciones = n_distinct(STATION_ID),
    .groups = "drop"
  ) %>%
  filter(estaciones > 1)

conteo_sat <- sat_long %>%
  group_by(curva) %>%
  summarise(
    pentadas = n_distinct(pentada),
    valores_validos = sum(!is.na(lluvia_sat_mm)),
    .groups = "drop"
  )

metadatos_total <- metadatos_total %>%
  left_join(conteo_sat, by = "curva") %>%
  mutate(
    estado = case_when(
      is.na(longitud) | is.na(latitud) ~ "Sin coordenadas",
      is.na(pixel_id) ~ "Fuera del raster",
      is.na(pentadas) | pentadas != 72 ~ "Faltan capas",
      valores_validos != 72 ~ "Pixel incompleto",
      TRUE ~ "Emparejado"
    )
  )

print(metadatos_total %>% count(estado))
print(pixeles_compartidos)

# Matrices finales alineadas
curvas_finales <- metadatos_total %>%
  filter(estado == "Emparejado") %>%
  pull(curva)

if (length(curvas_finales) == 0)
  stop("No hay curvas emparejadas")

Y_est <- Y_terrestre_total[, curvas_finales, drop = FALSE]

tabla_X <- sat_long %>%
  filter(curva %in% curvas_finales) %>%
  select(pentada, curva, lluvia_sat_mm) %>%
  pivot_wider(names_from = curva, values_from = lluvia_sat_mm) %>%
  arrange(pentada)

X_sat <- as.matrix(
  tabla_X[match(1:72, tabla_X$pentada),
          curvas_finales, drop = FALSE]
)

rownames(X_sat) <- as.character(1:72)

metadatos <- metadatos_total %>%
  filter(curva %in% curvas_finales) %>%
  slice(match(curvas_finales, curva))

stopifnot(
  identical(dim(Y_est), dim(X_sat)),
  identical(colnames(Y_est), colnames(X_sat)),
  identical(rownames(Y_est), rownames(X_sat)),
  identical(metadatos$curva, colnames(Y_est)),
  !anyNA(X_sat)
)

# EDA comparativo de ambas fuentes
comparacion <- tibble(
  pentada = rep(1:72, times = ncol(Y_est)),
  tierra = as.vector(Y_est),
  satelite = as.vector(X_sat)
) %>%
  filter(!is.na(tierra), !is.na(satelite)) %>%
  group_by(pentada) %>%
  summarise(
    tierra = mean(tierra),
    satelite = mean(satelite),
    .groups = "drop"
  ) %>%
  pivot_longer(c(tierra, satelite),
               names_to = "fuente",
               values_to = "precipitacion")

print(
  ggplot(comparacion,
         aes(pentada, precipitacion, color = fuente)) +
    geom_line(linewidth = 1) +
    labs(title = "Precipitacion media por fuente",
         x = "Pentada", y = "Precipitacion (mm)") +
    theme_bw()
)

# Resultados finales
cat("Matriz terrestre:", dim(Y_est), "\n")
cat("Matriz satelital:", dim(X_sat), "\n")
cat("Faltantes terrestres:", sum(is.na(Y_est)), "\n")
cat("Curvas con altitud:", sum(!is.na(metadatos$altitud)), "\n")

head(Y_est[, 1:min(5, ncol(Y_est)), drop = FALSE])
head(X_sat[, 1:min(5, ncol(X_sat)), drop = FALSE])
head(metadatos)


# CONSTRUCCION Y VERIFICACION DE MATRICES
# 1. Dominio temporal: 72 pentadas anuales
tiempos <- 1:72

# 2. Transponer matrices
# Filas: curvas estacion-anio
# Columnas: pentadas

Y_matriz <- t(Y_est)
X_matriz <- t(X_sat)

colnames(Y_matriz) <- paste0("P", sprintf("%02d", tiempos))
colnames(X_matriz) <- paste0("P", sprintf("%02d", tiempos))

# 3. Verificar dimensiones
cat("Dimensiones matriz terrestre:\n")
print(dim(Y_matriz))

cat("Dimensiones matriz satelital:\n")
print(dim(X_matriz))

# 4. Verificar correspondencia de curvas
stopifnot(
  identical(dim(Y_matriz), dim(X_matriz)),
  identical(rownames(Y_matriz), rownames(X_matriz)),
  identical(colnames(Y_matriz), colnames(X_matriz)),
  ncol(Y_matriz) == length(tiempos),
  identical(rownames(Y_matriz), metadatos$curva)
)

# 5. Verificar valores faltantes
cat("NA terrestres:", sum(is.na(Y_matriz)), "\n")
cat("NA satelitales:", sum(is.na(X_matriz)), "\n")

# 6. Separar entrenamiento y validacion por anio
indices_train <- which(
  metadatos$YEAR >= 1981 & metadatos$YEAR <= 1995
)

indices_val <- which(
  metadatos$YEAR >= 1996 & metadatos$YEAR <= 2000
)

Y_train <- Y_matriz[indices_train, , drop = FALSE]
X_train <- X_matriz[indices_train, , drop = FALSE]

Y_val <- Y_matriz[indices_val, , drop = FALSE]
X_val <- X_matriz[indices_val, , drop = FALSE]

# 7. Verificar separacion temporal
stopifnot(
  all(metadatos$YEAR[indices_train] <= 1995),
  all(metadatos$YEAR[indices_val] >= 1996),
  length(intersect(indices_train, indices_val)) == 0,
  nrow(Y_train) == nrow(X_train),
  nrow(Y_val) == nrow(X_val)
)

# 8. Resumen
cat("\nRESUMEN DE MATRICES\n")
cat("Curvas totales:", nrow(Y_matriz), "\n")
cat("Pentadas:", ncol(Y_matriz), "\n")
cat("Curvas entrenamiento:", nrow(Y_train), "\n")
cat("Curvas validacion:", nrow(Y_val), "\n")
cat("Faltantes terrestres:", sum(is.na(Y_matriz)), "\n")

# 9. Guardar resultados en formato RDS
dir.create("Resultados", showWarnings = FALSE)

saveRDS(
  list(
    tiempos = tiempos,
    Y_matriz = Y_matriz,
    X_matriz = X_matriz,
    metadatos = metadatos,
    Y_train = Y_train,
    X_train = X_train,
    Y_val = Y_val,
    X_val = X_val
  ),
  "Resultados/matrices_funcionales.rds"
)

cat("Matrices guardadas en Resultados/matrices_funcionales.rds\n")


matrices <- readRDS("Resultados/matrices_funcionales.rds")

Y_matriz <- matrices$Y_matriz
X_matriz <- matrices$X_matriz
metadatos <- matrices$metadatos


# DIAGNOSTICO DE VALORES FALTANTES EN LAS MATRICES
# Verificar alineacion
stopifnot(
  identical(dim(Y_matriz), dim(X_matriz)),
  identical(rownames(Y_matriz), rownames(X_matriz)),
  identical(rownames(Y_matriz), metadatos$curva),
  ncol(Y_matriz) == 72
)

# Contar valores faltantes por curva
faltantes_curva <- rowSums(is.na(Y_matriz))
observados_curva <- rowSums(!is.na(Y_matriz))

# Resumen general
diagnostico_faltantes <- data.frame(
  indicador = c(
    "Total de curvas",
    "Curvas completas",
    "Curvas con al menos un NA",
    "Total de valores NA",
    "Maximo de NA en una curva",
    "Minimo de pentadas observadas",
    "NA satelitales"
  ),
  cantidad = c(
    nrow(Y_matriz),
    sum(faltantes_curva == 0),
    sum(faltantes_curva > 0),
    sum(faltantes_curva),
    max(faltantes_curva),
    min(observados_curva),
    sum(is.na(X_matriz))
  )
)

print(diagnostico_faltantes)

# Distribucion de valores faltantes por curva
print(table(faltantes_curva))

# Curvas afectadas y sus respectivos anios
curvas_con_na <- metadatos %>%
  mutate(
    pentadas_observadas = observados_curva,
    pentadas_faltantes = faltantes_curva
  ) %>%
  filter(pentadas_faltantes > 0) %>%
  arrange(desc(pentadas_faltantes))

print(
  curvas_con_na %>%
    select(
      STATION_ID, YEAR, curva,
      pentadas_observadas, pentadas_faltantes
    ),
  n = Inf
)


# SUFICIENCIA DE OBSERVACIONES PARA MINIMOS CUADRADOS
K_evaluar <- c(6, 8, 10, 12, 15, 20, 25, 30)

suficiencia_K <- data.frame(
  K = K_evaluar,
  curvas_con_n_mayor_igual_K = sapply(
    K_evaluar,
    function(k) sum(observados_curva >= k)
  ),
  curvas_con_n_menor_K = sapply(
    K_evaluar,
    function(k) sum(observados_curva < k)
  ),
  minimo_observado = min(observados_curva)
)

print(suficiencia_K)


# VERIFICACION DEL RANGO DE LA MATRIZ B-SPLINE
# 1. Dominio y valores candidatos
tiempos <- 1:72
K_evaluar <- c(6, 8, 10, 12, 15, 20, 25, 30)
orden <- 4  # B-splines cubicas


# 2. CONSTRUCCION DE LA MATRIZ B-SPLINE CON FDA

# La libreria fda debe estar cargada
library(fda)

# Funcion para construir y evaluar bases B-spline
construir_bspline <- function(t, K, orden = 4) {
  
  if (K < orden) {
    stop("K debe ser mayor o igual al orden.")
  }
  
  # Dominio temporal comun: 72 pentadas
  rangval <- range(tiempos)
  
  # Crear base B-spline
  BasePrec <- create.bspline.basis(
    rangeval = rangval,
    nbasis = K,
    norder = orden
  )
  
  # Evaluar las funciones base en los tiempos indicados
  Phi <- eval.basis(t, BasePrec)
  
  colnames(Phi) <- paste0("B", seq_len(K))
  
  return(Phi)
}

# 3. COMPROBAR LA BASE COMPLETA
Phi_ejemplo <- construir_bspline(
  tiempos,
  K = 12,
  orden = orden
)

cat("Dimensiones de Phi:\n")
print(dim(Phi_ejemplo))

cat("Rango de Phi completa:\n")
print(qr(Phi_ejemplo, tol = 1e-9)$rank)

# Verificar particion de la unidad
stopifnot(
  ncol(Phi_ejemplo) == 12,
  all(abs(rowSums(Phi_ejemplo) - 1) < 1e-8)
)

# 4. RANGO POR CURVA, USANDO SOLO PENTADAS OBSERVADAS
resultados_rango <- list()

for (K in K_evaluar) {
  
  Phi <- construir_bspline(
    tiempos,
    K = K,
    orden = orden
  )
  
  rangos <- integer(nrow(Y_matriz))
  
  for (i in seq_len(nrow(Y_matriz))) {
    
    observados <- which(!is.na(Y_matriz[i, ]))
    
    Phi_i <- Phi[observados, , drop = FALSE]
    
    rangos[i] <- qr(
      Phi_i,
      tol = 1e-9
    )$rank
  }
  
  resultados_rango[[as.character(K)]] <- data.frame(
    curva = rownames(Y_matriz),
    K = K,
    n_observados = rowSums(!is.na(Y_matriz)),
    rango = rangos,
    rango_completo = rangos == K
  )
}

resultados_rango <- dplyr::bind_rows(
  resultados_rango
)

# 5. RESUMEN GENERAL
resumen_rango <- resultados_rango %>%
  dplyr::group_by(K) %>%
  dplyr::summarise(
    curvas_evaluadas = dplyr::n(),
    rango_completo = sum(rango_completo),
    rango_incompleto = sum(!rango_completo),
    rango_minimo = min(rango),
    .groups = "drop"
  )

print(resumen_rango)

# 6. IDENTIFICAR CURVAS PROBLEMATICAS
curvas_rango_incompleto <- resultados_rango %>%
  dplyr::filter(!rango_completo) %>%
  dplyr::arrange(K, rango)

cat("Casos con rango incompleto:",
    nrow(curvas_rango_incompleto), "\n")

print(
  tibble::as_tibble(curvas_rango_incompleto),
  n = Inf
)

# 7. VERIFICACION FINAL
if (all(resultados_rango$rango_completo)) {
  
  cat(
    "\nTodas las curvas tienen rango completo",
    "para los valores de K evaluados.\n"
  )
  
} else {
  
  cat(
    "\nExisten curvas con rango incompleto.",
    "Revisar resultados antes del ajuste.\n"
  )
}


# DIAGNOSTICO DE LA CURVA CON RANGO INCOMPLETO
curva_revision <- "5301400104_1988"

# Localizar la curva
i <- match(curva_revision, rownames(Y_matriz))

stopifnot(!is.na(i))

# Identificar pentadas observadas y faltantes
observadas <- which(!is.na(Y_matriz[i, ]))
faltantes <- which(is.na(Y_matriz[i, ]))

cat("Pentadas faltantes:\n")
print(faltantes)

# Comparar rangos y condicionamiento
diagnostico_rango <- lapply(c(12, 15, 20, 25, 30), function(K) {
  
  Phi <- construir_bspline(
    tiempos,
    K = K,
    orden = 4
  )
  
  Phi_i <- Phi[observadas, , drop = FALSE]
  
  data.frame(
    K = K,
    observaciones = nrow(Phi_i),
    rango = qr(Phi_i, tol = 1e-9)$rank,
    columnas_nulas = sum(colSums(abs(Phi_i)) < 1e-10),
    condicion = kappa(Phi_i, exact = TRUE)
  )
})

diagnostico_rango <- dplyr::bind_rows(diagnostico_rango)

print(diagnostico_rango)

# Verificacion global de los conteos
resultados_rango %>%
  dplyr::group_by(K) %>%
  dplyr::summarise(
    completas = sum(rango_completo, na.rm = TRUE),
    incompletas = sum(rango_completo == FALSE, na.rm = TRUE),
    sin_clasificar = sum(is.na(rango_completo)),
    .groups = "drop"
  ) %>%
  print()


# ============================================================
# CONSTRUCCION Y DIAGNOSTICO DE DATOS FUNCIONALES
# MCO, MCP, MCG Y MCG RESTRINGIDO
# ============================================================

# 1. DATOS Y PARAMETROS
stopifnot(identical(dim(Y_train), dim(X_train)),
          identical(rownames(Y_train), rownames(X_train)),
          identical(colnames(Y_train), colnames(X_train)),
          ncol(Y_train) == length(tiempos),
          exists("construir_bspline"), requireNamespace("fda", quietly = TRUE))

fuentes <- list(terrestre = Y_train, satelital = X_train)
K_evaluar <- c(6, 8, 10, 12, 15, 20, 25, 30)
orden <- 4
umbral_condicion <- 1e4
K_demo <- 12
K_mostrar <- c(6, 12, 20)
K_diagnostico <- c(6, 12, 20)
tolerancia <- 1e-8
malla <- seq(min(tiempos), max(tiempos), length.out = 721)


# 2. CONSTRUCCION FUNCIONAL POR MCO
# c = (Phi'Phi)^(-1) Phi'y
# C = (Phi'Phi)^(-1) Phi' ;  S = Phi C

ajustar_MCO <- function(M, K) {
  Phi <- construir_bspline(tiempos, K, orden)
  n <- nrow(M)
  coeficientes <- matrix(NA_real_, n, K,
                         dimnames = list(rownames(M), colnames(Phi)))
  predicciones <- residuos <- diagonal_S <- matrix(NA_real_, n, ncol(M),
                                                   dimnames = dimnames(M))
  metricas <- data.frame(curva = rownames(M), K = K,
                         n = rowSums(is.finite(M)), estado = "ok",
                         condicion = NA_real_, SSE = NA_real_, gl = NA_real_)
  
  for (i in seq_len(n)) {
    y <- as.numeric(M[i, ])
    obs <- which(is.finite(y))
    P <- Phi[obs, , drop = FALSE]
    if (length(obs) <= K || qr(P, tol = 1e-9)$rank < K) {
      metricas$estado[i] <- "rango_incompleto"
      next
    }
    condicion <- kappa(P, exact = TRUE)
    metricas$condicion[i] <- condicion
    if (!is.finite(condicion) || condicion > umbral_condicion) {
      metricas$estado[i] <- "mal_condicionada"
      next
    }
    C <- solve(crossprod(P), t(P))
    c_est <- as.vector(C %*% y[obs])
    pred <- as.vector(Phi %*% c_est)
    e <- y[obs] - pred[obs]
    coeficientes[i, ] <- c_est
    predicciones[i, ] <- pred
    residuos[i, obs] <- e
    diagonal_S[i, obs] <- rowSums(P * t(C))
    metricas$SSE[i] <- sum(e^2)
    metricas$gl[i] <- sum(diagonal_S[i, obs])
  }
  list(Phi = Phi, coeficientes = coeficientes, predicciones = predicciones,
       residuos = residuos, diagonal_S = diagonal_S, metricas = metricas)
}

resultados_MCO <- lapply(fuentes, function(M)
  setNames(lapply(K_evaluar, function(K) ajustar_MCO(M, K)), K_evaluar))

metricas_MCO <- do.call(rbind, lapply(names(fuentes), function(fuente)
  do.call(rbind, lapply(resultados_MCO[[fuente]], function(a)
    data.frame(fuente = fuente, a$metricas)))))
rownames(metricas_MCO) <- NULL
cat("\nAJUSTES MCO POR FUENTE Y K\n")
print(with(metricas_MCO, table(fuente, K, estado)))

# 2.1 VERIFICACION MATRICIAL DE C Y S
id_demo <- which(rowSums(is.finite(Y_train)) == ncol(Y_train))[1]
stopifnot(length(id_demo) == 1, !is.na(id_demo))
Phi <- resultados_MCO$terrestre[[as.character(K_demo)]]$Phi
y <- as.numeric(Y_train[id_demo, ])
C <- solve(crossprod(Phi), t(Phi))
y2cMap <- C  # Nombre utilizado en las clases de regresion funcional
S <- Phi %*% C
c_est <- as.vector(C %*% y)
y_est <- as.vector(Phi %*% c_est)
verificacion_S <- data.frame(
  curva = rownames(Y_train)[id_demo], K = K_demo,
  filas_C = nrow(C), columnas_C = ncol(C),
  filas_S = nrow(S), columnas_S = ncol(S),
  simetria = max(abs(S - t(S))),
  idempotencia = max(abs(S %*% S - S)),
  traza = sum(diag(S)), diferencia_K = abs(sum(diag(S)) - K_demo),
  diferencia_ajuste = max(abs(as.vector(S %*% y) - y_est)))
cat("\nPROPIEDADES DE LA MATRIZ SUAVIZADORA MCO\n")
print(verificacion_S)
stopifnot(verificacion_S$simetria < 1e-7,
          verificacion_S$idempotencia < 1e-7,
          verificacion_S$diferencia_K < 1e-7,
          verificacion_S$diferencia_ajuste < 1e-7)

# 2.2 COMPARACION CON fda::smooth.basis
base_demo <- fda::create.bspline.basis(range(tiempos), nbasis = K_demo,
                                       norder = orden)
comparacion_fda <- do.call(rbind, lapply(names(fuentes), function(fuente) {
  M <- fuentes[[fuente]]
  id <- which(rowSums(is.finite(M)) == ncol(M))[1]
  a <- fda::smooth.basis(tiempos, as.numeric(M[id, ]), base_demo)
  manual <- resultados_MCO[[fuente]][[as.character(K_demo)]]
  data.frame(fuente = fuente,
             diferencia_coef = max(abs(manual$coeficientes[id, ] - as.vector(a$fd$coefs))),
             diferencia_pred = max(abs(manual$predicciones[id, ] -
                                         as.vector(fda::eval.fd(tiempos, a$fd)))))
}))
cat("\nCOMPARACION MCO Y SMOOTH.BASIS\n")
print(comparacion_fda)

# 2.3 GRAFICA: EFECTO DEL NUMERO DE BASES
if (interactive()) {
  op <- par(mfrow = c(1, 2))
  for (fuente in names(fuentes)) {
    M <- fuentes[[fuente]]
    plot(tiempos, M[id_demo, ], pch = 16, cex = 0.5,
         xlab = "Pentada", ylab = "Precipitacion (mm)",
         main = paste("Ajuste funcional -", fuente))
    for (j in seq_along(K_mostrar))
      lines(tiempos,
            resultados_MCO[[fuente]][[as.character(K_mostrar[j])]]$predicciones[id_demo, ],
            col = j + 1, lwd = 2)
    legend("topright", legend = paste("K =", K_mostrar),
           col = seq_along(K_mostrar) + 1, lty = 1, bty = "n")
  }
  par(op)
}

# 3. DIAGNOSTICOS DE MCO: HETEROCEDASTICIDAD Y CORRELACION
# Se examinan los residuos MCO, no se asume que sean los errores
# originales: el suavizado tambien puede inducir correlaciones.
lluviosas <- c(13:30, 49:66)  # Indices de pentadas usados en el analisis
otras <- setdiff(seq_along(tiempos), lluviosas)

diagnostico_hetero <- do.call(rbind, lapply(names(fuentes), function(fuente) {
  do.call(rbind, lapply(K_diagnostico, function(K) {
    R <- resultados_MCO[[fuente]][[as.character(K)]]$residuos
    var_p <- apply(R, 2, var, na.rm = TRUE)
    v_lluviosa <- apply(R[, lluviosas, drop = FALSE], 1, var, na.rm = TRUE)
    v_otra <- apply(R[, otras, drop = FALSE], 1, var, na.rm = TRUE)
    validos <- is.finite(v_lluviosa) & is.finite(v_otra) & v_otra > 0
    razon <- v_lluviosa[validos] / v_otra[validos]
    data.frame(fuente = fuente, K = K, curvas = sum(validos),
               var_min = min(var_p, na.rm = TRUE),
               var_max = max(var_p, na.rm = TRUE),
               razon_mediana = median(razon),
               pct_mayor_var_lluviosa = 100 * mean(razon > 1))
  }))
}))
rownames(diagnostico_hetero) <- NULL
cat("\nHETEROCEDASTICIDAD: VARIANZA RESIDUAL POR TEMPORADA\n")
print(diagnostico_hetero)

# 3.1 GRAFICA: VARIANZA RESIDUAL POR PENTADA (K = 12)
if (interactive()) {
  op <- par(mfrow = c(1, 2))
  for (fuente in names(fuentes)) {
    R <- resultados_MCO[[fuente]][[as.character(K_demo)]]$residuos
    plot(tiempos, apply(R, 2, var, na.rm = TRUE), type = "l", lwd = 2,
         main = paste("Varianza residual MCO -", fuente),
         xlab = "Pentada", ylab = "Varianza residual")
  }
  par(op)
}

# 3.2 GRAFICA: CORRELACION RESIDUAL ENTRE PENTADAS (K = 12)
if (interactive()) {
  op <- par(mfrow = c(1, 2))
  for (fuente in names(fuentes)) {
    R <- resultados_MCO[[fuente]][[as.character(K_demo)]]$residuos
    Cor <- cor(R, use = "pairwise.complete.obs")
    contour(tiempos, tiempos, Cor, nlevels = 10,
            xlab = "Pentada s", ylab = "Pentada t",
            main = paste("Correlacion residual MCO -", fuente))
  }
  par(op)
}

# 4. DIAGNOSTICO DE PRECIPITACIONES NEGATIVAS EN MCO
diagnostico_negativos <- do.call(rbind, lapply(names(fuentes), function(fuente) {
  do.call(rbind, lapply(K_evaluar, function(K) {
    P <- resultados_MCO[[fuente]][[as.character(K)]]$predicciones
    data.frame(fuente = fuente, K = K,
               negativos = sum(P < -tolerancia, na.rm = TRUE),
               curvas_afectadas = sum(rowSums(P < -tolerancia, na.rm = TRUE) > 0),
               minimo = min(P, na.rm = TRUE))
  }))
}))
rownames(diagnostico_negativos) <- NULL
cat("\nVALORES NEGATIVOS EN MCO ANTES DE RESTRINGIR\n")
print(diagnostico_negativos)

# 5. COMPARACION MCO, MCP Y MCG
# MCP considera varianzas distintas (diagonal); MCG tambien
# incorpora covarianzas. Se comparan en pentadas reservadas.
set.seed(123)
estaciones <- sub("_[0-9]{4}$", "", rownames(Y_train))
est_ref <- sample(unique(estaciones), floor(0.7 * length(unique(estaciones))))
idx_ref <- estaciones %in% est_ref
idx_eval <- !idx_ref
pentadas_prueba <- seq(5, 72, by = 5)
alpha <- 0.5  # Regularizacion exploratoria de la covarianza

# Ajuste sin restricciones; Sigma = NULL equivale a MCO.
ajustar_metodo <- function(y, Phi, obs, Sigma = NULL) {
  tryCatch({
    P <- Phi[obs, , drop = FALSE]
    z <- y[obs]
    if (!is.null(Sigma)) {
      L <- t(chol(Sigma[obs, obs, drop = FALSE]))
      z <- forwardsolve(L, z)
      P <- forwardsolve(L, P)
    }
    if (nrow(P) <= ncol(P) || qr(P, tol = 1e-9)$rank < ncol(P)) return(NULL)
    condicion <- kappa(P, exact = TRUE)
    if (!is.finite(condicion) || condicion > umbral_condicion) return(NULL)
    coef <- as.vector(solve(crossprod(P), crossprod(P, z)))
    list(coef = coef, pred = as.vector(Phi %*% coef))
  }, error = function(e) NULL)
}

# MCG restringido: minimize ||L^(-1)(y-Phi*c)||^2, c >= 0.
# Una sola version definitiva: optim y alternativa nlminb.
ajustar_MCG_positivo <- function(y, Phi, obs, Sigma) {
  K <- ncol(Phi)
  if (length(obs) <= K) return(list(estado = "observaciones_insuficientes"))
  A <- tryCatch({
    L <- t(chol(Sigma[obs, obs, drop = FALSE]))
    list(P = forwardsolve(L, Phi[obs, , drop = FALSE]),
         z = forwardsolve(L, y[obs]))
  }, error = function(e) NULL)
  if (is.null(A)) return(list(estado = "error_covarianza"))
  Pw <- A$P
  zw <- A$z
  if (qr(Pw, tol = 1e-9)$rank < K)
    return(list(estado = "rango_incompleto"))
  condicion <- kappa(Pw, exact = TRUE)
  if (!is.finite(condicion) || condicion > umbral_condicion)
    return(list(estado = "mal_condicionada"))
  escala <- max(1, sum(zw^2))
  objetivo <- function(c) sum((zw - Pw %*% c)^2) / escala
  gradiente <- function(c) as.vector(2 * crossprod(Pw, Pw %*% c - zw) / escala)
  inicio <- tryCatch(pmax(0, as.vector(qr.solve(Pw, zw))),
                     error = function(e) NULL)
  if (is.null(inicio)) return(list(estado = "error_inicio"))
  ajuste <- tryCatch(
    optim(inicio, objetivo, gradiente, method = "L-BFGS-B",
          lower = rep(0, K), control = list(maxit = 3000)),
    error = function(e) NULL)
  if (is.null(ajuste) || ajuste$convergence != 0 ||
      any(!is.finite(ajuste$par))) {
    ajuste <- tryCatch(
      nlminb(inicio, objetivo, gradiente, lower = rep(0, K),
             control = list(iter.max = 3000, eval.max = 5000)),
      error = function(e) NULL)
  }
  if (is.null(ajuste) || ajuste$convergence != 0 ||
      any(!is.finite(ajuste$par))) return(list(estado = "sin_convergencia"))
  coef <- pmax(0, ajuste$par)
  g <- gradiente(coef)
  g[coef <= 1e-8 & g > 0] <- 0
  if (max(abs(g)) > 1e-4) return(list(estado = "optimalidad_no_verificada"))
  list(estado = "ok", coef = coef,
       pred = as.vector(Phi %*% coef),
       error_ponderado = sum((zw - Pw %*% coef)^2))
}

# 5.1 MATRICES DE COVARIANZA Y VALIDACION DE LOS TRES METODOS
resultados_comparacion <- list()
diagnostico_metodos <- list()
for (fuente in names(fuentes)) {
  M <- fuentes[[fuente]]
  resultados_comparacion[[fuente]] <- list()
  for (K in K_evaluar) {
    clave <- as.character(K)
    a <- resultados_MCO[[fuente]][[clave]]
    R <- a$residuos
    ref_ok <- idx_ref & rowSums(is.finite(R)) == ncol(R)
    R_ref <- R[ref_ok, , drop = FALSE]
    if (nrow(R_ref) <= ncol(R_ref)) stop("Curvas de referencia insuficientes")
    Sigma <- cov(R_ref)
    D <- diag(diag(Sigma), ncol(M))
    eps <- 0.01 * median(diag(Sigma))
    Sigma_MCP <- D + eps * diag(ncol(M))
    Sigma_MCG <- (1 - alpha) * Sigma + alpha * D + eps * diag(ncol(M))
    
    detalle <- list()
    ejemplo <- NULL
    # Negativos de los mismos ajustes validados; evita otro recorrido.
    conteo <- matrix(0, 3, 4, dimnames = list(c("MCO", "MCP", "MCG"),
                                              c("curvas", "negativos", "afectadas", "neg_prueba")))
    minimos <- rep(Inf, 3)
    for (i in which(idx_eval)) {
      y <- as.numeric(M[i, ])
      prueba <- intersect(which(is.finite(y)), pentadas_prueba)
      obs <- setdiff(which(is.finite(y)), prueba)
      if (!length(prueba)) next
      ajustes <- list(MCO = ajustar_metodo(y, a$Phi, obs),
                      MCP = ajustar_metodo(y, a$Phi, obs, Sigma_MCP),
                      MCG = ajustar_metodo(y, a$Phi, obs, Sigma_MCG))
      if (any(vapply(ajustes, is.null, logical(1)))) next
      detalle[[length(detalle) + 1]] <- data.frame(
        curva = rownames(M)[i],
        MCO = mean((y[prueba] - ajustes$MCO$pred[prueba])^2),
        MCP = mean((y[prueba] - ajustes$MCP$pred[prueba])^2),
        MCG = mean((y[prueba] - ajustes$MCG$pred[prueba])^2))
      for (j in seq_along(ajustes)) {
        p <- ajustes[[j]]$pred
        conteo[j, ] <- conteo[j, ] + c(1, sum(p < -tolerancia),
                                       as.integer(any(p < -tolerancia)),
                                       sum(p[prueba] < -tolerancia))
        minimos[j] <- min(minimos[j], p)
      }
      if (is.null(ejemplo))
        ejemplo <- c(list(y = y, prueba = prueba, curva = rownames(M)[i]),
                     lapply(ajustes, function(z) z$pred))
    }
    detalle <- if (length(detalle)) do.call(rbind, detalle) else
      data.frame(curva = character(), MCO = numeric(), MCP = numeric(), MCG = numeric())
    resultados_comparacion[[fuente]][[clave]] <- list(
      Sigma_MCP = Sigma_MCP, Sigma_MCG = Sigma_MCG,
      detalle = detalle, ejemplo = ejemplo)
    diagnostico_metodos[[paste(fuente, K)]] <- data.frame(
      fuente = fuente, K = K, metodo = rownames(conteo),
      curvas = conteo[, "curvas"], negativos = conteo[, "negativos"],
      curvas_afectadas = conteo[, "afectadas"],
      negativos_prueba = conteo[, "neg_prueba"],
      minimo = ifelse(is.finite(minimos), minimos, NA_real_))
  }
}

comparacion_metodos <- do.call(rbind, lapply(names(fuentes), function(fuente)
  do.call(rbind, lapply(K_evaluar, function(K) {
    d <- resultados_comparacion[[fuente]][[as.character(K)]]$detalle
    data.frame(fuente = fuente, K = K, curvas = nrow(d),
               ECM_MCO = if (nrow(d)) mean(d$MCO) else NA_real_,
               ECM_MCP = if (nrow(d)) mean(d$MCP) else NA_real_,
               ECM_MCG = if (nrow(d)) mean(d$MCG) else NA_real_)
  }))))
rownames(comparacion_metodos) <- NULL
comparacion_metodos$mejora_MCP <- with(comparacion_metodos,
                                       100 * (ECM_MCO - ECM_MCP) / ECM_MCO)
comparacion_metodos$mejora_MCG <- with(comparacion_metodos,
                                       100 * (ECM_MCO - ECM_MCG) / ECM_MCO)
cat("\nMCO, MCP Y MCG: ECM EN PENTADAS RESERVADAS Y MEJORA VS MCO\n")
print(comparacion_metodos)

diagnostico_metodos <- do.call(rbind, diagnostico_metodos)
rownames(diagnostico_metodos) <- NULL
diagnostico_metodos$porcentaje_afectadas <- with(diagnostico_metodos,
                                                 100 * curvas_afectadas / curvas)
cat("\nVALORES NEGATIVOS ANTES DE RESTRINGIR: MCO, MCP Y MCG\n")
print(diagnostico_metodos)

# 5.2 PROPIEDADES MATRICIALES DE MCG
Phi <- resultados_MCO$terrestre[[as.character(K_demo)]]$Phi
Sigma <- resultados_comparacion$terrestre[[as.character(K_demo)]]$Sigma_MCG
W <- solve(Sigma)
C_MCG <- solve(crossprod(Phi, W %*% Phi), t(Phi) %*% W)
S_MCG <- Phi %*% C_MCG
verificacion_MCG <- data.frame(
  idempotencia = max(abs(S_MCG %*% S_MCG - S_MCG)),
  traza = sum(diag(S_MCG)),
  simetria_ponderada = max(abs(t(S_MCG) %*% W - W %*% S_MCG)),
  simetria_ordinaria = max(abs(S_MCG - t(S_MCG))))
cat("\nPROPIEDADES DE LA MATRIZ SUAVIZADORA MCG\n")
print(verificacion_MCG)

# 5.3 GRAFICA: MCO, MCP Y MCG (K = 12)
if (interactive()) {
  op <- par(mfrow = c(1, 2))
  for (fuente in names(fuentes)) {
    a <- resultados_comparacion[[fuente]][[as.character(K_demo)]]$ejemplo
    if (is.null(a)) next
    plot(tiempos, a$y, pch = 16, cex = 0.5,
         xlab = "Pentada", ylab = "Precipitacion (mm)",
         main = paste("MCO, MCP y MCG -", fuente))
    lines(tiempos, a$MCO, col = "blue", lwd = 2)
    lines(tiempos, a$MCP, col = "red", lwd = 2)
    lines(tiempos, a$MCG, col = "darkgreen", lwd = 2)
    points(tiempos[a$prueba], a$y[a$prueba], pch = 4, cex = 1.2)
    legend("topright", c("MCO", "MCP", "MCG", "Pentadas reservadas"),
           col = c("blue", "red", "darkgreen", "black"),
           lty = c(1, 1, 1, 0), pch = c(NA, NA, NA, 4), bty = "n")
  }
  par(op)
}

# 6. CONSTRUCCION DE TODAS LAS CURVAS CON MCG RESTRINGIDO
# Coeficientes c >= 0 garantizan x(t) >= 0 porque las B-splines
# son no negativas. No se recortan las predicciones a posteriori.
# Esta restriccion puede aumentar el error: se comprueba con ECM.

resultados_MCG_restringido <- list()
for (fuente in names(fuentes)) {
  M <- fuentes[[fuente]]
  resultados_MCG_restringido[[fuente]] <- list()
  for (K in K_evaluar) {
    clave <- as.character(K)
    Phi <- resultados_MCO[[fuente]][[clave]]$Phi
    Sigma <- resultados_comparacion[[fuente]][[clave]]$Sigma_MCG
    coeficientes <- matrix(NA_real_, nrow(M), K,
                           dimnames = list(rownames(M), colnames(Phi)))
    predicciones <- pred_MCG <- matrix(NA_real_, nrow(M), ncol(M),
                                       dimnames = dimnames(M))
    metricas <- data.frame(curva = rownames(M), estado = "no_evaluado",
                           minimo_MCG = NA_real_, minimo_restringido = NA_real_,
                           negativos_MCG = NA_integer_, negativos_restringido = NA_integer_,
                           ECM_MCG = NA_real_, ECM_restringido = NA_real_)
    for (i in seq_len(nrow(M))) {
      y <- as.numeric(M[i, ])
      obs <- which(is.finite(y))
      mcg <- ajustar_metodo(y, Phi, obs, Sigma)
      positivo <- ajustar_MCG_positivo(y, Phi, obs, Sigma)
      metricas$estado[i] <- positivo$estado
      # Se conserva el restringido aunque falle el MCG ordinario.
      if (positivo$estado == "ok") {
        coeficientes[i, ] <- positivo$coef
        predicciones[i, ] <- positivo$pred
        metricas$minimo_restringido[i] <- min(positivo$pred)
        metricas$negativos_restringido[i] <- sum(positivo$pred < -tolerancia)
        metricas$ECM_restringido[i] <- mean((y[obs] - positivo$pred[obs])^2)
      }
      if (!is.null(mcg)) {
        pred_MCG[i, ] <- mcg$pred
        metricas$minimo_MCG[i] <- min(mcg$pred)
        metricas$negativos_MCG[i] <- sum(mcg$pred < -tolerancia)
        metricas$ECM_MCG[i] <- mean((y[obs] - mcg$pred[obs])^2)
      }
    }
    resultados_MCG_restringido[[fuente]][[clave]] <- list(
      Phi = Phi, Sigma = Sigma, coeficientes = coeficientes,
      predicciones = predicciones, pred_MCG = pred_MCG, metricas = metricas)
  }
}

# 6.1 VERIFICACION DE POSITIVIDAD PARA TODOS LOS K
verificacion_final <- do.call(rbind, lapply(names(fuentes), function(fuente)
  do.call(rbind, lapply(K_evaluar, function(K) {
    a <- resultados_MCG_restringido[[fuente]][[as.character(K)]]
    C <- a$coeficientes
    P <- a$predicciones
    data.frame(fuente = fuente, K = K,
               ajustadas = sum(a$metricas$estado == "ok"),
               pendientes = sum(a$metricas$estado != "ok"),
               coef_negativos = sum(C < -tolerancia, na.rm = TRUE),
               pred_negativas = sum(P < -tolerancia, na.rm = TRUE),
               minimo = min(P, na.rm = TRUE))
  }))))
rownames(verificacion_final) <- NULL
cat("\nVERIFICACION FINAL: COBERTURA Y NO NEGATIVIDAD\n")
print(verificacion_final)

curvas_pendientes <- do.call(rbind, lapply(names(fuentes), function(fuente)
  do.call(rbind, lapply(K_evaluar, function(K) {
    m <- resultados_MCG_restringido[[fuente]][[as.character(K)]]$metricas
    p <- m[m$estado != "ok", c("curva", "estado"), drop = FALSE]
    if (!nrow(p)) return(NULL)
    data.frame(fuente = fuente, K = K, p)
  }))))
cat("\nCASOS PENDIENTES DE AJUSTE RESTRINGIDO\n")
if (is.null(curvas_pendientes)) cat("Ninguno\n") else {
  rownames(curvas_pendientes) <- NULL
  print(curvas_pendientes)
}

# 6.2 VALIDACION MCG VS MCG RESTRINGIDO
# Solo estaciones de evaluacion; mismas pentadas reservadas.
# Incluye los ajustes que antes se recuperaban en otro bloque.
comparacion_validacion <- list()
for (fuente in names(fuentes)) {
  M <- fuentes[[fuente]]
  for (K in K_evaluar) {
    clave <- as.character(K)
    a <- resultados_MCG_restringido[[fuente]][[clave]]
    ecm <- list()
    for (i in which(idx_eval)) {
      y <- as.numeric(M[i, ])
      prueba <- intersect(which(is.finite(y)), pentadas_prueba)
      obs <- setdiff(which(is.finite(y)), prueba)
      if (!length(prueba)) next
      mcg <- ajustar_metodo(y, a$Phi, obs, a$Sigma)
      positivo <- ajustar_MCG_positivo(y, a$Phi, obs, a$Sigma)
      if (is.null(mcg) || positivo$estado != "ok") next
      ecm[[length(ecm) + 1]] <- c(
        MCG = mean((y[prueba] - mcg$pred[prueba])^2),
        restringido = mean((y[prueba] - positivo$pred[prueba])^2))
    }
    d <- if (length(ecm)) do.call(rbind, ecm) else matrix(numeric(), 0, 2)
    comparacion_validacion[[paste(fuente, K)]] <- data.frame(
      fuente = fuente, K = K, curvas = nrow(d),
      ECM_MCG = if (nrow(d)) mean(d[, 1]) else NA_real_,
      ECM_restringido = if (nrow(d)) mean(d[, 2]) else NA_real_)
  }
}
comparacion_validacion <- do.call(rbind, comparacion_validacion)
rownames(comparacion_validacion) <- NULL
comparacion_validacion$mejora_pct <- with(comparacion_validacion,
                                          100 * (ECM_MCG - ECM_restringido) / ECM_MCG)
cat("\nMCG VS MCG RESTRINGIDO: ECM EN PENTADAS RESERVADAS\n")
print(comparacion_validacion)

# 6.3 GRAFICA: MCG VS MCG RESTRINGIDO (K = 12)
if (interactive()) {
  op <- par(mfrow = c(1, 2))
  for (fuente in names(fuentes)) {
    M <- fuentes[[fuente]]
    a <- resultados_MCG_restringido[[fuente]][[as.character(K_demo)]]
    candidatos <- which(a$metricas$negativos_MCG > 0 &
                          a$metricas$estado == "ok")
    i <- if (length(candidatos)) candidatos[1] else
      which(a$metricas$estado == "ok")[1]
    if (is.na(i)) next
    base <- fda::create.bspline.basis(range(tiempos), K_demo, orden)
    pred_fina <- as.vector(fda::eval.basis(malla, base) %*% a$coeficientes[i, ])
    plot(tiempos, M[i, ], pch = 16, cex = 0.5,
         xlab = "Pentada", ylab = "Precipitacion (mm)",
         main = paste("MCG vs MCG restringido -", fuente))
    lines(tiempos, a$pred_MCG[i, ], col = "blue", lwd = 2)
    lines(malla, pred_fina, col = "red", lwd = 2)
    abline(h = 0, lty = 2)
    legend("topright", c("MCG", "MCG restringido"),
           col = c("blue", "red"), lty = 1, lwd = 2, bty = "n")
  }
  par(op)
}
