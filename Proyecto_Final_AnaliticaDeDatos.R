####################################################
# Proyecto 3 
# Por: Natalia Paredes, Mariana Cifuentes, John Franco
# y Estefanía Hernández
####################################################
##### LIBRERIAS
#install.packages("fda")
#install.packages("fda.usc")
#install.packages("readxl")
library(fda)
library(fda.usc)
library(readxl)
library(terra)
library(dplyr)
library(ggplot2)
library(tidyterra)
library(geodata)

set.seed(2026)

## DIRECTORIO DE TRABAJO
rm(list=ls())
setwd("C:/Users/maria/Downloads/proyecto_final_analiticaDeDatos")
dir.create("Figuras", showWarnings = FALSE)

##### CARGA DE DATOS: precipitacion satelite y tierra
p_tierra_1930_1979 <- read_excel("data/P_1930_1979.xlsx")
p_tierra_1980_1990 <- read.csv("data/P_1980_1990.csv")
p_tierra_1991_2003 <- read.csv("data/P_1991_2003.csv")

# de 1981 a 2013
p_satelite <- rast(list.files("data/ImagenesSatelitePentadal", 
                              pattern = "ppt.*\\.tif$", full.names = TRUE))

# de 1930 a 2003
p_tierra <- bind_rows(p_tierra_1930_1979, 
                      p_tierra_1980_1990, 
                      p_tierra_1991_2003)

p_tierra <- p_tierra[, -c(2, 7)]
p_tierra$fecha <- as.Date(sprintf("%04d-%02d-%02d", p_tierra$YEAR, p_tierra$MONTH, p_tierra$DAY))
p_tierra <- p_tierra[order(p_tierra$STATION_ID, p_tierra$fecha), ]

# Borde del Valle
departamentos_col <- gadm(country = "COL", level = 1, path = tempdir())  
valle <- departamentos_col[departamentos_col$NAME_1 == "Valle del Cauca", ]

ggplot() + geom_spatraster(data = p_satelite[[1]]) +
  geom_sf(data = valle, fill = NA, color = "black", linewidth = 0.5) +
  labs(title = "Precipitación 01-01-1981") 


# Valores imposibles
summary(p_tierra)

quantile(p_tierra$VALUE, probs = c(0.5, 0.9, 0.95, 0.99, 0.999, 0.99999), na.rm = TRUE)
tapply(p_tierra$VALUE, p_tierra$YEAR, summary)

## 1984: valor 859375
sum(p_tierra$VALUE > 300, na.rm = TRUE)

mayores_300 <- p_tierra[!is.na(p_tierra$VALUE) & p_tierra$VALUE > 300,
                        c("STATION_ID", "YEAR", "MONTH", "DAY", "VALUE")]
mayores_300 <- mayores_300[order(-mayores_300$YEAR), ]
print(mayores_300, n = 300)

## Valores menores al percentil 99.9 (120 mm)
hist(p_tierra$VALUE[p_tierra$VALUE < 120], 
     main = "Precipitación diaria",
     xlab = "Precipitación (mm)", ylab = "Frecuencia")

## Outliers en cada estacion y año
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

## Eliminar los valores aislados > 400
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
####################################################
# ACUMULADOS MENSUALES Y CEROS FALSOS
####################################################

## Datos del periodo disponible (1981-2003)
p_m <- p_tierra %>%
  filter(YEAR >= 1981, YEAR <= 2003) %>%
  arrange(STATION_ID, fecha) %>%
  mutate(fin_mes = format(fecha + 1, "%d") == "01")

cat("Registros del periodo:", nrow(p_m), "\n")
cat("Estaciones:", n_distinct(p_m$STATION_ID), "\n")
cat("Pares estacion-anio:", nrow(distinct(p_m, STATION_ID, YEAR)), "\n")


## Resumen de registros por estacion y mes
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

## Indicadores de posibles acumulados
resumen_mensual <- resumen_mensual %>%
  mutate(proporcion_ceros = ifelse(n_validos > 0, n_ceros / n_validos, NA_real_),
         proporcion_fin_mes = ifelse(total_observado > 0,
                                     valor_fin_mes / total_observado, NA_real_),
         # Regla exploratoria, no clasificacion definitiva
         candidato = !is.na(valor_fin_mes) &
           valor_fin_mes >= 100 &
           proporcion_fin_mes >= 0.80)


## Meses candidatos
candidatos_M <- resumen_mensual %>%
  filter(candidato) %>%
  arrange(STATION_ID, YEAR, MONTH)

cat("Meses candidatos:", nrow(candidatos_M), "\n")

print(candidatos_M %>%
        select(STATION_ID, YEAR, MONTH, n_registros, n_ceros, n_positivos,
               valor_fin_mes, total_observado, proporcion_fin_mes),
      n = 30)

## Patrones que se repiten en varios meses del mismo anio
patrones_M <- candidatos_M %>%
  group_by(STATION_ID, YEAR) %>%
  summarise(meses_sospechosos = n(), .groups = "drop") %>%
  arrange(desc(meses_sospechosos))

print(patrones_M, n = 30)

## Casos conocidos en detalle
caso_1981 <- resumen_mensual %>% filter(STATION_ID == 5421200201, YEAR == 1981)
caso_1986 <- resumen_mensual %>% filter(STATION_ID == 5300000201, YEAR == 1986)

print(caso_1981, n = 12)
print(caso_1986, n = 12)

## Ceros en los meses candidatos
# Son candidatos de revision, NO ceros falsos confirmados
ceros_M <- p_m %>%
  semi_join(candidatos_M, by = c("STATION_ID", "YEAR", "MONTH")) %>%
  filter(VALUE == 0)

cat("Ceros en meses candidatos:", nrow(ceros_M), "\n")

## Conteos del diagnostico
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


## Totales mensuales segun la estructura de registros
# Meses con una sola observacion valida ubicada en el ultimo dia del mes
mensuales_M <- resumen_mensual %>%
  filter(n_registros == 1, n_validos == 1, !is.na(valor_fin_mes)) %>%
  arrange(STATION_ID, YEAR, MONTH)

cat("Meses con un unico registro al final:", nrow(mensuales_M), "\n")

print(mensuales_M %>%
        select(STATION_ID, YEAR, MONTH, n_registros, valor_fin_mes),
      n = 40)

## Patrones repetidos por estacion y anio
patrones_mensuales <- mensuales_M %>%
  group_by(STATION_ID, YEAR) %>%
  summarise(meses_con_un_registro = n(), .groups = "drop") %>%
  arrange(desc(meses_con_un_registro))

print(patrones_mensuales, n = 30)


# Considerar todos los meses con un unico registro al final
alta_sospecha_M <- patrones_mensuales %>%
  filter(meses_con_un_registro >= 1)

cat("Pares estacion-anio con registro unico al final:", nrow(alta_sospecha_M), "\n")

print(alta_sospecha_M, n = 30)

## Estaciones problematicas a traves de los anios
problemas_estacion <- patrones_mensuales %>%
  group_by(STATION_ID) %>%
  summarise(anios_afectados = n_distinct(YEAR),
            meses_sospechosos = sum(meses_con_un_registro),
            anios = paste(sort(unique(YEAR)), collapse = ", "),
            .groups = "drop") %>%
  arrange(desc(meses_sospechosos))

print(problemas_estacion, n = Inf)


## Auditoria de acumulados: fechas exactas de los casos prioritarios
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

print(auditoria_acumulados_M, n = Inf)


## Diagnostico de ceros falsos
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

print(diagnostico_ceros_M, n = Inf)

## Registros cero para revision (no se modifican)
auditoria_ceros_M <- p_m %>%
  semi_join(diagnostico_ceros_M, by = c("STATION_ID", "YEAR", "MONTH")) %>%
  filter(VALUE == 0) %>%
  transmute(STATION_ID, YEAR, MONTH, fecha,
            valor_original = VALUE,
            motivo = "Cero en mes con patron sospechoso",
            decision = "Pendiente de revision")

cat("Ceros candidatos a revision:", nrow(auditoria_ceros_M), "\n")

print(auditoria_ceros_M, n = Inf)

## Dias con lluvia positiva dentro de esos meses
revision_ceros <- p_m %>%
  semi_join(diagnostico_ceros_M, by = c("STATION_ID", "YEAR", "MONTH")) %>%
  filter(!is.na(VALUE), VALUE > 0) %>%
  arrange(STATION_ID, YEAR, MONTH, fecha) %>%
  group_by(STATION_ID, YEAR, MONTH) %>%
  summarise(dias_con_lluvia = n(),
            fechas_lluvia = paste(format(fecha, "%d"), collapse = ", "),
            valores_lluvia = paste(VALUE, collapse = ", "),
            .groups = "drop")

print(revision_ceros, n = Inf)

## Contexto temporal: mes anterior, sospechoso y posterior
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

print(contexto_ceros_M, n = Inf)

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

print(comparacion_historica_M, n = Inf)


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

print(comparacion_ceros_historica_M, n = Inf)


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

print(comparacion_otras_estaciones_M, n = Inf, width = Inf)

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

print(historico_picos_M, n = Inf, width = Inf)

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

print(head(comparacion_todo_cero_M, 30), width = Inf)


# Limpieza final y trazabilidad
p_m_depurado <- p_m

# Los 36 registros con unica observacion al final del mes
acumulados_tratar_M <- auditoria_acumulados_M %>%
  select(STATION_ID, fecha, valor_original, motivo) %>%
  distinct(STATION_ID, fecha, .keep_all = TRUE)

# Verificar que los registros coincidan con la base
stopifnot(nrow(acumulados_tratar_M) == 36)

verificacion_M <- acumulados_tratar_M %>%
  left_join(p_m %>% select(STATION_ID, fecha, VALUE),
            by = c("STATION_ID", "fecha"))

stopifnot(nrow(verificacion_M) == 36,
          all(!is.na(verificacion_M$VALUE)),
          all(verificacion_M$valor_original == verificacion_M$VALUE))

# Marcar los acumulados como NA en la copia diaria
p_m_depurado <- p_m_depurado %>%
  left_join(acumulados_tratar_M %>%
              transmute(STATION_ID, fecha, tratar_M = TRUE),
            by = c("STATION_ID", "fecha")) %>%
  mutate(VALUE = if_else(!is.na(tratar_M), NA_real_, as.numeric(VALUE))) %>%
  select(-tratar_M)

# Registrar cada modificacion
trazabilidad_M <- acumulados_tratar_M %>%
  mutate(valor_nuevo = NA_real_,
         tratamiento = "Excluido de serie diaria",
         justificacion = "Unico registro valido del mes, ubicado en el ultimo dia")

# Comprobar que no quedan casos sin tratar
acumulados_pendientes_M <- mensuales_M %>%
  anti_join(acumulados_tratar_M %>%
              mutate(YEAR = as.integer(format(fecha, "%Y")),
                     MONTH = as.integer(format(fecha, "%m"))) %>%
              select(STATION_ID, YEAR, MONTH),
            by = c("STATION_ID", "YEAR", "MONTH"))

# Cinco meses con ceros sospechosos y precipitacion concentrada al final
meses_ceros_tratar_M <- tibble(
  STATION_ID = c(2623300101, 2611400103, 2612710106,
                 2622900301, 2610000105),
  YEAR = c(1987, 1987, 1987, 1992, 1986),
  MONTH = c(11, 9, 7, 8, 12),
  motivo = "Ceros sospechosos y precipitacion concentrada al final"
)

# Siete meses completamente en cero con anomalia historica y regional
meses_todo_cero_tratar_M <- tibble(
  STATION_ID = c(2612000113, 2612000110, 2612000112,
                 2612000202, 2612700203, 2612720101,
                 2622600201),
  YEAR = c(1997, 2003, 2003, 1990, 1990, 1992, 1994),
  MONTH = c(1, 5, 5, 2, 9, 10, 9),
  motivo = "Mes completamente en cero con anomalia historica y regional"
)

# Unificar los doce meses seleccionados
meses_adicionales_M <- bind_rows(
  meses_ceros_tratar_M,
  meses_todo_cero_tratar_M
)

stopifnot(nrow(meses_adicionales_M) == 12,
          nrow(distinct(meses_adicionales_M, STATION_ID, YEAR, MONTH)) == 12)

# Guardar los valores originales antes de modificarlos
trazabilidad_adicional_M <- p_m_depurado %>%
  inner_join(meses_adicionales_M,
             by = c("STATION_ID", "YEAR", "MONTH")) %>%
  filter(!is.na(VALUE)) %>%
  transmute(STATION_ID, YEAR, MONTH, fecha,
            valor_original = VALUE,
            valor_nuevo = NA_real_,
            tratamiento = "Excluido de serie diaria",
            justificacion = motivo)

# Verificar que los doce meses tengan registros
stopifnot(nrow(distinct(trazabilidad_adicional_M,
                        STATION_ID, YEAR, MONTH)) == 12)

# Convertir todos los valores de los meses seleccionados en NA
p_m_depurado <- p_m_depurado %>%
  left_join(meses_adicionales_M %>%
              transmute(STATION_ID, YEAR, MONTH, tratar_adicional_M = TRUE),
            by = c("STATION_ID", "YEAR", "MONTH")) %>%
  mutate(VALUE = if_else(!is.na(tratar_adicional_M), NA_real_, as.numeric(VALUE))) %>%
  select(-tratar_adicional_M)

# Unificar la trazabilidad de todas las modificaciones
trazabilidad_total_M <- bind_rows(
  trazabilidad_M %>%
    mutate(YEAR = as.integer(format(fecha, "%Y")),
           MONTH = as.integer(format(fecha, "%m"))),
  trazabilidad_adicional_M
) %>%
  arrange(STATION_ID, fecha)

# Verificar que las modificaciones coincidan con la trazabilidad
stopifnot(
  nrow(trazabilidad_total_M) ==
    sum(!is.na(p_m$VALUE)) - sum(!is.na(p_m_depurado$VALUE))
)

# Comparacion antes y despues
contar_pares <- function(datos) {
  datos %>%
    filter(!is.na(VALUE)) %>%
    distinct(STATION_ID, YEAR) %>%
    nrow()
}

comparacion_M <- data.frame(etapa = c("Antes de M", "Despues de M"),
                            registros = c(nrow(p_m), nrow(p_m_depurado)),
                            validos = c(sum(!is.na(p_m$VALUE)),
                                        sum(!is.na(p_m_depurado$VALUE))),
                            estaciones = c(n_distinct(p_m$STATION_ID[!is.na(p_m$VALUE)]),
                                           n_distinct(p_m_depurado$STATION_ID[!is.na(p_m_depurado$VALUE)])),
                            pares_estacion_anio = c(contar_pares(p_m),
                                                    contar_pares(p_m_depurado)))

print(comparacion_M)

cat("Acumulados tratados:", nrow(trazabilidad_M), "\n")
cat("Meses adicionales tratados:", nrow(meses_adicionales_M), "\n")
cat("Registros adicionales convertidos en NA:", nrow(trazabilidad_adicional_M), "\n")
cat("Total de registros modificados en M:", nrow(trazabilidad_total_M), "\n")
cat("Ceros modificados:",
    sum(p_m$VALUE == 0, na.rm = TRUE) -
      sum(p_m_depurado$VALUE == 0, na.rm = TRUE), "\n")

# Tablas de trazabilidad de M
print(trazabilidad_total_M, n = Inf)
print(acumulados_pendientes_M, n = Inf)


####################################################
# Dias faltantes y años incompletos (1981-2003)
####################################################

# Construir calendario completo por estacion-anio
datos_N <- p_m_depurado %>%
  select(STATION_ID, YEAR, fecha, VALUE) %>%
  arrange(STATION_ID, fecha)

stopifnot(nrow(datos_N) == nrow(distinct(datos_N, STATION_ID, fecha)),
          !any(is.na(datos_N$fecha)))

fechas_N <- tibble(fecha = seq(as.Date("1981-01-01"), as.Date("2003-12-31"), by = "day"),
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
