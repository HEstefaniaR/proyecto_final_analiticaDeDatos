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

# Casos con al menos 3 meses de este tipo
alta_sospecha_M <- patrones_mensuales %>%
  filter(meses_con_un_registro >= 3)

cat("Pares estacion-anio con patron repetido:", nrow(alta_sospecha_M), "\n")

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
            motivo = "Registro unico a fin de mes con patron repetido",
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


# Limpieza final y trazabilidad
p_m_depurado <- p_m

# Los 26 registros de patron mensual repetido
acumulados_tratar_M <- auditoria_acumulados_M %>%
  select(STATION_ID, fecha, valor_original, motivo) %>%
  distinct(STATION_ID, fecha, .keep_all = TRUE)

# Verificar que los registros coincidan con la base
stopifnot(nrow(acumulados_tratar_M) == 26)

verificacion_M <- acumulados_tratar_M %>%
  left_join(p_m %>% select(STATION_ID, fecha, VALUE),
            by = c("STATION_ID", "fecha"))

stopifnot(nrow(verificacion_M) == 26,
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
         justificacion = paste("Unico registro al cierre del mes",
                               "con patron repetido en el mismo anio"))

# Los 10 casos aislados se guardan sin modificar
acumulados_pendientes_M <- mensuales_M %>%
  anti_join(acumulados_tratar_M %>%
              mutate(YEAR = as.integer(format(fecha, "%Y")),
                     MONTH = as.integer(format(fecha, "%m"))) %>%
              select(STATION_ID, YEAR, MONTH),
            by = c("STATION_ID", "YEAR", "MONTH"))

# Verificar que los 313 ceros se conservaron
ceros_antes <- sum(p_m$VALUE == 0, na.rm = TRUE)
ceros_despues <- sum(p_m_depurado$VALUE == 0, na.rm = TRUE)

stopifnot(ceros_antes == ceros_despues)

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
cat("Casos aislados conservados:", nrow(acumulados_pendientes_M), "\n")
cat("Ceros modificados:", ceros_antes - ceros_despues, "\n")

# Tablas de trazabilidad
print(trazabilidad_M, n = Inf)
print(acumulados_pendientes_M, n = Inf)
