####################################################
# Proyecto 3 
# Por: Natalia Paredes, Mariana Cifuentes, John Franco
# y Estefanía Hernández
####################################################
##### LIBRERIAS
#install.packages("fda")
#install.packages("fda.usc")
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
setwd("~/Downloads/Proyecto Final")
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






