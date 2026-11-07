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

# Borde del Valle
departamentos_col <- gadm(country = "COL", level = 1, path = tempdir())  
valle <- departamentos_col[departamentos_col$NAME_1 == "Valle del Cauca", ]

ggplot() + geom_spatraster(data = p_satelite[[1]]) +
  geom_sf(data = valle, fill = NA, color = "black", linewidth = 0.5) +
  labs(title = "Precipitación 01-01-1981") 














