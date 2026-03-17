###################################
### Carga y validacion de datos ###
###################################

install.packages("tidyverse")
install.packages("iNEXT")

library(tidyverse)
library(vegan)
library(iNEXT)

glimpse(Diversidad_fichas)

stopifnot(!any(is.na(Diversidad_fichas)))

###################################
### Indices de diversidad alpha ###
###################################

#########################
### Indice de Shannon ###
#########################

Shannon <- function(abun){
  total <- sum(abun)
  frec <- abun/total
  mult <- (log(frec))*frec
  s<- -(sum(mult))
  return(s)
}


### Tabla con Sitio, Abundancia total y Shannon ###

tabla_shannon <- Diversidad_fichas %>% # Se crea un objeto, donde se guarara la tabla, se asigna el data frame, y %>% toma la base de datos y aplica el siguiente argumento
  group_by(Sitio) %>% # aqui se agrupa por el sitio. y se le aplica los filtros que siguen
  summarise( # crea el resumen 
    Abundancia_total = sum(Abundancias), # suma las abundanciass de cadda muestreo del sitio
    shannon = Shannon(Abundancias) # Aplica el indice de shannon con la funcion creada en clase 
  )

print(tabla_shannon) # se visualiza la tabla

#########################
### Indice de Simpson ###
#########################

simpsonSimplificado <- function(ab){
  tot <- sum(ab)
  cuadrados <- (ab/tot)^2
  adecion <- sum(cuadrados)
  return (adecion)
}

tabla_simpson <- Diversidad_fichas %>% # Se crea un objeto, donde se guarara la tabla, se asigna el data frame, y %>% toma la base de datos y aplica el siguiente argumento
  group_by(Sitio) %>% # aqui se agrupa por el sitio. y se le aplica los filtros que siguen
  summarise( # crea el resumen 
    Abundancia_total = sum(Abundancias), # suma las abundanciass de cadda muestreo del sitio
    simpson = simpsonSimplificado(Abundancias) # Aplica el indice de shannon con la funcion creada en clase 
  )

tabla_simpson

#################################
### Indice de Simpson inverso ###
#################################

inverSimpson <- function (abun){
  return(1/simpsonSimplificado(abun))
}

tabla_insimpson <- Diversidad_fichas %>% # Se crea un objeto, donde se guarara la tabla, se asigna el data frame, y %>% toma la base de datos y aplica el siguiente argumento
  group_by(Sitio) %>% # aqui se agrupa por el sitio. y se le aplica los filtros que siguen
  summarise( # crea el resumen 
    Abundancia_total = sum(Abundancias), # suma las abundanciass de cadda muestreo del sitio
    Simpson_inverso = inverSimpson(Abundancias) # Aplica el indice de shannon con la funcion creada en clase 
  )

tabla_insimpson

########################
### Indice de Pielou ###
########################

pielou <- function(abundancias){
  S <- sum(abundancias > 0)
  J <- Shannon(abundancias)/log(S)
  return (J)
}

Diversidad_fichas$Muestra[Diversidad_fichas$Muestra == "morado"] <- 1
Diversidad_fichas$Muestra[Diversidad_fichas$Muestra == "azul_turquesa"] <- 2
Diversidad_fichas$Muestra[Diversidad_fichas$Muestra == "azul_cielo"] <- 3
Diversidad_fichas$Muestra[Diversidad_fichas$Muestra == "azul_marino"] <- 4
Diversidad_fichas$Muestra[Diversidad_fichas$Muestra == "fucsia"] <- 5
Diversidad_fichas$Muestra[Diversidad_fichas$Muestra == "rosa_medio"] <- 6
Diversidad_fichas$Muestra[Diversidad_fichas$Muestra == "rosa_claro"] <- 7
Diversidad_fichas$Muestra[Diversidad_fichas$Muestra == "rojo"] <- 8
Diversidad_fichas$Muestra[Diversidad_fichas$Muestra == "verde_limon"] <- 9
Diversidad_fichas$Muestra[Diversidad_fichas$Muestra == "verde"] <- 10
Diversidad_fichas$Muestra[Diversidad_fichas$Muestra == "naranja"] <- 11
Diversidad_fichas$Muestra[Diversidad_fichas$Muestra == "amarillo"] <- 12

Diversidad_fichas$Muestra <- as.numeric(Diversidad_fichas$Muestra)

tabla_pielou <- Diversidad_fichas %>% # Se crea un objeto, donde se guarara la tabla, se asigna el data frame, y %>% toma la base de datos y aplica el siguiente argumento
  group_by(Sitio) %>% # aqui se agrupa por el sitio. y se le aplica los filtros que siguen
  summarise( # crea el resumen 
    Abundancia_total = sum(Abundancias), # suma las abundanciass de cadda muestreo del sitio
    Pielou = pielou(Abundancias) # Aplica el indice de shannon con la funcion creada en clase 
  )

tabla_pielou

indices_alpha <- tabla_pielou %>%
  left_join(tabla_insimpson, by = "Sitio") %>%
  left_join(tabla_shannon, by = "Sitio") %>%
  left_join(tabla_simpson, by = "Sitio")

indices_alpha

library(ggplot2)
install.packages("gridExtra")
library(gridExtra)

## Crear la tabla ##
tabla <- tableGrob(indices_alpha) # Transforma los objetos que forman la tabla en un plot (imagen)

## Crear un gráfico vacío y añadir la tabla ##
indices <- ggplot() +
  annotation_custom(tabla)
## Guardar el plot en resultados ##

# Supongamos que tu objeto gráfico se llama indices
ggsave("Figuras/idices_alpha.jpg", plot = indices, width = 12.5, height = 2.5)

##################################
####### Riqueza de especies ######
#################################

riqueza <- Diversidadd_fichas %>%
  group_by(Sitio) %>%
  summarise( riqueza = sum(Abundancias > 0))

riqueza


##################################
####### Matriz de especie-sitio ######
#################################

matriz <- Diversidadd_fichas %>%
  pivot_wider(names_from = Muestra,
              values_from =  Abundancias,
              values_fill = 0)
matriz <- as.data.frame(matriz)
rownames(matriz) <- matriz$Sitio 
matriz <- matriz [, -1]

matriz <- as.matrix(matriz)

chao1 <- estimateR (matriz)
tabla_chao1 <- data.frame(
  Sitio = rownames(matriz),
  Chaoo1 = chao1 ["S.chao1",]
)
rownames(tabla_chao1) <- NULL
tabla_chao1




