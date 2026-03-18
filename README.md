# Practica_2_Diversidad
Reporte y análisis de diversidad.
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

##################################
##### Curvas de rarefaccion #####
################################

lista_sit <- split(t(matriz), rownames(t(matriz)))
salir <- iNEXT(lista_sit, q = 0, datatype = "abundance")

grafico_rare <- ggiNEXT(salir, type = 1) +
  scale_color_viridis_d(option = "plasma") +
  labs( title = "Curvas de rarefaccion", 
        x = "Numero de individuos muestreados", 
        y= "Riqueza de especies") +
  theme_classic(base_size = 12)

grafico_rare


ggsave("Figuras/rarefaccion.jpg",
       grafico_rare,
       width = 8,
       height = 6)
    
B) Curvas de rarefacción

e) P5. ¿Cuál sitio alcanza la asíntota más rápidamente? ¿Qué te dice eso sobre la estructura
de la comunidad?
En la grafica observamos que varios sitios alcanza una estabilizacion relativa rapida en la riqueza de especies con un numero pequeño de individuos muestreadis. Los sitios representados por curvas que se aplanan cerca de los valores de riqueza entre aproximadamente 4 y 5 especies alcanzan la asintota con menor esfuerzo, esto nos puede sugerir que tienen una riqueza relativamente baja o una estructura dominada por pocas especies, en cambio algunas curvas continuan aumentando gradualmente con el incremento del numero de los individuos y esto nos indica que aun podria detectarse nuevas especies. 

f) P6. Dos sitios pueden tener la misma riqueza observada pero curvas con pendientes muy
distintas. ¿Cómo interpretas esa diferencia en términos de abundancias relativas?
Cuando ddos sitios podrian tener una riqueza similar al final de la curva, pero pendientes distintas, cuando la pendiente inicial es muy pronunciada esto indica que en kas primeras muestras se observaron mas especies dominantes, en cambio en una pendiente mas gradual las especies tienen abundancias mas equilibradas, esto nos muestra la probabilidad de encontrar nuevas especies a medida que aumenta el tamaño. De manera ecologica esto lo podriamos ver como que las diferencias nos reflejan variaciones en la estructura de abundancias relativas de las comunidades, donde algunas estan dominadas por pocas espcies y otras presentan mayor equitatividad. 

g) P7. ¿A qué tamaño de muestra estandarizado compararías los sitios? Justifica tu elección y
calcula la riqueza rarefactada para todos los sitios a ese nivel.
Como podemos observar en la grafica la mayoria de las curvas terminan alrededor de 30 a 40 individuos, por lo que un tamaño de muestta cercano a este rango seria apropiado para comparar los sitios de manera equitativa, al estandarizar el tamaño de muestta en ese nivel se puede estimar la riqueza rarefactada esperada para sitio bajo el mismo muestreo, esto nos permite realizar comparaciones mas justas entre comunidades, eliminando el efecto del tamaño de la muestra sobre la estimacion de la diversidad, para lo cua hicimos el siguiente codigo: 

###################################
### Riqueza rarefactada (P7) ###
###################################

# calculamos el número total de individuos por sitio
n_individuos <- rowSums(matriz) #Calculamos el numero total de individuos registrados en cada sitio utilizando la funcion rowSums
#  mínimo de muestra
n_estandar <- min(n_individuos) #Aqui seleccionamos el menor numero de individuos entre todos los sitios, utilizandolo como referencia para la rarefaccion, ya que nos ayuda a comparar la riqueza de especies.
n_estandar #Observamos el tamaño de la muestra utilizado para estandarizar la comparacion entre sitios.

# calculamos la riqueza rarefactada
riqueza_rarefactada <- rarefy(matriz, sample = n_estandar) #Utilizamos la funcion rarefy para calcular la riqueza rarefactada, donde nos estima el numero esperado de especies que se observa en cada sitio si tosos hubieran sido muestreados.

# sacamos la tabla 
tabla_rarefaccion <- data.frame(
  Sitio = rownames(matriz),
  Riqueza_rarefactada = riqueza_rarefactada
) #Aqui creamos un data.frame que contenga el nombre de cada sitio y su correspondiente valor de riqueza rarefactada

rownames(tabla_rarefaccion) <- NULL #Esto lo puse pq me salen repetidos los nombres y esto es para que los elimine 
# Imprimimos la tabla 
tabla_rarefaccion
#Imprimimos la tabla para observar los valores de cada uno.


##################################
#####Grafico rango-abundancia####
################################

rank_da <- Diversidadd_fichas %>%
  group_by(Sitio, Muestra) %>%
  summarise(n = sum(Abundancias)) %>%
  filter(n > 0) %>%
  group_by(Sitio) %>%
  mutate( prop = n/ sum(n),
          rango = rank(-prop, ties.method = "first")) %>%
  ungroup()
grafico_rank <- ggplot(rank_da, aes(x=rango, y = log10(prop), color = Sitio)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  labs( title = "Grafica de rango- abundancia", 
        x= "Rango de especie",
        y = "log10(Abundancia relativa)") +
  theme_classic()

grafico_rank

ggsave ("Figuras/rank_abundance.jpg",
        grafico_rank,
        width = 8,
        height = 6)
##################################
##########Diversidad beta########
################################

library(reshape2)

jac <- vegdist(matriz, method = "jaccard", binary = T)

brac <- vegdist( matriz, method = "bray")

plot_be <- function (dist_obj, titulo) {
  melt (as.matrix(dist_obj)) %>%
    ggplot(aes(Var1, Var2, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(value,2)), size = 3.5) +
    scale_fill_gradient(
      low = "#E8F5E9",
      high = "#1B5E20",
      name = "Disimilitud"
    ) +
    labs ( title = titulo, 
           x = NULL,
           y = NULL) +
    theme_minimal(base_size = 11)+
    theme( axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

plot_be(jac, "Disimilitud de jaccard")
plot_be(brac, "Disimilitud de Bray-Curtis")
