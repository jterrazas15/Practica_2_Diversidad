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
#### Matriz de especie-sitio ####
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
    
# B) Curvas de rarefacción

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
### Grafico rango-abundancia ###
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

# C. Gráficas de rango-abundancia
h) P8. Identifica la forma general de cada curva (log-normal, geométrica, etc.). ¿Qué modelo
ecológico subyace a cada forma? ¿Qué mecanismos comunitarios podrían explicarla?
El el caso del sitio A solo pudismos observar un punto en la parte superior de la grafica, esto nos indica que tiene una sola especie. 
El el caso del sitio B,C,D,E y G presentan forma geometrica, el modelo ecologico subyacente es el de particion secuencial de recursos, donde la especie dominante ocupa la mayor fraccion del recurso disponible, esto ocurre en comunidades simples con alta dominancia.
En el sitio F no se observa ningun punto o linea en la grafica. 

i) P9. ¿En qué sitio la especie más dominante representa la mayor fracción de individuos?
¿Cómo se relaciona ese resultado con el índice de Simpson?
En la grafica podemos observar que la especie mas dominante representa la mayor fraccion de individuos en el sitio B, ya que su primera especie tiene el valor mas alto de log10 cercano a 0, que en la escala original equivale a una abundancia relativa de 1 o sea casi el 100% de los individuos.
Se relaciona con el indice de Simpson ya que mide la probabilidad de que dos individuos seleccionados al azar sean de la misma especie (dominancia), en este caso el sitio B tendra un valor alto en el indice de Simpson ya que su primera esoecie como lo mencione representa mas de la mitad de los individuos, pero tendra un valor bajo en el indice se Simpson inverso, ya qie representa el numero efectivo de especies equivalentes. 

j) P10. Imagina que en un monitoreo de largo plazo una curva log-normal se convierte en
geométrica. ¿Qué perturbación ecológica podría explicarlo?
El cambio de este nos indica una fuerte perturbacion que simplifica el ecosistema, algunos ejemplos de perturbaciones que podrian ser la causante de esto serian: 
-Perturbacion antropogenica severa, como la deforestacion, contaminacion extrema o indencios forestaes, los cuales eliminan a la mayoria de espedies y permiten la proliferacion de unas pocas. 
-Introduccion de una nueva especie invasora muy competitiva, como alguna especie exotica puede ser hiperdominante lo cual pudiese desplazar a las demas y pueden acaparar los recursos. 
-Estres ambiental cronico, Cambio en las condiciones ambientales como lo es el aumento de temperatura, que supera los limites de tolerancia de la mayoria de las especies, dejando solo a las mas resientes. 



##################################
### Diversidad beta ###
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

# D. Diversidad beta
k) P11. Analiza la matriz de Bray-Curtis. ¿Qué par de sitios es más diferente? ¿Es el mismo
par más diferente según Jaccard? Si no coinciden, ¿a qué se debe?
En este caso la matriz de Bray-Curtis, el par mas distinro es el sitio_A y sitio_D los cuales tienen un valor de disimilitud de 0.98, en cambio en el Jaccard el par de sitios mas diferente es el sitio_C y sitio_A con un valor de 0.92. 
No es el mismo par, esto se debe a que jaccard es un indice cualitativo es decir esta basado solo en ausencia o presencia de especies, el alto valor entre C y A nos indica que comparten muy pocas especies, en cambio Bray-Curtis es un indice cuantitativo ya que considera las abundancias de las especies, el alto valor entre A y D nos indica que comparten pocas especies, las abundancias que comparten son diferentes. 

l) P12. ¿Qué componentes del recambio beta (reemplazo vs. anidamiento) podrían explicar
los patrones que observas? (Pista: investiga el paquete betapart en R.)
El componente de reemplazo es muy fuerte entre ciertos pares de sitios, como lo es en el caso del sitio A (tiene una sola especie) vs el sitio C muestran un alto reemplazo, ya que sus especies son completamente diferentes, esto nos indica un cambio en la identidad de las especies a lo largo de un gradiente ambiental. 
El componente de anidamiento: Los sitios con forma geometrica y menor riqueza podrian ser subconjuntos anidados de sitios mas diversos como el B. Esto sugeriria un proceso de filtrado ambiental o perdida de especies, para confirmarlo tendriamos que usar el paquete betapart el cual descompone la disimilitud total (Jaccard) en: 
-βJTU=reemplazo (cambio de especies)
-βJNE=anidamiento (perdida de especie)

m) P13. Si estos sitios representaran puntos a lo largo de un gradiente altitudinal, ¿qué patrón
de beta-diversidad esperarías según la teoría del dominio medio (mid-domain effect)?
¿Coincide con tus resultados?
El patron esperado predice que la riqueza de especies sera maxima en el centro del gradiente y disminuira hacia los extremos, lo cual nos generaria un patron de anidamiento, donde las comunidades de los extremos son subconjuntos depauperados de la comunidad del centro, con poco reemplazo de especies. Nuestros resultados no coinciden completamente, ya que aunque podemos observar un cierto anidamiento (C,D,E,G podrian ser subconjuntos de B), tambien hay un fuerte componente de reemplazo, especialmente entre el sitio A y los demas, lo que nos sugiere que ademas de la perdida de especies hay in cambio de identidad de especies a lo  largo del gradiente. 
