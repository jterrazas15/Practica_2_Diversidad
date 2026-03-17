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

### Indice de Shannon ###

Diversidad <- Diversidad_fichas %>% column_to_rownames("Sitio") 

alfa <- data.frame( 
  sitio    = rownames(Diversidad), 
  riqueza  = specnumber(Diversidad), 
  shannon  = diversity(Diversidad, index = "shannon"), 
  simpson  = diversity(Diversidad, index = "simpson"), 
  inv_simp = diversity(Diversidad, index = "invsimpson"), 
  pielou   = diversity(Diversidad, "shannon") / log(specnumber(sp_matrix)) 
)
