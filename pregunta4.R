

# Leer el archivo corregido
datos <- read_excel("datos_cafe_7rep_final.xlsx")


# Seleccionar solo variables numéricas desde ese mismo objeto
X <-  datos[c(2,3,4,5,6,7)]

# Verificación
cat("nrow(X):", nrow(X), "\n")              # debe dar 49
cat("length(Hibrido):", length(datos$Hibrido), "\n")  # debe dar 49

# Ahora sí correr Box's M
resultado <- boxM(X, datos$Hibrido)
print(resultado)

# Variables respuesta:
# rendimiento, altura_planta, diametro_fruto, largo_fruto, peso_grano, contenido_cafeina
# Factor: Hibrido

# =====================================================
# SUPUESTOS
# =====================================================
library(readxl)
datos <- read_excel("datos_cafe_7rep.xlsx")
View(datos_cafe_7rep)
# Normalidad multivariada por grupo
library(tidyverse)
library(mvnormtest)
summary(datos)
library(tidyverse)
trat1 = datos %>% filter(Hibrido == "Hibrido1") %>%
  dplyr::select(rendimiento, altura_planta, diametro_fruto, largo_fruto, peso_grano, contenido_cafeina)
trat2 = datos %>% filter(Hibrido == "Hibrido2") %>%
  dplyr::select(rendimiento, altura_planta, diametro_fruto, largo_fruto, peso_grano, contenido_cafeina)
trat3 = datos %>% filter(Hibrido == "Hibrido3") %>%
  dplyr::select(rendimiento, altura_planta, diametro_fruto, largo_fruto, peso_grano, contenido_cafeina)
trat4 = datos %>% filter(Hibrido == "Hibrido4") %>%
  dplyr::select(rendimiento, altura_planta, diametro_fruto, largo_fruto, peso_grano, contenido_cafeina)
trat5 = datos %>% filter(Hibrido == "Hibrido5") %>%
  dplyr::select(rendimiento, altura_planta, diametro_fruto, largo_fruto, peso_grano, contenido_cafeina)
trat6 = datos %>% filter(Hibrido == "Hibrido6") %>%
  dplyr::select(rendimiento, altura_planta, diametro_fruto, largo_fruto, peso_grano, contenido_cafeina)
trat7 = datos %>% filter(Hibrido == "Hibrido7") %>%
  dplyr::select(rendimiento, altura_planta, diametro_fruto, largo_fruto, peso_grano, contenido_cafeina)

#Con cada hibrido no hay normalidad entonces deberian plantearse 
#transformaciones o usar las pruebas
#robustas que se describen más abajo.

# Supuesto de homogeneidad de matrices variancia covariancia

library(heplots)
datos$Hibrido <- as.factor(datos$Hibrido)

res <- boxM(X, datos$Hibrido)
res
summary(res)

heplots::boxM(cbind(rendimiento, altura_planta, diametro_fruto, largo_fruto, peso_grano, contenido_cafeina) ~ Hibrido, data = datos)

library(biotools)

resultado <- boxM(X, datos$Hibrido)
print(resultado)

library(covTestR)

cafe <- unique(datos$Hibrido)
cafe1 <- lapply(cafe,
                function(x){as.matrix(datos[datos$Hibrido == x, 2:7])}
)

names(cafe1) <- cafe
Ahmad2017(cafe1)

## Prueba Wrapper
homogeneityCovariances(datos, group = Hibrido, covTest = BoxesM)

library(DFA.CANCOR)
HOMOGENEITY(data = datos,groups = 'Hibrido', 
            variables = c('rendimiento','altura_planta','diametro_fruto','largo_fruto','peso_grano',"contenido_cafeina"))

# Supuesto de variables dependientes correlacionadas. 
# Prueba de esfericidad de Bartlett

library(psych)
options(scipen = 0)
cortest.bartlett(cor(datos[, -1]), n = nrow(datos[, -1]))

library(MVTests)
res1 = Bsper(datos[, -1])
summary(res1)

#Trabajando con el modelo de MANOVA en DCA
modelo = manova(cbind(rendimiento, altura_planta, diametro_fruto, largo_fruto, peso_grano, contenido_cafeina) ~ Hibrido, data = datos)

#Determinación de la matriz residual y la matriz factorial del MANOVA.
str(modelo)
Matrices = summary(modelo)$SS
F = Matrices$Hibrido
W = Matrices$Residuals

#Variabilidad explicada por el factor (Tratamientos). Matriz suma de 
#cuadrados y productos cruzados del factor (SCOCF)
F

#Variabilidad residual. Matriz suma de cuadrados y productos cruzados
#del residual (SCOCR)
W

#Variabilidad Total. Matriz suma de cuadrados y productos cruzados total (SCOCT)
#del factor 
T = F + W
T

#Bondad de ajuste. Un valor proximo a 1 indica que la mayor parte
#de la variabilidad total puede atribuirse al factor, mientras que un 
#valor proximo a 0 significa que el factor explica muy poco
#de esa variabilidad total.

eta2 = 1 - det(W)/det(T)
eta2
det(F)/det(T) # La razon de determinantes no es aditiva como las
#suma de cuadrados.

# Pruebas de hipotesis del modelo. Calculo de contrastes del modelo

#Contrastes del modelo en relacion a los supuestos. Todos los 
#estadisticos son bastante robustos ante violaciones de normalidad,
#y la prueba de Roy es muy sensible a violaciones de la hipotesis
#de la matriz de covariancias. Cuando las muestras son iguales por
#grupo, la prueba de Pillai es el estadistico más robusto ante
#violaciones de los supuestos.

k = 7 #numero de grupos
p = 6 #numero de variables
n = 7 #numero de observaciones por grupo
datosc1 = datos
head(datosc1)
datosc1$Hibrido <- as.numeric(datosc1$Hibrido)
str(datosc1)

VMPG = matrix(NA, k, p) #vector de medias por grupo
for(i in 1:k){
  VMPG[i,]=colMeans(datos[datosc1$Hibrido == i, -1])
}
VMPG #cada fila es un vector de medias

#Computar B
(B=n*(k-1)*cov(VMPG))
# Tambien
n*(t(VMPG)-colMeans(VMPG))%*%t(t(VMPG)-colMeans(VMPG))

# Computar W
W = (n-1)*cov(datos[datosc1$Hibrido == 1, -1])
for(i in 2:k){
  W = W + (n-1)*cov(datos[datosc1$Hibrido == 1, -1])
}
W
W+B

## autovalores de W^{-1}B
(lambdas = eigen(solve(W)%*%B)$values)

##traza de pillai
sum(lambdas/(1+lambdas))
#tambien
sum(diag(solve(W+B)%*%B))

#El software R obtiene un valor calculado de pillai diferente al que se calculó 
#manualmente, esto se debe a que hay varias formas de calcular este valor
summary(modelo, test = "Pillai")
1-pf(2.2491, 36,252)

##Lambda de Wilks
det(W)/det(W + B)
#tambien
prod(1/(1 + lambdas))

#El software R obtiene un valor calculado diferente al que se calculó 
#manualmente, esto se debe a que hay varias formas de calcular este valor
summary(modelo, test = "Wilks")
1-pf(2.5341, 36,165.24)

## Lawley Hotelling
LH = sum(lambdas)
LH

summary(modelo, test = "Hotelling-Lawley")
1-pf(2.681, 36,212)

## Raiz mayor de Roy
lambdas[1]/(1 + lambdas[1])
summary(modelo, test = "Roy")
1-pf( 9.0125, 6, 42)

# Se observa que las variables rendimiento y altura de planta resultaron altamente significativas, lo que indica que contribuyen de manera importante al rechazo de la hipótesis nula que planteaba que los vectores de medias de los híbridos eran iguales. En consecuencia, se concluye que al menos uno de los vectores de medias difiere entre híbridos. Por otro lado, las variables diámetro de fruto y contenido de cafeína mostraron una tendencia a la significancia, mientras que largo de fruto y peso de grano no evidenciaron diferencias estadísticamente significativas entre los híbridos evaluados.

summary.aov(modelo)

#Comparar dos hibridos

modelo1 = manova(cbind(rendimiento, altura_planta, diametro_fruto, largo_fruto, peso_grano, contenido_cafeina) ~ Hibrido, data = datos,
                 subset = Hibrido %in% c("Hibrido1","Hibrido2"))

summary(modelo1,test="Pillai")
summary(modelo1,test="Wilks")
summary(modelo1,test="Hotelling-Lawley")
summary(modelo1,test="Roy")
summary.aov(modelo1)

#---------------------------------------------------#
#           Comparación General por pares           #
#---------------------------------------------------#

# H0:los 2 vectores de promedios son iguales
# H1:los 2 vectores de promedios difieren

Prueba <- c("Hibrido1", "Hibrido2", "Hibrido3", "Hibrido4", "Hibrido5", "Hibrido6","Hibrido7")
comb<-t(combn(length(Prueb) , 2))

for(i in 1:nrow(comb)){
  modelo.comp = manova(cbind(rendimiento, altura_planta, diametro_fruto, largo_fruto, peso_grano, contenido_cafeina) ~ Hibrido, data=datos,
                       subset = Hibrido %in% Prueba[comb[i,]])
  print(paste("Prueba: ", Prueb[comb[i,]][1], "y", Prueba[comb[i,]][2]))
  print(summary(modelo.comp, test = "Pillai"))
  cat("\n")
  
}
#Se concluye que no todos los híbridos son iguales entre sí, ya que varias comparaciones por pares resultaron significativas. En particular, el Híbrido4 presenta diferencias claras con la mayoría de los otros híbridos, mientras que otros (como Híbrido6 y Híbrido7) no mostraron diferencias notorias respecto a la mayoría.