###############################
#Dornbusch Overshooting
#Guillermo Arroyo
#Facultad de Econom?a-UNAM
#Link: https://kevinkotze.github.io/ts-6-tut/
###############################

rm(list = ls())
graphics.off()

#Cargamos las librerias:
library(vars)
#library(ts)
library(urca)
library(car)
library(lmtest)
library(MTS)
library(aod)
library(graphics)
library(rugarch)
library(aTSA)
library(tidyverse)
library(lubridate)
library(timetk)
library(readxl)
library(dplyr)

#Cambiar el directorio de trabajo
setwd("D:/Doctorado Economía/Papers/Blockchain/Data")

#Cargar Datos: Bitcoin, Etherium, Chile, Mexico. 
#A<- read.csv("Mexico.csv", encoding="UTF-8-BOM")
A<-read_excel("Bitcoin.xlsx")

#Gr?fica:
plot_time_series(.data=A,.date_var=A$date, .value=A$open,.smooth_degree=2, .interactive = FALSE, .plotly_slider = TRUE, .title="US/Bitcoin", .y_lab = "US dollars")

#volatility
volatility <- sd(A$open, na.rm = FALSE)*sqrt(12/(length(A$open)-1))
volatility

#Basta cambiar la letra para cambiar el pa?s
TC <- A$open
date<- A$date
summary(TC)
#View(TC)

# Establecer que las variables son series de tiempo
TC <- ts(TC, frequency=12, start=c(2015,08))
lTC<-log(TC)
dlTC<-diff(lTC)

TC.desc = decompose(TC)
plot(TC.desc, xlab='A?o')

TC.desc$random
TC.desc$random <-na.omit( TC.desc$random)

plot(TC, plot.type = "single")
legend("bottomright", legend = c("A"), lty = 1, 
       col = c("blue"), bty = "n")

#Funci?n de autocorrelaci?n:
#acf(TC, max.lag = 20,plot =FALSE)
#plot(acf)
#pacf(TC, max.lag = 20,plot =FALSE)
#plot(pacf)

#Link: https://kevinkotze.github.io/ts-6-tut/
#H0: la prueba no es estacionaria, tiene una ra?z unitaria
#HA: la prueba es estacionario, no tiene ra?z unitaria
#Si |t*|<= |tau 95%|, entonces, Acepte H0, rechace a HA. Serie no estacionaria.
#Si |t*|> |tau 95%|, entonces, Acepte a HA, rechace H0. serie estacionaria.
#Test Dickey-Fuller:
dfc.ger <- ur.df(lTC, lags=12, type = "drift", selectlags = c("AIC"))
summary(dfc.ger)
plot(dfc.ger)

dfc.ger <- ur.df(TC.desc$random, lags=12, type = "none", selectlags = c("AIC"))
summary(dfc.ger)
plot(dfc.ger)



#===================Modelos EGARCH=======================
#Asimetr?a POSITIVA 
#?C?mo afecta al TC las BUENAS NOTICIAS?
#Se tiene que asegurar que la serie sea estacionaria, por ello tomamos la primera diferencia. 
#gamma 1 pondera la asimetr?a cuando estamos en la parte positiva
#Ho: Si el par?metro es estad?sticamente igual a cero: NO hay asimetr?a positiva
#H1: Si el par?metro NO es estad?sticamente igual a cero: hay asimetr?a positiva
#Criterios: M?x varosimilitud y m?n AIC

spec <- ugarchspec(mean.model=list(armaOrder=c(2,2)), 
                   variance.model= list(model="eGARCH",garchOrder=c(2,4)),
                   distribution.model = "sstd")
mod_egarch <- ugarchfit(dlTC, spec = spec)
print(mod_egarch)
coef(mod_egarch)
#Par?metro de sesgo gamma1 (en rugarch: skew): cuando gamma1 = 1: simetr?a. Cuando gamma1 <1: sesgo negativo. Para gamma1> 1: sesgo positivo.
likelihood(mod_egarch)
infocriteria(mod_egarch)
stdmsftret <- residuals(mod_egarch, standardize = TRUE)
Box.test(abs(stdmsftret), 22, type = "Ljung-Box")
#Regla de oro: el valor de p inferior al 5% indica que el modelo utilizado no es v?lido.


#===================Modelos TGARCH=======================
#Asimetr?a NEGATIVA
#?C?mo afecta al TC las MALAS NOTICIAS?
#eta11 es el par?metro de la dummmy
#Ho: Si eta11 es estad?sticamente igual a cero: NO hay asimetr?as negativas
#H1: Si eta11 NO es estad?sticamente igual a cero: s? hay asimetr?as negativas
#Criterios: M?x varosimilitud y m?n AIC
spec <- ugarchspec(mean.model=list(armaOrder=c(2,2)), 
                   variance.model= list(model="fGARCH",garchOrder=c(2,4),
                                        submodel = "TGARCH"),
                   distribution.model = "sstd")
mod_tgarch <- ugarchfit(dlTC, spec = spec)
print(mod_tgarch)
coef(mod_tgarch)
#Par?metro de sesgo ?? (en rugarch: skew): cuando ?? = 1: simetr?a. Cuando ?? <1: sesgo negativo. Para ??> 1: sesgo positivo.
likelihood(mod_tgarch)
infocriteria(mod_tgarch)
stdmsftret <- residuals(mod_tgarch, standardize = TRUE)
Box.test(abs(stdmsftret), 22, type = "Ljung-Box")
#Regla de oro: el valor de p inferior al 5% indica que el modelo utilizado no es v?lido.


#====================  Modelos ARCH-M =======================
#Estimacion del modelo GARCH(1,1)
#Mide el efecto del riesgo sobre el comportamiento promedio del rendimiento
#arch-m es la varianza que se estim? con un proceso garch y luego lo 
#incorpor? como variable ex?gena en el modelo de rendimiento
#Ho: Si arch-m es estad?sticamente igual a cero: NO hay efecto del riesgo sobre el rendimiento
#H1: Si arch-m NO es estad?sticamente igual a cero: s? hay efecto del riesgo sobre el rendimiento
#Criterios: M?x varosimilitud y m?n AIC
spec <- ugarchspec(mean.model=list(armaOrder=c(1,1),archm = TRUE), 
                   variance.model= list(model="sGARCH",garchOrder=c(1,1)),
                   distribution.model = "sstd")
mod_arch_m <- ugarchfit(dlTC, spec = spec)
print(mod_arch_m)
coef(mod_arch_m)
#Par?metro de sesgo ?? (en rugarch: skew): cuando ?? = 1: simetr?a. Cuando ?? <1: sesgo negativo. Para ??> 1: sesgo positivo.
likelihood(mod_arch_m)
infocriteria(mod_arch_m)
stdmsftret <- residuals(mod_arch_m, standardize = TRUE)
Box.test(abs(stdmsftret), 22, type = "Ljung-Box")
#Regla de oro: el valor de p inferior al 5% indica que el modelo utilizado no es v?lido.






















