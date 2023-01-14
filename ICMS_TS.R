#########################################################
#########################################################
##               Time Series Analysis                  ##
#########################################################
#########################################################

### Packages used
require(fpp2)
require(xts)
library(ggplot2)
library(tseries)
library(forecast)
library(astsa)

### Load data ###
ipeadata <- read.csv("C:/Users/RodrigoSalles/Desktop/ICMS_TS/data_3")
View(ipeadata)

### Delete first column ###
ipeadata$X <- NULL
View(ipeadata) 

### Initial Graphic ###
ipeadata.ts = ts(ipeadata$valor, frequency = 12, start = c(1993,1))
plot(ipeadata.ts, col=1, lwd=2, main = 'ICMS - Minas Gerais(01/1993 - 12/2015 ', xlab = 'Years', ylab = 'Values')
grid (NULL,NULL, lty = 6,  col = "cornsilk2") 

### Initial Exploratory Analysis ###
summary(ipeadata)
sd(ipeadata$valor)


### Stationarity ###
# ACF
acf(ipeadata$valor, main = 'Autocorrelation function')

### Decomposition of the Time Series ###
ipeadata_dec =decompose(ipeadata.ts)
plot(ipeadata_dec)

### Augmented Dickey Fuller Test : raiz unitária ###
adf.test(ipeadata$valor) #p-value = 0.5561 - Non-stationary series

### KPSS Test ###
kpss.test(ipeadata$valor) # p-value = 0.01 - Non-stationary series

### Number of Differentiations ###
ndiffs(ipeadata$valor) # Only one differentiation is needed

### Number of Seasonal Differentiations ###
nsdiffs(ipeadata.ts) # No need for seasonal differentiations



### Treatment to Achieve Stationarity ###
# Applying one Differentiation
ipea_dif_1 = diff(ipeadata$valor, differences = 1)

# Checking ACF
acf(ipea_dif_1, main = 'Autocorrelation function')
pacf(ipea_dif_1, main = 'Partial autocorrelation function')

ndiffs(ipea_dif_1) # No need for new differentiations

# Augmented Dickey Fuller Test : unit root , 1 dif
adf.test(ipea_dif_1) # p-value = 0.01, with one differentiation the series is stationary.

# KPSS Test
kpss.test(ipea_dif_1) # p-value = 0.1, with one differentiation the series is stationary.


#### Models ####
# initial Model
mod = auto.arima(ipeadata$valor)
mod # ARIMA(1,1,1) with drift

# Arima(1,1,1) model with drift
mod = Arima(ipeadata.ts, order = c(1,1,1), include.drift = T)

# Check model parameters and residuals of arima(1,1,1)
sarima(ipeadata$valor, 1,1,1)

### Testing the models ###
sarima(ipeadata$valor, 2,1,1 )
sarima(ipeadata$valor, 1,1,2 )
sarima(ipeadata$valor, 3,1,1 )
sarima(ipeadata$valor, 3,1,2)
sarima(ipeadata$valor, 3,1,3 )
sarima(ipeadata$valor, 1,1,2 )
sarima(ipeadata$valor, 1,1,3 )


# Model with the best result
mod1 <- Arima(ipeadata.ts,order = c(3,1,3),include.drift=TRUE)

# Check model parameters and residuals of arima(3,1,3)
sarima(ipeadata$valor, 3,1,3 )


# Testing the model with seasonality 
mod_saz = sarima(ipeadata$valor, 1,1,1, 2,0,1,12)
mod_saz = Arima(ipeadata$valor,order=c(1,1,1), seasonal=c(2,0,1), include.drift = T)
mod_saz$ttable

# Testing the models with Box Cox transformation
# Lambda ideal value
BoxCox.lambda(ipeadata.ts) # 0.5993851
# Checking the efect of the transformatiom
plot.ts(BoxCox(ipeadata.ts, lambda = 0.5993851))
# Residue analysis 
acf(box_trans$residuals)



#### Predictions ####
# Forecast for the next 12 months: mod
forecast(mod, h = 12, level = c(80,95))
mod_prev = forecast(mod, h = 12, level = c(80,95))
# Plotting the graphs with predictions mod
plot(forecast(mod, h = 12, level = c(80,95)),xlab = 'Years', ylab = 'Values')

# Forecast for the next 12 months: mod1
forecast(mod1, h = 12, level = c(80,95))
# Plotting the graphs with predictions mod1
plot(forecast(mod1, h = 12, level = c(80,95)),xlab = 'Years', ylab = 'Values')

# Forecast for the next 12 months: mod_saz
forecast(mod_saz, h = 12, level = c(80,95))
# Plotting the graphs with predictions mod_saz
plot(forecast(mod_saz, h = 12, level = c(80,95)),xlab = 'Years', ylab = 'Values')

# Forecast for the next 12 months: box_trans
box_trans = forecast(auto.arima(ipeadata.ts, lambda=0.5993851, biasadj=TRUE),h=12)
forecast(box_trans, h = 12, level = c(80,95))
# Plotting the graphs with predictions box_trans
plot(forecast(box_trans, h = 12, level = c(80,95)),xlab = 'Years', ylab = 'Values')



###### Artificial Neural Networks - ANN ######
fit <- nnetar(ipeadata.ts, lambda=0)
autoplot(forecast(fit,h=12), xlab = 'Anos', ylab = 'Values', main = 'Forecast - RNA - 2016')
fit

fcast <- forecast(fit, PI=TRUE, h=12)
autoplot(fcast, xlab = 'Years', ylab = 'Values', main = 'Forecast - ANN - 2016')

# Residue analysis 
acf(fcast$residuals[which(!is.na(fcast$residuals))], main = 'Autocorrelation function - ANN')



#### Combination of Models ####

train <- ipeadata.ts

ETS <- forecast(ets(train), h=12)
ARIMA <- forecast(auto.arima(train, lambda=0, biasadj=TRUE),h=12)
STL <- stlf(train, lambda=0, h=12, biasadj=TRUE)
NNAR <- forecast(nnetar(train), h=12)
TBATS <- forecast(tbats(train, biasadj=TRUE), h=12)
Combination <- (ETS[["mean"]] + ARIMA[["mean"]] + NNAR[["mean"]] + TBATS[["mean"]])/4
Combination


# Residue analysis 
par(mfrow=c(2,2))

acf(ETS$residuals)

acf(ARIMA$residuals)

acf(STL$residuals)

acf(TBATS$residuals)

acf(na.remove(fit$residuals))


