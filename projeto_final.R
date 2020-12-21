
#-- Removendo variáveis carregadas anteriormente
rm(list = ls())
graphics.off()

#--- Instalação e carregamento das bibliotecas
#install.packages("XLConnect")
#install.packages("rJava")
#system("java -version")
#Sys.setenv(JAVA_HOME="/usr/lib/jvm/java-8-openjdk-amd64/bin/")

library("DBI")
library("rJava")
library(zoo)
library(XLConnect)
library(lubridate)
library(ggplot2)

#--- Pré-prossesamento
#----- Dados do IBCBR
#Importando os dados e renomeando as colunas
setwd("~/Documentos/USP/2_sem/SeriesTemps/trabalho_final/")
wb <- loadWorkbook("IBCBR.xls")
Matriz <- readWorksheet(wb, sheet = "Séries", startRow = 0, startCol = 0)

#alterando o vetor de datas para o formato correto   
names(Matriz) <- c("Data", "IBCBr")
Matriz$Data <- parse_date_time(Matriz$Data, "%Y.%m")

#transformando dados em time-series     
tsIBCBr <- zoo(Matriz$IBCBr, Matriz$Data)
inf <- ts(Matriz$IBCBr, start = c(2003, 1), frequency = 12)

#gerando gráficos mais elaborados      
lc <- qplot(Data, IBCBr, data = Matriz, geom = "line") + 
  geom_line(colour = "darkgrey", size = .7)
lc + ggtitle("IBC-Br (2002=100, Banco Central do Brasil)") + 
  theme(plot.title = element_text(lineheight=.9, face="bold"))

#------ Modelo linear dinâmico

#devtools::install_github("KevinKotze/tsm")
#install.packages("dlm", repos = "https://cran.rstudio.com/", 
#                 dependencies = TRUE)
library(dlm)
library(tsm)


### Estimação de parâmetros do modelo
fn <- function(parm) {
  dlmModPoly(order = 1, dV = exp(parm[1]), dW = exp(parm[2]))
}
fit <- dlmMLE(inf, rep(0, 2), build = fn, hessian = TRUE)
(conv <- fit$convergence) # teste para ver se convergio
  
loglik <- dlmLL(inf, dlmModPoly(1))
n.coef <- 2
r.aic <- (2 * (loglik)) + 2 * (sum(n.coef))  #dlmLL caculates the neg. LL
r.bic <- (2 * (loglik)) + (log(length(inf))) * (n.coef)

# Criação do modelo com parâmetros estimados
mod <- fn(fit$par)
obs.error.var <- V(mod)
state.error.var <- W(mod)

# Kalman filter
filtered <- dlmFilter(inf, mod = mod)
smoothed <- dlmSmooth(filtered)

resids <- residuals(filtered, sd = FALSE)
mu <- dropFirst(smoothed$s)
mu.1 <- mu[1]
mu.end <- mu[length(mu)]

#-- prev
conf.tmp <- unlist(dlmSvd2var(smoothed$U.S, smoothed$D.S))
#length(conf.tmp)
conf <- ts(conf.tmp[-1], start = c(2003, 1), frequency = 12)
#length(conf)
wid <- qnorm(0.05, lower = FALSE) * sqrt(conf)
#length(wid)
conf.pos <- mu + wid
conf.neg <- mu - wid

comb.state <- cbind(mu, conf.pos, conf.neg)

forecast <- dlmForecast(filtered, nAhead = 6) # prevendo os 3 ultimos meses de 2020
var.2 <- unlist(forecast$Q)
wid.2 <- qnorm(0.05, lower = FALSE) * sqrt(var.2)
comb.fore <- cbind(forecast$f, forecast$f + wid.2, forecast$f - 
                     wid.2)

result <- ts(rbind(comb.state, comb.fore), start = c(2003, 1), frequency = 12)

##-- Critérios de informaçaõ e Testes de hipóteses
cat("AIC", r.aic)
cat("BIC", r.bic)
cat("V.variance", obs.error.var)
cat("W.variance", state.error.var)

Box.test(resids, lag = 12, type = "Ljung", fitdf = 2)
shapiro.test(resids)  # normality


#----- plot
#par(mfrow = c(2, 1), mar = c(2.2, 2.2, 1, 1), cex = 0.8) # dois plots na fig
par(mfrow = c(2, 1), mar = c(2.2, 2.2, 1, 1), cex = 0.8) # um plot em cada fig
plot.ts(inf, col = "darkgrey", xlab = "", ylab = "", lwd = 1.5)
lines(mu, col = "black")
legend("bottomright", legend = c("IBC-Br", "DLM polinomial de ordem 1"),
       lwd = c(2, 1), col = c("darkgrey", "black"), bty = "n")

plot.ts(resids, ylab = "", xlab = "", col = "darkgrey", 
        lwd = 1.5)
abline(h = 0)
legend("topleft", legend = "Resíduos", lwd = 1.5, col = "darkgrey", 
       bty = "n")

#--- Histograma
par(mfrow = c(1, 1), mar = c(2.2, 2.2, 1, 1), cex = 0.8)
hist(resids, prob = TRUE, col = "grey", 
     main = "Residuos", breaks = seq(-4.5, 4, length.out = 30))

#-- grafico de previsão
par(mfrow = c(1, 1), mar = c(2.2, 2.2, 1, 1), cex = 0.8)
plot.zoo(result, 
         plot.type = "single", 
         col = c("black", "red", "red"),
         xlab = "", ylab = "")
lines(inf, col = "darkgrey", lwd = 1.5)
abline(v = 2020.7 , col = "blue", lwd = 1, lty = 3)
legend("topleft", 
       legend = c("IBC-Br", "DLM Polinomial de ordem 1", "Intervalo de confiança"), 
       lwd = c(1.5, 1, 1), col = c("darkgrey", "black", "red"), bty = "n")


#-- ACF e PACF
ac(resids)

#--- Modelo hibrido
seasonSize = 4
fn <- function(parm) {
  mod <- dlmModPoly(order = 1) + dlmModSeas(frequency = seasonSize)
  V(mod) <- exp(parm[1])
  diag(W(mod))[1:2] <- exp(parm[2:3])
  return(mod)
}
# Encontrando par?metros do modelo
fit <- dlmMLE(inf, rep(0, seasonSize-1), build = fn, hessian = TRUE)
(conv <- fit$convergence)  

loglik <- dlmLL(inf, dlmModPoly(1) + dlmModSeas(seasonSize))
n.coef <- 3
r.aic <- (2 * (loglik)) + 2 * (sum(n.coef))  #dlmLL caculates the neg. LL
r.bic <- (2 * (loglik)) + (log(length(inf))) * (n.coef)

mod <- fn(fit$par) #dlmHibridModel <- buildFunHib(fit$par)
obs.error.var <- V(mod) #drop(V(dlmHibridModel))
state.error.var <- diag(W(mod)) #diag(W(dlmHibridModel))

# Suavização
filtered <- dlmFilter(inf, mod = mod) 
smoothed <- dlmSmooth(filtered) #hybSmooth <- dlmSmooth(IBCBr_ts, mod = dlmHibridModel)
resids <- residuals(filtered, sd = FALSE)

mu <- dropFirst(smoothed$s[, 1])
gammas <- dropFirst(smoothed$s[,2])
trend <-  dropFirst(smoothed$s[,3])
mu.1 <- mu[1]
mu.end <- mu[length(mu)]
gammas.1 <- gammas[1]
gammas.end <- gammas[length(mu)]


#--- plot
alpha <- mu + gammas
par(mfrow = c(2, 1), mar = c(2.2, 2.2, 1, 1), cex = 0.8)
plot.ts(inf, col = "darkgrey",main='', xlab = "", ylab = "", lwd = 2)
lines(alpha, col = "black")
legend("bottomright", legend = c("IBC-Br", "Modelo Híbrido"), 
       lwd = c(2, 1), col = c("darkgrey", "black"), bty = "n")
plot.ts(resids, ylab = "", xlab = "", col = "darkgrey", 
        lwd = 2)
abline(h = 0)
legend("bottomleft", legend = "Residuos", lwd = 2, col = "darkgrey", 
       bty = "n")

#--- Histograma
par(mfrow = c(1, 1), mar = c(2.2, 2.2, 1, 1), cex = 0.8)
hist(resids, prob = TRUE, col = "grey", 
     main = "Residuos", breaks = seq(-4.5, 4, length.out = 30))


ac(resids)  # acf

#--- estatisticas e testes de hipoteses
cat("AIC", r.aic)
cat("BIC", r.bic)
cat("V.variance", obs.error.var)
cat("W.variance", state.error.var)
Box.test(resids, lag = 12, type = "Ljung", fitdf = 2)  # joint autocorrelation
shapiro.test(resids)  # normality
