rm(list=ls())           # Limpieza de variables del ambiente, �til cada vez que empezamos un nuevo an�lisis

library(gplots)
library(ggplot2)
library(reshape)
#library(FLEDA)
library(mgcv)
library(lattice)
library(nlme)
library(plotrix)

setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crust�ceos/2020/stock_CAS/rfiles")
getwd()

source("file_1.R")

setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crust�ceos/2020/stock_CAS/output")

S1<-reptoRlist("base.rep") # salidas base 2016
names(S1)

std<-read.table("base.std",header=T,sep="",na="NA",fill=T)

table(std$name)
names(std)


#####DENSIDAD DE BDo###########


#define population mean and standard deviation

population_mean <- std[std$name=="SSBo",,]$value
population_sd <- std[std$name=="SSBo",,]$std.dev

#define upper and lower bound
lower_bound <- population_mean - population_sd
upper_bound <- population_mean + population_sd

#Create a sequence of 1000 x values based on population mean and standard deviation
x <- seq(-4, 4, length = 1000) * population_sd + population_mean

#create a vector of values that shows the height of the probability distribution
#for each value in x
y <- dnorm(x, population_mean, population_sd)

BDo = 0.4*9618.5

#plot normal distribution with customized x-axis labels

png(file="C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crust�ceos/2020/stock_CAS/figures/figure_4.png",
    width=1200, height=700)

#windows()
plot(x,y, type = "l", lwd = 2, axes = F, xlim=c(3000, 18000), xlab = "BDo", ylab = "Probabilidad",
     main = "Distribuci�n de probabilidad de BDo")
abline(v=BDo, col="red")

sd_axis_bounds = 5
axis_bounds <- seq(-sd_axis_bounds * population_sd + population_mean,
                   sd_axis_bounds * population_sd + population_mean,
                   by = population_sd)

axis(side = 1, at = seq(3000, 18000, by=500), pos = 0, xlab="Probabilidad")
axis(side = 2, at = seq(0, 1, by=0.0001), pos = 3000)

legend("topright", c("Densdidad", "40%BDo"), lty=c(1,1), 
       pch=c(NA,NA), col=c("black", "red"), 
       bty="n", cex=1.3, lwd=c(2,2))


dev.off()


setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crust�ceos/2020/stock_CAS/rfiles")

save.image(file = "file_5.RData")


