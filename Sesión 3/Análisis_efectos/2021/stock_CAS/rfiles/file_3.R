rm(list=ls())           # Limpieza de variables del ambiente, �til cada vez que empezamos un nuevo an�lisis

library(gplots)
library(ggplot2)
library(reshape)
#library(FLEDA)
library(mgcv)
library(lattice)
library(nlme)
library(plotrix)

setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crust�ceos/2021/stock_CAS/rfiles")
getwd()

source("file_1.R")

setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crust�ceos/2021/stock_CAS/output")

S1<-reptoRlist("base.rep") # salidas base 2016
names(S1)

std<-read.table("base.std",header=T,sep="",na="NA",fill=T)

table(std$name)
names(std)


#####DENSIDAD DE R0###########


#define population mean and standard deviation

population_mean <- exp(std[std$name=="log_Ro",,]$value)
population_sd <- exp(std[std$name=="log_Ro",,]$std.dev)

#define upper and lower bound
lower_bound <- population_mean - population_sd
upper_bound <- population_mean + population_sd

#Create a sequence of 1000 x values based on population mean and standard deviation
x <- seq(-4, 4, length = 1000) * population_sd + population_mean

#create a vector of values that shows the height of the probability distribution
#for each value in x
y <- dnorm(x, population_mean, population_sd)

#plot normal distribution with customized x-axis labels

png(file="C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crust�ceos/2021/stock_CAS/figures/figure_2.png",
    width=1200, height=700)


#windows()
plot(x,y, type = "l", lwd = 2, axes = F, xlab = "Ro", ylab = "Probabilidad",
     main = "Distribuci�n de probabilidad de Ro")
sd_axis_bounds = 5
axis_bounds <- seq(-sd_axis_bounds * population_sd + population_mean,
                   sd_axis_bounds * population_sd + population_mean,
                   by = population_sd)
axis(side = 1, at = axis_bounds, pos = 0, xlab="Probabilidad")
axis(side = 2, at = seq(0, 1, by=0.1), pos = 1138.481)

dev.off()


setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crust�ceos/2021/stock_CAS/rfiles")

save.image(file = "file_3.RData")



