rm(list=ls())           # Limpieza de variables del ambiente, útil cada vez que empezamos un nuevo análisis

library(gplots)
library(ggplot2)
library(reshape)
#library(FLEDA)
library(mgcv)
library(lattice)
library(nlme)
library(plotrix)

setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2019/stock_CAN/rfiles")
getwd()

source("file_1.R")

setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2019/stock_CAN/output")

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

BDo = 0.4*5901.8

#plot normal distribution with customized x-axis labels

png(file="C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2019/stock_CAN/figures/figure_3.png",
    width=1200, height=700)

plot(x,y, type = "l", lwd = 2, axes = F, xlab = "BDo", ylab = "Probabilidad",
     main = "Distribución de probabilidad de BDo")
#abline(v=BDo, col="red")

sd_axis_bounds = 5
axis_bounds <- seq(-sd_axis_bounds * population_sd + population_mean,
                   sd_axis_bounds * population_sd + population_mean,
                   by = population_sd)

axis(side = 1, at = axis_bounds, pos = 0, xlab="Probabilidad")
axis(side = 2, at = seq(0, 1, by=0.001), pos = 4879.04)

dev.off()


setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2019/stock_CAN/rfiles")

save.image(file = "file_4.RData")


