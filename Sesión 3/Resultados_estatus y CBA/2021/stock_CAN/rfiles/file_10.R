rm(list=ls())           # Limpieza de variables del ambiente, útil cada vez que empezamos un nuevo análisis

library(gplots)
library(ggplot2)
library(reshape)
#library(FLEDA)
library(mgcv)
library(lattice)
library(nlme)
library(plotrix)

setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2021/stock_CAN/rfiles")
getwd()

source("file_1.R")


setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2021/stock_CAN/output")

S1<-reptoRlist("base.rep") # salidas base 2016
names(S1)

std<-read.table("base.std",header=T,sep="",na="NA",fill=T)

setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2020/stock_CAN/output")

S2<-reptoRlist("base.rep") # salidas base 2016
names(S2)

std2<-read.table("base.std",header=T,sep="",na="NA",fill=T)


setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2019/stock_CAN/output")

S3<-reptoRlist("base.rep") # salidas base 2016
names(S3)
std3<-read.table("base.std",header=T,sep="",na="NA",fill=T)

setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2018/stock_CAN/output")

S4<-reptoRlist("base.rep") # salidas base 2016
names(S4)
std4<-read.table("base.std",header=T,sep="",na="NA",fill=T)


table(std$name)
names(std)


#####DENSIDAD DE R0###########


#define population mean and standard deviation 2020


population_mean_2020 <- exp(std[199,"value"])
population_sd_2020 <- std[199,"std.dev"]

class(std)

#define population mean and standard deviation 2019

population_mean_2019 <-  exp(std2[196,"value"])
population_sd_2019 <- std2[196,"std.dev"]


#define population mean and standard deviation 2018

population_mean_2018 <-  exp(std3[193,"value"])
population_sd_2018 <- std3[193,"std.dev"]

#define population mean and standard deviation 2017

population_mean_2017 <- exp(std4[190,"value"])
population_sd_2017 <- std4[190,"std.dev"]



#########Plot distribuciones#########


png(file="C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2021/stock_CAN/figures/figure_12.png",
    width=1200, height=700)

#windows()


plot.function(x = function(t) dnorm(x = t, mean = population_mean_2018, sd = population_sd_2018),
              from = -0.5,
              to = 0.8,
              col = "orange",
              xlab = "Mortalidad por pesca (F)",
              ylab = "Densidad",
              main = "Densidades de F")
plot.function(x = function(t) dnorm(x = t, mean = population_mean_2019, sd = population_sd_2019),
              from = -0.5,
              to = 0.8,
              col = "blue",
              add = TRUE)

plot.function(x = function(t) dnorm(x = t, mean = population_mean_2020, sd = population_sd_2020),
              from =-0.5,
              to = 0.8,
              col = "brown4",
              add = TRUE)

plot.function(x = function(t) dnorm(x = t, mean = population_mean_2017, sd = population_sd_2017),
              from =-0.5,
              to = 0.8,
              col = "green4",
              add = TRUE)

legend(x = "topright",
       legend = paste("CBA_", 2021:2020:2019:2018), cex = 1.3, 
       col = c("brown4", "blue", "orange", "green4"),bty="n",
       lty = 1, lwd = 2)



dev.off()


setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2021/stock_CAN/rfiles")

save.image(file = "file_10.RData")



