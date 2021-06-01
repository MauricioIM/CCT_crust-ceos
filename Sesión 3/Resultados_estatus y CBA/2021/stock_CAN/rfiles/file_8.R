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

population_mean_2020 <- exp(std[std$name=="log_Ro",,]$value)
population_sd_2020 <- exp(std[std$name=="log_Ro",,]$std.dev)

#define upper and lower bound 2020
lower_bound_2020 <- 500
upper_bound_2020 <- 800


#define population mean and standard deviation 2019

population_mean_2019 <- exp(std2[std2$name=="log_Ro",,]$value)
population_sd_2019 <- exp(std2[std2$name=="log_Ro",,]$std.dev)

#define upper and lower bound 2019

lower_bound_2019 <- population_mean_2019 - population_sd_2019
upper_bound_2019 <- population_mean_2019 + population_sd_2019


#define population mean and standard deviation 2018

population_mean_2018 <- exp(std3[std3$name=="log_Ro",,]$value)
population_sd_2018 <- exp(std3[std3$name=="log_Ro",,]$std.dev)

#define upper and lower bound 2018

lower_bound_2018 <- population_mean_2018 - population_sd_2018
upper_bound_2018 <- population_mean_2018 + population_sd_2018

#define population mean and standard deviation 2017

population_mean_2017 <- exp(std4[std4$name=="log_Ro",,]$value)
population_sd_2017 <- exp(std4[std4$name=="log_Ro",,]$std.dev)

#define upper and lower bound 2017

lower_bound_2017 <- population_mean_2017 - population_sd_2017
upper_bound_2017 <- population_mean_2017 + population_sd_2017


#########Plot distribuciones#########


png(file="C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2021/stock_CAN/figures/figure_10.png",
    width=1200, height=700)


plot.function(x = function(t) dnorm(x = t, mean = population_mean_2020, sd = population_sd_2020),
              from = 570,
              to = 770,
              col = "brown4",
              xlab = "Ro",
              ylab = "Densidad",
              main = "Densidades de Ro")
plot.function(x = function(t) dnorm(x = t, mean = population_mean_2019, sd = population_sd_2019),
              from = 720,
              to = 750,
              col = "blue",
              add = TRUE)

plot.function(x = function(t) dnorm(x = t, mean = population_mean_2018, sd = population_sd_2018),
              from = 720,
              to = 750,
              col = "orange",
              add = TRUE)

plot.function(x = function(t) dnorm(x = t, mean = population_mean_2017, sd = population_sd_2017),
              from = 600,
              to = 650,
              col = "green4",
              add = TRUE)

legend(x = "top",
       legend = paste("CBA_", 2021:2020:2019:2018), cex = 1.3, 
       col = c("brown4", "blue", "orange", "green4"),
       lty = 1)



dev.off()


setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2021/stock_CAN/rfiles")

save.image(file = "file_8.RData")



