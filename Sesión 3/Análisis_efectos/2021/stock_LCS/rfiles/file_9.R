rm(list=ls())           # Limpieza de variables del ambiente, útil cada vez que empezamos un nuevo análisis

library(gplots)
library(ggplot2)
library(reshape)
#library(FLEDA)
library(mgcv)
library(lattice)
library(nlme)
library(plotrix)

setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2021/stock_LCS/rfiles")
getwd()

source("file_1.R")


setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2021/stock_LCS/output")

S1<-reptoRlist("base.rep") # salidas base 2016
names(S1)

std<-read.table("base.std",header=T,sep="",na="NA",fill=T)

setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2020/stock_LCS/output")

S2<-reptoRlist("base.rep") # salidas base 2016
names(S2)

std2<-read.table("base.std",header=T,sep="",na="NA",fill=T)


setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2019/stock_LCS/output")

S3<-reptoRlist("base.rep") # salidas base 2016
names(S3)
std3<-read.table("base.std",header=T,sep="",na="NA",fill=T)

setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2018/stock_LCS/output")

S4<-reptoRlist("base.rep") # salidas base 2016
names(S4)
std4<-read.table("base.std",header=T,sep="",na="NA",fill=T)


table(std$name)
names(std)



#####DENSIDAD DE BDo###########


#define population mean and standard deviation 2020

population_mean_2020 <- std[std$name=="SDo2",,]$value
population_sd_2020 <- std[std$name=="SDo2",,]$std.dev


#define population mean and standard deviation 2019

population_mean_2019 <- std2[std2$name=="SDo2",,]$value
population_sd_2019 <- std2[std2$name=="SDo2",,]$std.dev


#define population mean and standard deviation 2018

population_mean_2018 <- std3[std3$name=="SDo2",,]$value
population_sd_2018 <- std3[std3$name=="SDo2",,]$std.dev

#define population mean and standard deviation 2017

population_mean_2017 <- std4[std4$name=="SDo2",,]$value
population_sd_2017 <- std4[std4$name=="SDo2",,]$std.dev



#########Plot distribuciones#########


png(file="C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2021/stock_LCS/figures/figure_11.png",
    width=1200, height=700)



#windows()

plot.function(x = function(t) dnorm(x = t, mean = population_mean_2017, sd = population_sd_2017),
              from = 20000,
              to = 80000,
              col = "green4",
              xlab = "Biomasa desovante virgen (t)",
              ylab = "Densidad",
              main = "Densidades de BDo")
plot.function(x = function(t) dnorm(x = t, mean = population_mean_2019, sd = population_sd_2019),
              from = 20000,
              to = 80000,
              col = "blue",
              add = TRUE)

plot.function(x = function(t) dnorm(x = t, mean = population_mean_2018, sd = population_sd_2018),
              from = 20000,
              to = 80000,
              col = "orange",
              add = TRUE)

plot.function(x = function(t) dnorm(x = t, mean = population_mean_2020, sd = population_sd_2020),
              from = 20000,
              to = 80000,
              col = "brown4",
              add = TRUE)

legend(x = "topleft",
       legend = paste("CBA_", 2021:2020:2019:2018), cex = 1.3, 
       col = c("brown4", "blue", "orange", "green4"), bty="n",
       lty = 1, lwd = 2)



dev.off()


setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2021/stock_LCS/rfiles")

save.image(file = "file_9.RData")



