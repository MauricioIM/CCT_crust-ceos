#################################################################################
##                                                                            ##  
##                       SALIDAS MODELO  ANUAL EN TALLAS                      ##
##                              PARA CAMARÓN NAILON                           ##
##                           MODIFICADO DE ELSON LEAL                         ##
##                                     IFOP                                   ##
################################################################################

## Cargar librerias ##
rm(list=ls())           # Limpieza de variables del ambiente, útil cada vez que empezamos un nuevo análisis

library(gplots)
library(ggplot2)
library(reshape)
#library(FLEDA)
library(mgcv)
library(lattice)
library(nlme)
library(plotrix)


## cargan salidas reporte modelo ##

setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2018/stock_CAN/rfiles")
getwd()

source("file_1.R")

setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2018/stock_CAN/output")

S1<-reptoRlist("base.rep") # salidas base 2016
names(S1)


std<-read.table("base.std",header=T,sep="",na="NA",fill=T)

table(std$name)
names(std)


###### Seteo variables extra a utilizar en los plots############

#---------------------------------------------#
ano <- seq(1961,2017,1)
#---------------------------------------------#
anos <- seq(1961,2017,1)
years <- seq(1961,2017,1)

############# Reemplazo valores 0  por NA, para graficar ########## 

#class(R0prom)
length(anos)
names(S1)

########## BIOMASAS ########## 

#----------------#
#Cargar librerias#
#----------------#
library(reshape)
#library(FLEDA)
library(mgcv)
library(lattice)
#trellis.device(color=F)
library(nlme)

#-------------------------------#
# cargan salidas reporte modelo #
#-------------------------------#

#--------------------------#
# Elimino celdas sin datos #
#--------------------------#

table(std$name)
names(std)

R0 <- exp(std[std$name=="log_Ro",,]$value)

class(R0)
length(anos)

#-------------------------------------------------#
# Se abre el archivo STD arrojado por el ADMB #

B <- read.table("base.std",header=T,sep="",na="NA",fill=T)


#--------------------------------------#
# Elementos para Gráfico BIOMASA Total #
#--------------------------------------#

Bs    <- B[B[,2]=="BT",,]; Bs1 <- matrix(Bs$value,ncol=1) # Biomasa anual
Bs_Li <- matrix((Bs$value-1.96*Bs$std),ncol=1)                               # Límite inferior Biomasa anual
Bs_Ls <- matrix((Bs$value+1.96*Bs$std),ncol=1)                               # Límite superior Biomasa anual
xx1 <- c(years,rev(years))
yy1 <- c(Bs_Li,rev(Bs_Ls))

##------------------------------------##
## Elementos para Grafico BIOMASA DESOVANTE          ##
##------------------------------------##
Bd    <- B[B[,2]=="BD",,]; Bd1  <- matrix(Bd$value,ncol=1)                   # Biomasa desovante
Bd_Li <- matrix((Bd$value-1.96*Bd$std),ncol=1)                               # Límite inferior Biomasa desovante
Bd_Ls <- matrix((Bd$value+1.96*Bd$std),ncol=1)                               # Límite superior Biomasa desovante
xx2    <- c(years,rev(years))
yy2    <- c(Bd_Li,rev(Bd_Ls))

##------------------------------------##
## Elementos para Grafico Reclutamiento          ##
##------------------------------------##
Rh    <- B[B[,2]=="RecH",,]; Recl  <- matrix(Rh$value,ncol=1)                   # Biomasa desovante
Rh_Li <- matrix((Rh$value-Rh$std),ncol=1)                               # Límite inferior Biomasa desovante
Rh_Ls <- matrix((Rh$value+Rh$std),ncol=1)                               # Límite superior Biomasa desovante
xx4    <- c(years,rev(years))
yy4    <- c(Rh_Li,rev(Rh_Ls))

##------------------------------------ ##
## Elementos para Grafico Mortalidad por pesca          ##
##------------------------------------##
Fh <- std[std$name=="log_Fh",,]$value
Fh <- exp(Fh)
Fh_std <- std[std$name=="log_Fh",,]$std.dev
xx<-c(anos,rev(anos)) #polígonos
yy5 <- c((Fh-Fh_std*1.28), rev(Fh+Fh_std*1.28))

################################################################
#   GRAFICOS DESVIACION ESTANDAR DE VARIABLES POBLACIONALES    #
################################################################

png(file="C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2018/stock_CAN/figures/figure_1.png",
    width=1200, height=700)

#####BIOMASA TOTAL#######

par(mfrow=c(2,2),mar=c(4,8,1,0.5)+1)
plot   (xx1, yy1, type="n", ylim=c(5000,max(yy1)),
        ylab="Biomasa (t)",las=1,xlab="Año",cex.lab=1.4, cex.axis=1, main = "Biomasa total (t)")
polygon(xx1, yy1, col="grey", border = "grey");lines  (years,Bs1)


#####BIOMASA DESOVANTE#######

plot(xx2, yy2, type="n", ylim=c(0,max(yy2)),
     ylab="Biomasa (ton)",las=1,xlab="Año",cex.lab=1.4, main = "Biomasa desovante (t)")
polygon(xx2, yy2, col="gray", border = "gray");lines  (years,Bd1)


#####RECLUTAMIENTO HEMBRAS#######

plot(xx4, yy4, type="n", ylim=c(0,max(yy4)),
     ylab="Reclutas (millones de hembras)",las=1,xlab="Año",cex.lab=1.4, main = "Reclutamiento")
polygon(xx4, yy4, col="gray", border = "gray");lines  (years,Recl)
abline(h=R0,col="black",lwd=1.5)


#####MORTALIDAD POR PESCA HEMBRAS#######

plot(xx, yy5, type="n", ylim=c(0,max(yy5)),
     ylab="Mortalidad por pesca (F)",las=1,xlab="Año",cex.lab=1.4, main = "Mortalidad por pesca")
polygon(xx, yy5, col="gray", border = "gray");lines  (years,Fh)
#abline(h=R0,col="black",lwd=1.5)

dev.off()

setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2018/stock_CAN/rfiles")

save.image(file = "file_2.RData")

#names(S1)

