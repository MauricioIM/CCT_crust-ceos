##Directorio
rm(list=ls()) # Borra los objetos creados

##Directorio
setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2021/stock_CAS/rfiles")
list.files()

source("file_1.R")

setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2021/stock_CAS/output")

#LEE EL REPORTE del modelo
dat1 <- reptoRlist("base.rep")
std1<-read.table("base.std",header=T,sep="",na="NA",fill=T)

setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2020/stock_CAS/output")

#LEE EL REPORTE del modelo
dat2 <- reptoRlist("base.rep")
std2<-read.table("base.std",header=T,sep="",na="NA",fill=T)

setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2019/stock_CAS/output")
## Cargar librerias ##

#LEE EL REPORTE del modelo
dat3 <- reptoRlist("base.rep")
std3<-read.table("base.std",header=T,sep="",na="NA",fill=T)


setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2018/stock_CAS/output")
## Cargar librerias ##

#LEE EL REPORTE del modelo
dat4 <- reptoRlist("base.rep")
std4<-read.table("base.std",header=T,sep="",na="NA",fill=T)


ano_2020 <- seq(1945,2020,1)
ano_2019 <- seq(1945,2019,1)
ano_2018 <- seq(1945,2018,1)
ano_2017 <- seq(1945,2017,1)

names(std1)

list(std1)

########BIOMASA DESOVANTE###########

BD1 <- std1[std1$name=="BD",,]$value
BDstd1 <- std1[std1$name=="BD",,]$std.dev
#BDi <- BD+BDstd*1.96
#BDs <- BD-BDstd*1.96
yy2 <- c((BD1-BDstd1*1.96), rev(BD1+BDstd1*1.96))

xx<-c(ano_2020,rev(ano_2020)) #polígonos


BD2 <- std2[std2$name=="BD",,]$value
BD3 <- std3[std3$name=="BD",,]$value
BD4 <- std4[std4$name=="BD",,]$value

length(ano_2020)
length(BD1)

png(file="C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2021/stock_CAS/figures/figure_6.png",
    width=1200, height=700)

#windows()
plot(ano_2020, BD1, type="n", ylab="Biomasa desovante (t)",
     xlab="Año", ylim=c(0,max(BD1)*1.7), main="Biomasa desovante", cex.lab=1,
     xaxt = 'n', cex=1)#, xaxp=c(1979,2015,25))
polygon(xx,yy2,col="darkgray",border="darkgray");lines(ano_2020,BD1,lwd=3, col="black")
axis(1, at = 1945:2020, labels = 1945:2020, cex=0.8)  

lines(ano_2020,BD1, pch=19, col="brown4", lwd=3)
lines(ano_2019,BD2, pch=19, col="blue", lwd=2)
lines(ano_2018,BD3, pch=19, col="orange", lwd=2)
lines(ano_2017,BD4, pch=19, col="green4", lwd=2)


legend("topright", c("Ago_2020 (CBA 2021)","Oct_2019 (CBA 2020)",
                        "Oct_2018 (CBA 2019)","Ago_2017 (CBA 2018)"), 
                    lty=c(1,1,1,1), pch=c(NA,NA,NA,NA), lwd=c(3,2,2,2), 
                    col=c("brown4","blue", "orange","green4"),
                    bty="n", cex=0.8)

dev.off()  


########BIOMASA TOTAL###########

BT1 <- std1[std1$name=="BT",,]$value
BTstd1 <- std1[std1$name=="BT",,]$std.dev
#BDi <- BD+BDstd*1.96
#BDs <- BD-BDstd*1.96
yy2 <- c((BT1-BTstd1*1.96), rev(BT1+BTstd1*1.96))

xx<-c(ano_2020,rev(ano_2020)) #polígonos


BT2 <- std2[std2$name=="BT",,]$value
BT3 <- std3[std3$name=="BT",,]$value
BT4 <- std4[std4$name=="BT",,]$value

length(ano_2020)
length(BT1)

png(file="C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2021/stock_CAS/figures/figure_7.png",
    width=1200, height=700)


#windows()
plot(ano_2020, BT1, type="n", ylab="Biomasa total (t)",
     xlab="Año", ylim=c(5000,max(BT1)*1.4), main="Biomasa total", cex.lab=1,
     xaxt = 'n', cex=1)#, xaxp=c(1979,2015,25))
polygon(xx,yy2,col="darkgray",border="darkgray");lines(ano_2020,BT1,lwd=3, col="black")
axis(1, at = 1945:2020, labels = 1945:2020, cex=0.8)  

lines(ano_2020,BT1, pch=19, col="brown4", lwd=3)
lines(ano_2019,BT2, pch=19, col="blue", lwd=2)
lines(ano_2018,BT3, pch=19, col="orange", lwd=2)
lines(ano_2017,BT4, pch=19, col="green4", lwd=2)


legend("bottomright", c("Ago_2020 (CBA 2021)","Oct_2019 (CBA 2020)",
                        "Oct_2018 (CBA 2019)","Ago_2017 (CBA 2018)"), 
       lty=c(1,1,1,1), pch=c(NA,NA,NA,NA), lwd=c(3,2,2,2), 
       col=c("brown4","blue", "orange","green4"),
       bty="n", cex=0.8)

dev.off()  


##########Mortalidad por pesca ########

F_h1 <- std1[std1$name=="log_Fh",,]$value
F_h1 <- exp(F_h1)
Fh_std1 <- std1[std1$name=="log_Fh",,]$std.dev

yy5 <- c((F_h1-Fh_std1*1.28), rev(F_h1+Fh_std1*1.28))

F_h2 <- std2[std2$name=="log_Fh",,]$value
F_h2 <- exp(F_h2)

F_h3 <- std3[std3$name=="log_Fh",,]$value
F_h3 <- exp(F_h3)

F_h4 <- std4[std4$name=="log_Fh",,]$value
F_h4 <- exp(F_h4)


png(file="C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2021/stock_CAS/figures/figure_8.png",
    width=1200, height=700)

#X11()

plot(ano_2020,F_h1,type="l",ylab="Mortalidad por pesca (F)",xlab="Año", 
     main="Mortalidad por pesca", ylim=c(0,3.5), cex.lab=1,
     xaxt = 'n', cex=1)#,las=1,xaxp=c(1979,2015,25))
polygon(xx,yy5,col="darkgray",border="darkgray");lines(ano_2020,F_h1,lwd=3, col="brown4")
axis(1, at = 1945:2020, labels = 1945:2020, cex=0.8) 

lines(ano_2020,F_h1, pch=19, col="brown4", lwd=3)
lines(ano_2019,F_h2, pch=19, col="blue", lwd=2)
lines(ano_2018,F_h3, pch=19, col="orange", lwd=2)
lines(ano_2017,F_h4, pch=19, col="green4", lwd=2)


legend("topleft", c("Ago_2020 (CBA 2021)","Oct_2019 (CBA 2020)",
                        "Oct_2018 (CBA 2019)","Ago_2017 (CBA 2018)"), 
       lty=c(1,1,1,1), pch=c(NA,NA,NA,NA), lwd=c(3,2,2,2), 
       col=c("brown4","blue", "orange","green4"),
       bty="n", cex=0.8)

dev.off()  


########RECLUTAS########


Rec1 <- std1[std1$name=="RecH",,]$value
Recstd1 <- std1[std1$name=="RecH",,]$std.dev
xx4<-c(ano_2020,rev(ano_2020))
yy4 <- c((Rec1-Recstd1 *1.96),rev(Rec1+Recstd1 *1.96))

Rec2 <- std2[std2$name=="RecH",,]$value
Rec3 <- std3[std3$name=="RecH",,]$value
Rec4 <- std4[std4$name=="RecH",,]$value


png(file="C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2021/stock_CAS/figures/figure_9.png",
    width=1200, height=700)


#x11()

plot(ano_2020,Rec1,type="l",ylab="Reclutas",xlab="Año", cex.lab=1.2,
     main="Reclutamiento", ylim=c(0,max(Rec1)*2.1), xaxt="n", cex.main=1.3)#,las=1,xaxp=c(1979,2015,25))
polygon(xx4,yy4,col="gray82",border="gray82");lines(ano_2020,Rec1,lwd=2)
axis(1, at = 1945:2020, labels = 1945:2020, cex=1.3)


lines(ano_2020,Rec1, pch=19, col="brown4", lwd=3)
lines(ano_2019,Rec2, pch=19, col="blue", lwd=2)
lines(ano_2018,Rec3, pch=19, col="orange", lwd=2)
lines(ano_2017,Rec4, pch=19, col="green4", lwd=2)


legend("topright", c("Ago_2020 (CBA 2021)","Oct_2019 (CBA 2020)",
                    "Oct_2018 (CBA 2019)","Ago_2017 (CBA 2018)"), 
       lty=c(1,1,1,1), pch=c(NA,NA,NA,NA), lwd=c(3,2,2,2), 
       col=c("brown4","blue", "orange","green4"),
       bty="n", cex=0.8)

dev.off() 

setwd("C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2021/stock_CAS/rfiles")
save.image(file = "file_7.RData")
