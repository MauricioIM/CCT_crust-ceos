rm(list=ls())


df1 <- data.frame(Año=seq(2018,2021,1), min=c(803,594,246,1112), 
                 max=c(955,712,314,1329))

Rango = seq(0,1500,500)


library(ggplot2)
library(scales)

png(file="C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2021/stock_LCN/figures/figure_5.png",
    width=1200, height=700)

#windows()
 ggplot(df1, aes(x=Año, y=Rango))+
  geom_linerange(data=df1, aes(ymin=min,ymax=max),linetype=2,color="red")+
  geom_point(data=df1, aes(y=min),size=3,color="red")+
  geom_point(data=df1, aes(y=max),size=3,color="red")

 
 dev.off()  
  
  #xlab("AÑO") + ylab("RANGO DE CBA (ton)")+
  
  
  

