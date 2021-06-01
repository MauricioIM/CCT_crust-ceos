rm(list=ls())


df1 <- data.frame(Año=seq(2018,2021,1), min=c(4569,3173,2491,2647), 
                 max=c(5609,3999,3411,3462))

Rango = seq(0,6000,2000)


library(ggplot2)
library(scales)

png(file="C:/Users/mauricio.ibarra/Documents/IFOP/Proyectos/CTP/Orden_carpetas_crustáceos/2021/stock_CAS/figures/figure_5.png",
    width=1200, height=700)

#windows()
 ggplot(df1, aes(x=Año, y=Rango))+
  geom_linerange(data=df1, aes(ymin=min,ymax=max),linetype=2,color="red")+
  geom_point(data=df1, aes(y=min),size=3,color="red")+
  geom_point(data=df1, aes(y=max),size=3,color="red")

 
 dev.off()  
  
  #xlab("AÑO") + ylab("RANGO DE CBA (ton)")+
  
  
  

