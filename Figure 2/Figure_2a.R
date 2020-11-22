rm(list=ls())
library(fields); library(maps); library(ncdf4);library(abind)
setwd("~/Documents")
source('Function/get_geo.R')
source('Function/get_met.R')

plot.box=function(loc,col=1,lwd=1){
arrows(x0=loc[1],y0=loc[3],x1=loc[2],code=0,col=col,lwd= lwd)
arrows(x0=loc[1],y0=loc[4],x1=loc[2],code=0,col=col,lwd= lwd)
arrows(x0=loc[1],y0=loc[3],y1=loc[4],code=0,col=col,lwd= lwd)
arrows(x0=loc[2],y0=loc[3],y1=loc[4],code=0,col=col,lwd= lwd)
}

fun=function(){
	
plot.box(c(-99.4,-98.2,25.4,26.1),lwd=0.5)
plot.box(c(-98.4,-97.5,21.8,22.7),lwd=0.5)
plot.box(c(-98.0,-96.9,20.3,21.2),lwd=0.5)
plot.box(c(-96.5,-95.5,18,19),lwd=0.5)
plot.box(c(-94.5,-91.9,17.5,18.6),lwd=0.5)
plot.box(c(-93.0,-91.7,18.6,19.8),lwd=0.5)
plot.box(c(-100+0.5,-99+0.3,19.2,19.8),lwd=0.5)	
	
}

#=============== functions ====================
setwd("/Users/lu/Documents/CH4_Mexico/Figures/Figure 2")
xlim=c(-106,-86);ylim=c(14,30)
mai=c(0.2,0.2,0.5,0.2);ps=15

pdf("Figure_2a.pdf",width=7,height=2.5)
# dev.new(width=7,height=2.5)

m <- rbind(	       c(0, 0.67, 0.7,1),
		c(0,0.33,0.23,0.9),c(0.33,0.67,0.23,0.9),
		c(0,0.33,0.1,0.50),c(0.33,0.67,0.1,0.5)
         )

# dev.new(width=7.5,height=2.8)         
close.screen(all.screens = TRUE)
split.screen(m)
screen(1)
par(mai=c(0.2, 0.2, 0.25, 0.2), mgp=c(1.3, 0.4, 0), tcl=-0.2, ps=10)
mtext(bquote("Oil/gas production, flaring radiant heat and "*NO[2]*" in easterm Mexico"),side=3,line=0)

screen(2)
par(mai=c(0.1, 0.1, 0.2, 0.1), mgp=c(1.3, 0.2, 0), tcl=-0.2, ps=9)
ss=load("~/Documents/CH4_Mexico_archive/Data/Well_production/CNH_historical/grid_oil_production.Rdata")#gas unit is 10^9 m3 per month; oil unit: 10^3 barrels month

ind1=(lon.out>=xlim[1] & lon.out<=xlim[2])
ind2=(lat.out>=ylim[1] & lat.out<=ylim[2])
avg=apply(monthly_production,c(1,2),mean)
zlim=c(0,3)
ap=log10(avg[ind1,ind2]);ap[ap>=zlim[2]]=zlim[2];ap[ap<=zlim[1]]=zlim[1]
image(lon.out[ind1],lat.out[ind2],ap,zlim=zlim,xaxt="n",yaxt="n",xlab="",ylab="",col=tim.colors(32))
box(); fun()
map("world",add=T)
title("(a) Oil production",font.main=1,cex.main=1,line=+0.3)

screen(4)
par(mai=c(0.1, 0.1, 0.2, 0.1), mgp=c(1.3, 0.2, 0), tcl=-0.2, ps=9)
image.plot(NA,legend.only=T,col=tim.colors(32), zlim=zlim, horizontal = T, axis.args=list(mai=c(0.1, 0.1, 0.2, 0.1),mgp=c(0,0.2,0), tcl=-0.2))
# text(x=0.1,y=-0.08,bquote("[ "*10^3~barrels*day^-1*" ]" )  )
text(x=0.1,y=-0.08, "[ Mbd ]")
#------
screen(3)
par(mai=c(0.1, 0.1, 0.2, 0.1), mgp=c(1.3, 0.2, 0), tcl=-0.2, ps=9)
ss=load("~/Documents/CH4_Mexico_archive/Data/Well_production/CNH_historical/grid_gas_production.Rdata")#gas unit is 10^9 m3 per month; oil unit: 10^3 barrels month

ind1=(lon.out>=xlim[1] & lon.out<=xlim[2])
ind2=(lat.out>=ylim[1] & lat.out<=ylim[2])
avg=apply(monthly_production,c(1,2),mean)
zlim=c(0,3)
ap=log10(avg[ind1,ind2]);ap[ap>=zlim[2]]=zlim[2];ap[ap<=zlim[1]]=zlim[1]
image(lon.out[ind1],lat.out[ind2],ap,zlim=zlim,xaxt="n",yaxt="n",xlab="",ylab="",col=tim.colors(32))
box(); fun()
map("world",add=T)
title("(b) Gas production",font.main=1,cex.main=1,line=+0.3)

screen(5)
par(mai=c(0.1, 0.1, 0.2, 0.1), mgp=c(1.3, 0.2, 0), tcl=-0.2, ps=9)
image.plot(NA,legend.only=T,col=tim.colors(32), zlim=zlim, horizontal = T, axis.args=list(mai=c(0.1, 0.1, 0.2, 0.1),mgp=c(0,0.2,0), tcl=-0.2))
# text(x=0.5,y=-0.08,bquote("[ "*10^6~m^3*month^-1*" ]" )  )
text(x=0.5,y=-0.08,"[ MMcfd ]")

dev.off()
