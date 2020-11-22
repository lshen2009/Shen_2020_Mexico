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

pdf("Figure_2.pdf",width=7,height=2.5)
# dev.new(width=7,height=2.5)

m <- rbind(	       c(0, 1, 0.7,1),
		c(0,0.33,0.23,0.9),c(0.33,0.67,0.23,0.9),c(0.67,0.99,0.23,0.9),
		c(0,0.33,0.1,0.50),c(0.33,0.67,0.1,0.5),c(0.67,0.99,0.1,0.5)
         )

# dev.new(width=7.5,height=2.8)         
close.screen(all.screens = TRUE)
split.screen(m)

screen(3)
par(mai=c(0.1, 0.1, 0.2, 0.1), mgp=c(1.3, 0.2, 0), tcl=-0.2, ps=9)
ss=load("~/Documents/CH4_Mexico_archive/Data/flaring/Daily_flaring.Rdata")
ind1=(lon.out>=xlim[1] & lon.out <=xlim[2])
ind2=(lat.out >=ylim[1] & lat.out<=ylim[2])
avg=apply(daily_flaring[ind1,ind2,],c(1,2),sum,na.rm=T)/dim(daily_flaring)[3]
avg[avg<=0.05]=NA
zlim=c(-2,2)
ap=log10(avg);ap[ap>=zlim[2]]=zlim[2];ap[ap<=zlim[1]]=zlim[1]
image(lon.out[ind1],lat.out[ind2],ap,zlim=zlim,xaxt="n",yaxt="n",xlab="",ylab="",col=tim.colors(32))
box(); fun()
map("world",add=T)
title("(c) Flaring radiant heat",font.main=1,cex.main=1,line=+0.3)

screen(6)
par(mai=c(0.1, 0.1, 0.2, 0.1), mgp=c(1.3, 0.2, 0), tcl=-0.2, ps=9)
image.plot(NA,legend.only=T,col=tim.colors(32), zlim=zlim, horizontal = T, axis.args=list(mai=c(0.1, 0.1, 0.2, 0.1),mgp=c(0,0.2,0), tcl=-0.2))
text(x=0.5,y=-0.08,bquote("[ MW ]" )  )

screen(4)
par(mai=c(0.1, 0.1, 0.2, 0.1), mgp=c(1.3, 0.2, 0), tcl=-0.2, ps=9)
datafile=nc_open("~/Documents/TROPOMI/NO2/CPU_North_America/Step3_daily/Daily_data_qa_0.75.nc")
lon=ncvar_get(datafile,varid="lon")
lat=ncvar_get(datafile,varid="lat")
NO2 =ncvar_get(datafile,varid="NO2")#mol m-2
CF=6.02e23/(1e4)/1e15
NO2=NO2*CF
date=ncvar_get(datafile,varid="date")
nc_close(datafile)
NO2[122,40,]=NO2[122,40,]*1.5
avg=apply(NO2,c(1,2),mean,na.rm=T)
ind1=(lon>=xlim[1] & lon <=xlim[2])
ind2=(lat>=ylim[1] & lat<=ylim[2])
zlim=c(0,2.1)
ap=(avg[ind1,ind2]);ap[ap>=zlim[2]]=zlim[2];ap[ap<=zlim[1]]=zlim[1]
image(lon[ind1],lat[ind2],ap,zlim=zlim,xaxt="n",yaxt="n",xlab="",ylab="",col=tim.colors(32))
# plot.field(avg[ind1,ind2],lon[ind1],lat[ind2],legend.mar=5,type="def",zlim=c(0,4),mai=mai,ps=ps)
map("world",add=T)
box(); fun()
title(bquote("(d) TROPOMI "*NO[2]*""), font.main=1, cex.main=1,line=0.3)

screen(7)
par(mai=c(0.1, 0.1, 0.2, 0.1), mgp=c(1.3, 0.2, 0), tcl=-0.2, ps=9)
image.plot(NA,legend.only=T,col=tim.colors(32), zlim=zlim, horizontal = T, axis.args=list(mai=c(0.1, 0.1, 0.2, 0.1),mgp=c(0,0.2,0), tcl=-0.2))
text(x=0.92,y=-0.08,bquote("[ "*10^15*" molec "*cm^-2*" ]" )  )

dev.off()