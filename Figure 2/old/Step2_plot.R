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

#=============== functions ====================
setwd("/Users/lu/Documents/CH4_Mexico/Figures/Figure 2")
xlim=c(-106,-86);ylim=c(14,31)
mai=c(0.2,0.2,0.5,0.2);ps=15

# pdf("Flaring.pdf",width=7.5,height=2.5)

dev.new(width=7.5,height=2.5)
par(mar=c(3,2,1,3))
par(mfrow=c(1,3))

ss=load("~/Documents/CH4_Mexico/Data/Well_production/CNH_historical/grid_gas_production.Rdata")
ind1=(lon.out>=xlim[1] & lon.out<=xlim[2])
ind2=(lat.out>=ylim[1] & lat.out<=ylim[2])
avg=apply(monthly_production,c(1,2),mean)
plot.field(log10(avg[ind1,ind2]),lon.out[ind1], lat.out[ind2],type="def",zlim=c(0,3),mai=mai,ps=ps)
plot.box(c(-99.4,-98.2,25.4,26.1),lwd=0.5)
plot.box(c(-98.4,-97.5,21.8,22.7),lwd=0.5)
plot.box(c(-98.0,-96.9,20.3,21.2),lwd=0.5)
plot.box(c(-96.5,-95.5,18,19),lwd=0.5)
plot.box(c(-94.5,-91.9,17.5,19.8),lwd=0.5)
plot.box(c(-100+0.5,-99+0.3,19.2,19.8),lwd=0.5)
title("(a) Gas production",font.main=1,cex.main=1,line=+0.3)


ss=load("~/Documents/CH4_Mexico/Data/flaring/Daily_flaring.Rdata")
ind1=(lon.out>=xlim[1] & lon.out <=xlim[2])
ind2=(lat.out >=ylim[1] & lat.out<=ylim[2])
avg=apply(daily_flaring[ind1,ind2,],c(1,2),sum,na.rm=T)/dim(daily_flaring)[3]*365
avg[avg<=10]=NA
plot.field(log10(avg),lon.out[ind1],lat.out[ind2],type="def",zlim=c(0,4),mai=mai,ps=ps)
plot.box(c(-99.4,-98.2,25.4,26.1),lwd=0.5)
plot.box(c(-98.4,-97.5,21.8,22.7),lwd=0.5)
plot.box(c(-98.0,-96.9,20.3,21.2),lwd=0.5)
plot.box(c(-96.5,-95.5,18,19),lwd=0.5)
plot.box(c(-94.5,-91.9,17.5,19.8),lwd=0.5)
plot.box(c(-100+0.5,-99+0.3,19.2,19.8),lwd=0.5)
title("(b) Flaring radiant heat",font.main=1,cex.main=1,line=+0.3)

datafile=nc_open("~/Documents/TROPOMI/NO2/CPU_North_America/Step3_daily/long_term/Daily_data.nc")
lon=ncvar_get(datafile,varid="lon")
lat=ncvar_get(datafile,varid="lat")
NO2 =ncvar_get(datafile,varid="NO2")#mol m-2
CF=6.02e23/(1e4)/1e15
NO2=NO2*CF
date=ncvar_get(datafile,varid="date")
nc_close(datafile)
avg=apply(NO2,c(1,2),mean,na.rm=T)
ind1=(lon>=xlim[1] & lon <=xlim[2])
ind2=(lat>=ylim[1] & lat<=ylim[2])
plot.field(avg[ind1,ind2],lon[ind1],lat[ind2],legend.mar=5,type="def",zlim=c(0,4),mai=mai,ps=ps)
map("world",add=T)
# title(bquote("(b) TROPOMI "*NO[2]*""),font.main=1,cex.main=1,line=0.3)
plot.box(c(-99.4,-98.2,25.4,26.1),lwd=0.5)
plot.box(c(-98.4,-97.5,21.8,22.7),lwd=0.5)
plot.box(c(-98.0,-96.9,20.3,21.2),lwd=0.5)
plot.box(c(-96.5,-95.5,18,19),lwd=0.5)
plot.box(c(-94.5,-91.9,17.5,19.8),lwd=0.5)
plot.box(c(-100+0.5,-99+0.3,19.2,19.8),lwd=0.5)
# mtext(bquote("column density ("*10^15*" molec "*cm^-2*")"),side=4, line=3.5,cex=0.9)
title("(c) TROPOMI NO2", font.main=1, cex.main=1,line=0.3)

mtext("Gas production, flaring and NO2 in eastern Mexico",side=3,outer=T,line=-2,cex=0.75)
dev.off()