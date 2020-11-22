rm(list=ls())
library(fields); library(maps); library(ncdf4);library(abind)
library(geosphere);library(mapdata);library(Hmisc)
setwd("~/Documents/CH4")
source('Function/get_geo.R')
source('Function/get_met.R')

#=== load global background ===
setwd("~/Documents/CH4/Data/TROPOMI/CPU_global")
ss1=load("Daily_CH4.Rdata")

convert=function(x)return((paste(x,collapse="-")))
strdate=as.character()
for(k in 1:dim(TROPOMI.date)[1])strdate=abind(strdate,convert(TROPOMI.date[k,]))
strdate=as.Date(strdate)

ind1=(TROPOMI.lon>=-180 & TROPOMI.lon<=180)
ind2=(TROPOMI.lat>=0 & TROPOMI.lat<=60)
bkgd=apply(TROPOMI.data[ind1,ind2,],c(3),mean,na.rm=T)
count=apply(!is.na(TROPOMI.data[ind1,ind2,]),c(3),sum,na.rm=T)
bkgd[count<=200]=NA
bkgd[bkgd>=1900]=NA

dev.new(width=5,height=2.8)
par(mai=c(0.3, 0.6, 0.25, 0.2), mgp=c(1.4, 0.4, 0), tcl=-0.2, ps=13)
plot(strdate,bkgd,pch=16,cex=0.5,axe=F,ylab="CH4 (ppb)")
axis(1,at=c(as.Date("2018-07-15"),as.Date("2019-01-15"),as.Date("2019-06-15")),labels=c("2018/07","2019/01","2019/06"))
axis(2)# minor.tick(nx=2, ny=5)
box()
bkgd_smooth=mov.avg(bkgd,1:length(bkgd),15,15)
lines(strdate, bkgd_smooth,col=2,lwd=2)
title("Global CH4 background",font.main=1,cex.main=1)
legend("topleft","0-60°N",bty="n",cex=0.8)

ind=(strdate<=as.Date("2019-04-30"))
adjusted_bkgd= bkgd_smooth-mean(bkgd_smooth[ind],na.rm=T)
dev.new(width=5,height=2.8)
par(mai=c(0.3, 0.6, 0.25, 0.2), mgp=c(1.4, 0.4, 0), tcl=-0.2, ps=13)
plot(strdate, adjusted_bkgd,pch=16,cex=0.5,axe=F,ylab="CH4 (ppb)",type="l",lwd=2)
axis(1,at=c(as.Date("2018-07-15"),as.Date("2019-01-15"),as.Date("2019-06-15")),labels=c("2018/07","2019/01","2019/06"))
axis(2)# minor.tick(nx=2, ny=5)
box()
abline(h=0,lty=2,v=as.Date("2019-04-30"))
title("Adjusted CH4 background",font.main=1,cex.main=1)

#=== load NA CH4 data ===
setwd("~/Documents/CH4/Data/TROPOMI/CPU_NorthAmerica_01x01/data_highQA")
ss=load("Daily_CH4.Rdata")
TROPOMI.data[TROPOMI.data<=1750]=NA

# for(t in 1:length(strdate)){
	# TROPOMI.data[,,t]= TROPOMI.data[,,t]-adjusted_bkgd[t]
# }

avg=apply(TROPOMI.data,c(1,2),mean,na.rm=T)
count=apply(!is.na(TROPOMI.data),c(1,2),sum,na.rm=T)
ind1=(TROPOMI.lon>=-130 & TROPOMI.lon<=-60)
ind2=(TROPOMI.lat>=5 & TROPOMI.lat<=55)

ind1=(TROPOMI.lon>=-120 & TROPOMI.lon<=-85)
ind2=(TROPOMI.lat>=10 & TROPOMI.lat<=35)

setwd("/Users/lu/Documents/CH4_Mexico/Figures/CH4_TROPOMI_0.1x0.1")
pdf(file="CH4_posteriri_Mexico.pdf",width=5,height=3)
par(mai=c(0.3, 0.6, 0.25, 0.2), mgp=c(1.6, 0.4, 0), tcl=-0.2, ps=13)
ap=avg;ap[count<=10]=NA
plot.field(ap[ind1,ind2], TROPOMI.lon[ind1], TROPOMI.lat[ind2],type="def",zlim=c(1800,1900),legend.mar=5)
title("CH4 a posteriori (2018-2019, seasonality adjusted)",font.main=1,cex.main=1)
dev.off()

#==========
avg_altitude=apply(TROPOMI.altitude,c(1,2),mean,na.rm=T)
x1=as.numeric(avg)
x2=as.numeric(avg_altitude)
MAT=na.omit(cbind(x1,x2))
fit<-lm(MAT[,1]~MAT[,2])

bias= TROPOMI.altitude*coef(fit)[2]
bias= TROPOMI.altitude*(-7/1000)
TROPOMI.data2 =  TROPOMI.data - bias

avg_corrected=apply(TROPOMI.data2,c(1,2),mean,na.rm=T)
ind1=(TROPOMI.lon>=-130 & TROPOMI.lon<=-60)
ind2=(TROPOMI.lat>=5 & TROPOMI.lat<=55)

ind1=(TROPOMI.lon>=-120 & TROPOMI.lon<=-85)
ind2=(TROPOMI.lat>=10 & TROPOMI.lat<=35)

pdf(file="CH4_posteriri_NA_altitude_adjusted.pdf",width=5,height=3)
par(mai=c(0.3, 0.6, 0.25, 0.2), mgp=c(1.6, 0.4, 0), tcl=-0.2, ps=13)
plot.field(avg_corrected[ind1,ind2], TROPOMI.lon[ind1], TROPOMI.lat[ind2],type="def",zlim=c(1800,1905),legend.mar=5)
title("CH4 a posteriori (2018-2019, bkgd/altitude adjusted)",font.main=1,cex.main=1)
dev.off()


# pdf(file="Altitude.pdf",width=5,height=3)
# par(mai=c(0.3, 0.6, 0.25, 0.2), mgp=c(1.6, 0.4, 0), tcl=-0.2, ps=13)
# plot.field(avg_altitude[ind1,ind2], TROPOMI.lon[ind1], TROPOMI.lat[ind2],type="def",zlim=c(0,4000),legend.mar=5)
# # title("CH4 a posteriori (2018-2019, bkgd/altitude adjusted)",font.main=1,cex.main=1)
# dev.off()
save(avg, count, avg_corrected, avg_altitude, TROPOMI.lon, TROPOMI.lat, file="TROPOMI_XCH4_01x01.Rdata"  )