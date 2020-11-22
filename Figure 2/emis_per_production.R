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

all_loc=rbind(
c(-99.4,-98.2,25.4,26.1),
c(-98.4,-97.5,21.8,22.7),
c(-98.0,-96.9,20.3,21.2),
c(-96.5,-95.5,18,19),
c(-94.5,-91.9,17.3,18.6),
c(-93.0,-91.7,18.6,20),
c(-104,-85,10,35)
)

methane_contents=c(0.90, 0.77, 0.77, 0.96, 0.70, 0.48, 0.60)
#=============== functions ====================
setwd("/Users/lu/Documents/CH4_Mexico/Figures/Figure 2")

all_percent=NULL
result=NULL

for(index in 1:7){
	
ss=load("~/Documents/CH4_Mexico_archive/Data/Well_production/CNH_historical/grid_gas_production.Rdata")#gas unit is 2.832x10^4 m3 per day
avg=apply(monthly_production,c(1,2),mean,na.rm=T)*365
loc= all_loc[index,]
# loc=c(-104,-85,10,35)
ind1=(lon.out>= loc[1] & lon.out<= loc[2])
ind2=(lat.out>= loc[3] & lat.out<= loc[4])
# production=sum(avg[ind1,ind2],na.rm=T)*1000/22.4*16#Gg a-1
# production=sum(avg[ind1,ind2],na.rm=T)*0.02832*1e6*1000/22.4*16/1e9*0.69#Gg a-1
production=sum(avg[ind1,ind2],na.rm=T)*0.02832*1e6*1000/22.4*16/1e9* methane_contents[index]#Gg a-1

ss=load("/Users/lu/Documents/CH4_Mexico/Data/Flaring/prior_posterior_emis.Rdata")
ind1=(lon>= loc[1] & lon <= loc[2])
ind2=(lat>= loc[3] & lat<= loc[4])
emis1=sum(prior_emis[ind1,ind2],na.rm=T)
emis2= sum(prior_emis[ind1,ind2]*grid_ratio[ind1,ind2],na.rm=T)

ratio1= emis1/production
ratio2= emis2/production
all_percent=rbind(all_percent, c(ratio1,ratio2))

result=rbind(result, c(production, emis1, emis2))
}


setwd("/Users/lu/Documents/CH4_Mexico/Figures/Figure 2")
# dev.new(width=5,height=2.8)
pdf("ratio.pdf",width=5,height=2.8)
par(mai=c(0.3, 0.6, 0.25, 0.2), mgp=c(1.4, 0.3, 0), tcl=-0.2, ps=11)
centers=barplot(t(all_percent*100),beside=TRUE,col=c("#E69F00", "#56B4E9"),ylim=c(0, 16),axe=F,ylab="Relative leakage (%)")
x=apply(centers,c(2),mean)
axis(1,at=x,labels=c("R1","R2","R3","R4","R5","R6"))
# axis(2,at=c(0,0.005,0.01,0.015),labels=c("0%","0.005%","0.01%","0.015%"))
axis(2)
box()
legend("topleft", c("Prior","Posterior"), cex=1.0, bty="n",fill=c("#E69F00","#56B4E9"),border=1)
title("Relative leakage in R1-R6",font.main=1,cex.main=1)
dev.off()