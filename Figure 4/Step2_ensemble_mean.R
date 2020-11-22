rm(list=ls())
library(fields); library(maps); library(ncdf4);library(abind)
library(geosphere)
setwd("~/Documents")
source('Function/get_geo.R')
source('Function/get_met.R')
source('Function/read_met.R')
source('Function/read_method.R')

ss=load("/Users/lu/Documents/CH4_Mexico/Figures/Figure 4/multiple_choices/all_figs.Rdata")

dirnames=c("GC_format_orig","GC_format_offshoreX13%","GC_format_offshoreX13%_mexicox50%","GC_format_offshoreX13%_onshorex100%")

result=NULL
all_figa=NULL
all_figb=NULL
figa_AK=NULL
all_figa_error=NULL
all_figb_error=NULL
all_anthro_emis=NULL
all_natural_emis=NULL
for(ilist in 1:length(all_figs)){
	met=all_figs[[ilist]]
	if(met$dirname=="GC_format_orig")next
	all_figa=abind(all_figa, met$figa, along=3)
	all_figa_error=abind(all_figa_error, t(met$figa_error), along=3)
	all_figb=abind(all_figb, met$figb, along=3)	
	all_figb_error=abind(all_figb_error, t(met$figb_error), along=3)
	all_anthro_emis =rbind(all_anthro_emis, met$anthro_emis)
	all_natural_emis =rbind(all_natural_emis, met$nature_emis)	
	figa_AK=rbind(figa_AK, met$figa_AK)
}

#=== replace figa and figb with "GC_format_offshoreX13%__0.25"
met=all_figs[["GC_format_offshoreX13%__0.25"]]
figa=met$figa
figb=met$figb

set.seed(10)

ind=sample(1:9, 1000, replace=T)
avg=all_figa;error=all_figa_error
all_sample=array(NA,c(2,7,1000))
for(i in 1:2){
	for(j in 1:7){
		for(k in 1:1000){
			all_sample[i,j,k]=rnorm(1, mean=avg[i,j,ind[k]], sd= error[i,j,ind[k]])
		}
	}
}
figa_error=apply(all_sample,c(1,2),sd)

avg=all_figb;error=all_figb_error
all_sample=array(NA,c(2,7,1000))
for(i in 1:2){
	for(j in 1:7){
		for(k in 1:1000){
			all_sample[i,j,k]=rnorm(1, mean=avg[i,j,ind[k]], sd= error[i,j,ind[k]])
		}
	}
}
figb_error=apply(all_sample,c(1,2),sd)

avg=all_anthro_emis[,1];error= all_anthro_emis[,2]
all_sample=array(NA,c(1000))
for(k in 1:1000){
	all_sample[k]=rnorm(1, mean=avg[ind[k]],sd=error[ind[k]])
}
c(mean(all_sample),sd(all_sample))/1000


avg= all_natural_emis[,1];error= all_natural_emis[,2]
all_sample=array(NA,c(1000))
for(k in 1:1000){
	all_sample[k]=rnorm(1, mean=avg[ind[k]],sd=error[ind[k]])
}
c(mean(all_sample),sd(all_sample))/1000
# 1.48+/0.09

#====================================================================
#====================================================================
pdf(file="/Users/lu/Documents/CH4_Mexico/Figures/Figure 4/Figure_4_ensemble.pdf",width=6.5,height=5)

# dev.new(width=6.5, height=5)
par(mar=c(3,2,1,3))
par(mfrow=c(2,1))

par(mai=c(0.3, 0.6, 0.25, 0.2), mgp=c(1.2, 0.2, 0), tcl=-0.2, ps=10)
barCenters=barplot(height= figa,main="", ylab="Emissions (Gg)", beside=TRUE,col=rep(c("#E69F00", "#56B4E9"),7),bty="n",ylim=c(0,2300),border=1)
means= figa[2,];standardErrors= figa_error[2,]
arrows(barCenters[2,], means-standardErrors*2, barCenters[2,], means+standardErrors*2, lwd=1, angle=90, code=3,length=0.02)
box()
title("(a) Prior and Posterior emissions from different sectors",font.main=1,cex.main=1)

par(mai=c(0.3, 0.6, 0.25, 0.2), mgp=c(1.2, 0.2, 0), tcl=-0.2, ps=10)
barCenters=barplot(height= figb,main="", ylab="Emissions (Gg)", beside=TRUE,col=rep(c("#E69F00", "#56B4E9"),7),bty="n",ylim=c(0,1300),border=1)
means= figb[2,];standardErrors= figb_error[2,]
arrows(barCenters[2,], means-standardErrors*2, barCenters[2,], means+standardErrors*2, lwd=1, angle=90, code=3,length=0.02)

legend("topright", c("Prior","Posterior"), cex=1.0, bty="n",fill=c("#E69F00","#56B4E9"),border=1)
# legend(x=18,y=1050, c("  Aircraft"), cex=1.0, bty="n",pch=16,text.col=2,border=1,col=2)
# legend(x=18,y=950, c("  Flaring"), cex=1.0, bty="n",pch=16,text.col=4,border=1,col=4)
box()
title("(b) Prior and posterior emissions in major oil/gas basins",font.main=1,cex.main=1)

# Flaring_emis=c(10.2, 15.2, 20.3)*3
# arrows(18.25, 30.6, 18.25, 60.9, lwd=1, angle=90, code=3,length=0.02,col=2)

# Flaring_emis=c(18, 90)
# arrows(15.25, 18, 15.25, 90, lwd=1.5, angle=90, code=3,length=0.02,col=4)

dev.off()