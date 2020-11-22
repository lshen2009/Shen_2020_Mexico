rm(list=ls())
library(fields); library(maps); library(ncdf4);library(abind)
library(geosphere)
setwd("~/Documents")
source('Function/get_geo.R')
source('Function/get_met.R')
source('Function/read_met.R')
source('Function/read_method.R')

#========================================================
#========================================================
setwd("~/Documents/CH4_Mexico/Data")

dirnames=c("GC_format_orig","GC_format_offshoreX13%","GC_format_offshoreX13%_mexicox50%","GC_format_offshoreX13%_onshorex100%")

NP=c(-93.122,17.864)
# loc=find.lon.lat(NP[1],NP[2],lon,lat)
# lon[loc[1]]
# lat[loc[2]]
# Clusters[loc[1],loc[2]]
# sum(WW[,935])
NP_emis=39.7

result=NULL
ss=load("/Users/lu/Documents/CH4_Mexico/Figures/Figure 4/multiple_choices/all_figs.Rdata")
ilist=1
for(ilist in 1:length(all_figs)){
	met=all_figs[[ilist]]
	result=rbind(result, c(met$onshore, met$gamma, match(met$dirname, dirnames)))
}

colnames(result)=c("Prior","Posterior","Gamma","Inversions")
result2=cbind(#c(result[1,"Prior"],result[1:3,"Posterior"]+ NP_emis),
c(result[4,"Prior"],result[4:6,"Posterior"]+ NP_emis),
c(result[7,"Prior"],result[7:9,"Posterior"]+ NP_emis),
c(result[10,"Prior"],result[10:12,"Posterior"]+ NP_emis))

setwd("~/Documents/CH4_Mexico/Figures/Figure 4/multiple_choices")
dev.new(width=7.5, height=2.8)
# pdf("Figure_onshore.pdf",width=7.5, height=2.8)

m <- rbind(	       c(0,1,0.7,1),
		c(0,0.75,0,0.95),c(0.75,1.0,0,0.95),
		c(0.75,0.98,0,0.5)
		
         )

close.screen(all.screens = TRUE)
split.screen(m)
screen(1)
par(mai=c(0.2, 0.2, 0.25, 0.2), mgp=c(1.3, 0.4, 0), tcl=-0.2, ps=11)
mtext("Emissions in onshore regions using different priors and γ",side=3,line=0)

screen(2)
par(mai=c(0.1, 0.5, 0.15, 0), mgp=c(1.2, 0.2, 0), tcl=-0.2, ps=10)
centers=barplot(result2,main="", ylab="Emissions (Gg)", beside=TRUE,col=rep(c("#999999", "#E69F00", "#56B4E9", "#009E73"),7),bty="n",ylim=c(0,800),border=1,xlim=c(1.0,18))
box()
abline(v=c(5.5,10.5,15.5) ,lty=2)

emis_Aircraft=c(32000-7740,37000-7740, 42000-7740)*1000*24*365/1e9 #Gg a-1
arrows(21.5-5, emis_Aircraft[1], 21.5-5, emis_Aircraft[3], lwd=2, angle=90, code=3,length=0.04,col=2)
points(21.5-5, emis_Aircraft[2],pch=16,col=2)

emis_TROPOMI=c(20000,49000,78000)*1000*24*365/1e9 #Gg a-1
arrows(22.5-5, emis_TROPOMI[1], 22.5-5, emis_TROPOMI[3], lwd=2, angle=90, code=3,length=0.04,col=4)
points(22.5-5, emis_TROPOMI[2],pch=16,col=4)

# text(x=3,y=750,"Scarpelli (2020)")
# text(x=3,y=750,"R6 x 0.13")
# text(x=8,y=750,"R6 x 0.13")
text(x=8,y=750,"Mexico x 1.5")
# text(x=13,y=750,"R6 x 0.13")
text(x=13,y=750,"R5 x 2.0")
# text(x=22.5,y=750,"")
# mtext(bquote("(a) Standard"),side=3,line=0.2,cex=1.1)

screen(3)
par(mai=c(0.1, 0, 0.15, 0.0), mgp=c(1.2, 0.2, 0), tcl=-0.2, ps=9)
plot(NA,NA,xlim=c(0,1),ylim=c(0,1),axe=F,xlab="",ylab="")
#legend(x=0,y=1,c("Scarpelli","OffshoreX13%","OffshoreX13%_MexicoX152%","OffshoreX13%_onshoreX205%"),fill=c("#999999", "#E69F00", "#56B4E9", "#009E73"),bty="n")
legend(x=0,y=1,c("Prior","γ=0.10","γ=0.25","γ=0.50"),fill=c("#999999", "#E69F00", "#56B4E9", "#009E73"),bty="n",y.intersp=0.8)
legend(x=0,y=0.67,c("Aircraft","TROPOMI (Mass balance)"),col=c(2,4),text.col=c(2,4),bty="n",lty=1,pch=16,lwd=2,y.intersp=0.8,cex=0.9)

screen(4)
par(mai=c(0.1, 0.1, 0.1, 0.1), mgp=c(1.2, 0.2, 0), tcl=-0.2, ps=10)
plot(NA,NA,xlim=c(-97,-87),ylim=c(13,22),xlab="",ylab="",axe=F)
map("world",add=T)
ss=load("/Users/lu/Documents/CH4_Mexico/Data/Aircraft/coords_onshore.Rdata")
lines(lonlat[,1], lonlat[,2],col=2,lwd=2)
points(-93.122,17.864,cex=0.6,pch=16,col=2)

# dev.off()