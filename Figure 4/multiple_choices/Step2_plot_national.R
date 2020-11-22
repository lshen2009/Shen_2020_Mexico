rm(list=ls())
library(fields); library(maps); library(ncdf4);library(abind)
library(geosphere)
setwd("~/Documents")
source('Function/get_geo.R')
source('Function/get_met.R')
source('Function/read_met.R')
source('Function/read_method.R')

plot.box=function(loc,lwd=lwd){
arrows(x0=loc[1],y0=loc[3],x1=loc[2],code=0,col=2,lwd=lwd)
arrows(x0=loc[1],y0=loc[4],x1=loc[2],code=0,col=2,lwd=lwd)
arrows(x0=loc[1],y0=loc[3],y1=loc[4],code=0,col=2,lwd=lwd)
arrows(x0=loc[2],y0=loc[3],y1=loc[4],code=0,col=2,lwd=lwd)
}

#========================================================
#========================================================
setwd("~/Documents/CH4_Mexico/Data")

dirnames=c("GC_format_orig","GC_format_offshoreX13%","GC_format_offshoreX13%_mexicox50%","GC_format_offshoreX13%_onshorex100%")

result=NULL
ss=load("/Users/lu/Documents/CH4_Mexico/Figures/Figure 4/multiple_choices/all_figs.Rdata")
ilist=1
for(ilist in 1:length(all_figs)){
	met=all_figs[[ilist]]
	result=rbind(result, c(met$figa[1,1], met$figa[2,1], met$gamma, match(met$dirname, dirnames)))
}

colnames(result)=c("Prior","Posterior","Gamma","Inversions")
result2=cbind(#c(result[1,"Prior"],result[1:3,"Posterior"]),
c(result[4,"Prior"],result[4:6,"Posterior"]),
c(result[7,"Prior"],result[7:9,"Posterior"]),
c(result[10,"Prior"],result[10:12,"Posterior"]))


# dev.new(width=5, height=2.8)
# par(mai=c(0.3, 0.6, 0.25, 0.2), mgp=c(1.2, 0.2, 0), tcl=-0.2, ps=10)
# centers=barplot(result2,main="", ylab="Emissions (Gg)", beside=TRUE,col=rep(c("#999999", "#E69F00", "#56B4E9", "#009E73"),7),bty="n",ylim=c(0,2300),border=1,xlim=c(1.0,20))
# box()
# abline(v=c(5.5,10.5,15.5) ,lty=2)
# h=2200
# text(x=3,y=h,"(A)")
# text(x=8,y=h,"(B)")
# text(x=13,y=h,"(C)")
# text(x=18,y=h,"(D)")



setwd("~/Documents/CH4_Mexico/Figures/Figure 4/multiple_choices")
dev.new(width=7.5, height=2.8)
# pdf("Figure_EastMexico.pdf",width=7.5, height=2.8)

m <- rbind(	       c(0,1,0.7,1),
		c(0,0.75,0,0.95),c(0.75,1.0,0,0.95),
		c(0.75,0.98,0,0.5)
		
         )

close.screen(all.screens = TRUE)
split.screen(m)
screen(1)
par(mai=c(0.2, 0.2, 0.25, 0.2), mgp=c(1.3, 0.4, 0), tcl=-0.2, ps=11)
mtext("Emissions in eastern Mexico using different priors and γ",side=3,line=0)

screen(2)
par(mai=c(0.1, 0.5, 0.15, 0), mgp=c(1.2, 0.2, 0), tcl=-0.2, ps=10)
centers=barplot(result2,main="", ylab="Emissions (Gg)", beside=TRUE,col=rep(c("#999999", "#E69F00", "#56B4E9", "#009E73"),7),bty="n",ylim=c(0,2400),border=1,xlim=c(1.0,15))
box()
abline(v=c(5.5,10.5) ,lty=2)

#text(x=3,y=2300,"Scarpelli (2020)")
# text(x=3,y=2300,"R6 x 0.13")
# text(x=8,y=2300,"R6 x 0.13")
text(x=8,y=2300,"Mexico x 1.5")
# text(x=13,y=2300,"R6 x 0.13")
text(x=13,y=2300,"R5 x 2.0")

screen(3)
par(mai=c(0.1, 0, 0.15, 0.1), mgp=c(1.2, 0.2, 0), tcl=-0.2, ps=9)
plot(NA,NA,xlim=c(0,1),ylim=c(0,1),axe=F,xlab="",ylab="")
legend(x=0,y=1,c("Prior","γ=0.10","γ=0.25","γ=0.50"),fill=c("#999999", "#E69F00", "#56B4E9", "#009E73"),bty="n",y.intersp=0.8)

screen(4)
par(mai=c(0.1, 0.1, 0.1, 0.1), mgp=c(1.2, 0.2, 0), tcl=-0.2, ps=10)
plot(NA,NA,xlim=c(-106,-86),ylim=c(13,31),xlab="",ylab="",axe=F)
map("world",add=T)
plot.box(c(-102,-86.5,15,30),lwd=1)

# dev.off()