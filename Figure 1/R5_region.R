rm(list=ls())
library(fields); library(maps); library(ncdf4);library(abind)
library(geosphere);library(mapdata);library(Hmisc)
setwd("~/Documents/CH4")
source('Function/get_geo.R')
source('Function/get_met.R')

plot.box=function(loc,lwd=lwd){
arrows(x0=loc[1],y0=loc[3],x1=loc[2],code=0,col=1,lwd=lwd)
arrows(x0=loc[1],y0=loc[4],x1=loc[2],code=0,col=1,lwd=lwd)
arrows(x0=loc[1],y0=loc[3],y1=loc[4],code=0,col=1,lwd=lwd)
arrows(x0=loc[2],y0=loc[3],y1=loc[4],code=0,col=1,lwd=lwd)
}

plot.field=function(spdata, lon.map, lat.map, type=NULL, same=FALSE, zlim=NULL, col=NULL, nlevel=32, mai=c(0.2, 0.2, 0.2, 0.2), mgp=c(1.4, 0.5, 0), tcl=-0.3, ps=12, legend.mar=3, legend.width=1.2, xaxt="n", yaxt="n", map.region='world', Pacific.centric=FALSE, custom.breaks=NULL, map.xlim=NULL, map.ylim=NULL) {
   
   # This function plots spatial field data using function "image.plot" from package "fields".
   # Packages "fields" and "maps" must be pre-loaded.
   # "spdata" is a matrix or an array of spatial data.
   # "type" is a vector of intended types of presentation, as explained below.
   # Five types of presentation are supported:
   # 1. "sign": variable that has both positive and negative values, good for comparing signs of correlations/effects/deviations.
   # 2. "frac": variable that is a fraction or percentage, good for comparing proportions.
   # 3. "abs": variable that has absolute values only, good for comparing magtitudes.
   # 4. "def": user-defined scale and z-limits. (If chosen, must define the same same scale/limits for all plots.) 
   # 5. Default: no specification for display.
   # If you want all plots to have the same type, simply enter a string scalar for "type".
   
	if (lat.map[length(lat.map)] < lat.map[1]) {
		# 'lat.map' is in decreasing order. Need to reverse it to proceed further.
		lat.map = rev(lat.map)
		spdata = spdata[,length(lat.map):1]
	}	
	par(mai=mai, mgp=mgp, tcl=tcl, ps=ps)
	if (is.null(custom.breaks)) {
		if (is.null(type)) {
			zlim = c(min(na.omit(as.vector(spdata))), max(na.omit(as.vector(spdata))))
		} else if (type == 'sign') {
			if (is.null(zlim)) zlim = c(-max(abs(na.omit(as.vector(spdata)))), max(abs(na.omit(as.vector(spdata))))) else zlim = zlim
			if (is.null(col)) col = rwb.colors(nlevel)
		} else if (type == 'frac') {
			zlim = c(0,1)
		} else if (type == 'abs') {
			zlim = c(0, max(na.omit(as.vector(spdata))))
		} else if (type == 'def') {
			zlim = zlim
		} else {
			zlim = c(min(na.omit(as.vector(spdata))), max(na.omit(as.vector(spdata))))
		}
		if (is.null(col)) col = tim.colors(nlevel)
		spdata[which(spdata > zlim[2])] = zlim[2]
		spdata[which(spdata < zlim[1])] = zlim[1]
		image.plot(lon.map, lat.map, spdata, zlim=zlim, xlab='', ylab='', axis.args=list(mgp=c(0,0.5,0), tcl=-0.2), legend.mar=legend.mar, legend.width=legend.width, col=col, nlevel=nlevel, xaxt=xaxt, yaxt=yaxt)
	} else {
		nbreaks = length(custom.breaks)
		nlevel = nbreaks - 1
		if (!is.null(type)) {
			if (type == 'sign') col = rwb.colors(nlevel) else col = tim.colors(nlevel)
		} else col = tim.colors(nlevel)
		zval.breaks = 0:nlevel
		zval.center = 0:(nlevel - 1) + 0.5
		spdata.new = spdata
		for (n in 1:nlevel) spdata.new[which(spdata >= custom.breaks[n] & spdata <= custom.breaks[n + 1])] = zval.center[n]
		spdata.new[which(spdata < custom.breaks[1])] = zval.center[1]
		spdata.new[which(spdata > tail(custom.breaks, 1))] = tail(zval.center, 1)
		spdata = spdata.new
		zlim = c(0, nlevel)
		image.plot(lon.map, lat.map, spdata, zlim=zlim, xlab='', ylab='', axis.args=list(mgp=c(0,0.5,0), tcl=-0.2), legend.mar=legend.mar, legend.width=legend.width, col=col, nlevel=nlevel, xaxt=xaxt, yaxt=yaxt, breaks=zval.breaks, lab.breaks=custom.breaks)
	}
	if (Pacific.centric) map('world2', add=TRUE, xlim=map.xlim, ylim=map.ylim,col=1,lwd=0.5) else map(map.region, add=TRUE, xlim=map.xlim, ylim=map.ylim,col=1,lwd=0.5)
}


ss=load("/Users/lu/Documents/CH4_Mexico/Data/Emis_Tia/emis_per_area.Rdata")
ind11=(lon.in>=-93.0 & lon.in <=-91.7)
ind21=(lat.in>=18.6 & lat.in <=19.8)
emis_per_area[ind11,ind21]=emis_per_area[ind11,ind21]*0.1


#=== load global background ===
ss=load("/Users/lu/Documents/CH4_Mexico/Figures/Figure 1/TROPOMI_XCH4_01x01.Rdata")
# xlim=c(-120,-100);ylim=c(30,45)
# ind1=(TROPOMI.lon>=xlim[1] & TROPOMI.lon<=xlim[2])
# ind2=(TROPOMI.lat>=ylim[1] & TROPOMI.lat<=ylim[2])
# plot.field(avg_corrected[ind1,ind2], TROPOMI.lon[ind1], TROPOMI.lat[ind2],type="def",zlim=c(1830,1890),legend.mar=6,mai=c(0.1,0.2,0.5,0.2),ps=10)
# map("state",add=T)


#======================
setwd("/Users/lu/Documents/CH4_Mexico/Figures/Figure 1")

pdf(file="CH4_posteriri_Mexico.pdf",width=6.5,height=2.5)
# dev.new(width=6.5,height=2.5)
xlim=c(-106,-86);ylim=c(14,30)

par(mar=c(3,2,1,3))
par(mfrow=c(1,2))

ind1=(TROPOMI.lon>=xlim[1] & TROPOMI.lon<=xlim[2])
ind2=(TROPOMI.lat>=ylim[1] & TROPOMI.lat<=ylim[2])
ap= avg_corrected;ap[count<=1]=NA
plot.field(ap[ind1,ind2], TROPOMI.lon[ind1], TROPOMI.lat[ind2],type="def",zlim=c(1830,1890),legend.mar=6,mai=c(0.1,0.2,0.5,0.2),ps=10)
title("(a) TROPOMI XCH4,  elevation corrected",font.main=1,cex.main=0.9,line=+0.3)
plot.box(c(-99.4,-98.2,25.4,26.1),lwd=0.5)
plot.box(c(-98.4,-97.5,21.8,22.7),lwd=0.5)
plot.box(c(-98.0,-96.9,20.3,21.2),lwd=0.5)
plot.box(c(-96.5,-95.5,18,19),lwd=0.5)
plot.box(c(-94.5,-91.9,17.5,18.6),lwd=0.5)
plot.box(c(-93.0,-91.7,18.6,19.8),lwd=0.5)
plot.box(c(-100+0.5,-99+0.3,19.2,19.8),lwd=0.5)
text(x=-84.5,y=31,"ppb",xpd=NA,cex=0.9)
text(x=-101,y=15,"May 2018 - Dec 2019",cex=0.8)

ind1=(lon.in>=xlim[1] & lon.in <=xlim[2])
ind2=(lat.in>=ylim[1] & lat.in <=ylim[2])
ap= emis_per_area
ap[ap==0]=NA
plot.field(log10(ap[ind1,ind2]),lon.in[ind1],lat.in[ind2],legend.mar=5,type="def",zlim=c(7,13),mai=c(0.1,0.2,0.5,0.2),ps=10)
title("(b) Prior oil/gas emissions in 2015",cex.main=0.9,font.main=1,line=+0.3)
plot.box(c(-99.4,-98.2,25.4,26.1),lwd=0.5)
plot.box(c(-98.4,-97.5,21.8,22.7),lwd=0.5)
plot.box(c(-98.0,-96.9,20.3,21.2),lwd=0.5)
plot.box(c(-96.5,-95.5,18,19),lwd=0.5)
plot.box(c(-94.5,-91.9,17.3,18.6),lwd=0.5)
plot.box(c(-93.0,-91.7,18.6,20),lwd=0.5)
plot.box(c(-100+0.5,-99+0.3,19.2,19.8),lwd=0.5)

text(x=-100.2,y=25.8,"R1",cex=0.8)
text(x=-96.7,y=22.3,"R2",cex=0.8)
text(x=-96.1,y=20.7,"R3",cex=0.8)
text(x=-95.9,y=17.3,"R4",cex=0.8)
text(x=-93.1,y=16.8,"R5",cex=0.8)
text(x=-92.4,y=20.3,"R6",cex=0.8)
text(x=-99.1,y=18.6,"R7",cex=0.8)

mtext("TROPOMI XCH4 and Oil/gas emissions in eastern Mexico",side=3,outer=T,line=-1.5)

dev.off()


lonlat=rbind(c(-97.81278, 22.27944),
c(-98.48250, 26.01528),
c(-93.19209, 17.90342),
c(-92.48750, 17.87500),
c(-94.35816, 18.09292),
c(-94.04417, 18.09056),
c(-96.40139, 18.82250),
c(-93.12111, 17.85861),
c(-97.47167, 20.52111))
 
 
ss=load("/Users/lu/Documents/CH4_Mexico/Figures/Figure 1/TROPOMI_XCH4_01x01.Rdata") 
pdf(file="R5.pdf",width=4.5,height=3)
xlim=c(-94.5,-91.9);ylim=c(17.3,18.6)
ind1=(TROPOMI.lon>=xlim[1]-1 & TROPOMI.lon<=xlim[2]+1)
ind2=(TROPOMI.lat>=ylim[1]-1 & TROPOMI.lat<=ylim[2]+1)
ap= avg_corrected;ap[count<=1]=NA
plot.field(ap[ind1,ind2], TROPOMI.lon[ind1], TROPOMI.lat[ind2],type="def",zlim=c(1830,1890),legend.mar=6,mai=c(0.1,0.2,0.5,0.2),ps=12)
title("(a) TROPOMI XCH4 in R5,  elevation corrected",font.main=1,cex.main=0.9,line=+0.3)
# lonlat=all_facility[[4]]
points(lonlat[,1],lonlat[,2],cex=0.5)
plot.box(c(-94.5,-91.9,17.3,18.6),lwd=0.5)

dev.off()