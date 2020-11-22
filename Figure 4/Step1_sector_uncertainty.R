rm(list=ls())
library(fields); library(maps); library(ncdf4);library(abind)
library(geosphere)
setwd("~/Documents")
source('Function/get_geo.R')
source('Function/get_met.R')
source('Function/read_met.R')
source('Function/read_method.R')

#===== define functions =====
filename="/Users/lu/Documents/CH4_Mexico_archive/Data/Emissions_new/emis_Flexgrid/create_new_mask/East_mexico_clusters_v2.nc"
datafile=nc_open(filename)
lon=ncvar_get(datafile,varid="lon")
lat=ncvar_get(datafile,varid="lat")
Clusters=ncvar_get(datafile,varid="Clusters")
nc_close(datafile)

setwd("/Users/lu/Documents/CH4_Mexico/Data/GC_format_offshoreX13%/inversion")

datafile=nc_open("inversion_result.nc")
all_part1 =ncvar_get(datafile,varid="all_part1")
all_part2 =ncvar_get(datafile,varid="all_part2")
ratio =ncvar_get(datafile,varid="ratio")
nc_close(datafile)

emis_error=rep(0.5^2,1207)
inv_Sa =diag(1/emis_error)
gamma=0.25

ratio=solve(gamma*all_part1+inv_Sa)%*%(gamma*all_part2)+1
temp=(ratio-1)^2/emis_error
S_posterior= solve(gamma*all_part1+inv_Sa)
A=diag(1207)-S_posterior%*%inv_Sa
sum(diag(A))

ss=load("/Users/lu/Documents/CH4_Mexico/Data/GC_format_offshoreX13%/Flexgrid/WW.Rdata")
emis_by_sector=array(NA,c(10,1207))

cal_error=function(loc){
	WW2=array(0,dim(WW))
	WW2[,loc]=WW[,loc]
	S_red= WW2%*%S_posterior%*%t(WW2)
	error=cbind(sqrt(diag(WW2%*%diag(emis_error)%*%t(WW2))),sqrt(diag(S_red)))
	colnames(error)=c("prior","posterior")
	return(error)
}

fun=function(){
	# ind=which.max(avg)
	# return(error[ind])
	f=avg/sum(avg)
	return(sum(f*error))
}
#Coal, 1B
avg=c(275.6,2.5);error=c(37.6, 66.0); fun()
#Oil/gas, 1B2
avg=c(185.1, 317.6, 251.9);error=c(56.5, 39.7, 50);fun()
#Landfill
avg=c(607.4,87.7,87.7,4.2,0.1,22.1);error=c(3.5,2.4,0.8,104.4,6.5,97.6);fun()
#Livestock, 3A1, 3A2
avg=c(1790, 43.6, 43.5, 31.4, 284.6, 148.3,10.4);error=c(6.5,9.7,10.1,10,6.4,7.3,10);fun()
#Wastewater, 4D
avg=c(133.1,596.7);error=c(7.6,5.2);fun()
#Rice, 3C7
avg=c(5.9);error=c(55.8);fun()
#Biomass, 3C1
avg=c(5.8, 24.0, 0.9);error=c(50,241.1,50);fun()

# 1 "Oil_gas"    "Coal"       "Landfill"   "Livestock"  "Wastewater"
# 6 "Wetland"    "Rice"       "Biomass"    "Termites"   "Others"    
base_nation_error=c(0.379, 0.472, 0.062, 0.067, 0.056, 1.4, 0.558, 2, 0.5, 0.4)/2
all_figa=NULL
all_figb=NULL

iter=1
for(iter in 1:1000){
print(iter)	
nation_error= base_nation_error
for(k in 1:10){nation_error[k]=rnorm(1, mean= base_nation_error[k], sd=base_nation_error[k]*0.25)}
nation_error[nation_error<=0]= base_nation_error[nation_error<=0]

igrid=1
for(igrid in 1:1207){
emis=WW[,igrid]
f0 = ratio[igrid]
g0  = f0-1
alpha=emis/sum(emis)

gamma=sum(emis)^2*0.5^2/sum(nation_error^2*emis^2)
sigma=sqrt(gamma)* nation_error
beta=sum(alpha^2*sigma^2)/g0
g=alpha*sigma^2/beta
f=g+1
emis_by_sector[,igrid]=f
}

grid_ratio=array(1,c(length(lon),length(lat)))
grid_ratio2=array(1,c(length(lon),length(lat)))
for(k in 1:1199)grid_ratio[Clusters==k]=ratio[k]
for(k in 1:1199) grid_ratio2[Clusters==k]=emis_by_sector[1,k]

#============================================
all_loc=NULL
loc=c(-99.4,-98.2,25.4,26.1); all_loc=rbind(all_loc,loc)
loc=c(-98.4,-97.5,21.8,22.7); all_loc=rbind(all_loc,loc)
loc=c(-98.0,-96.9,20.3,21.2); all_loc=rbind(all_loc,loc)
loc=c(-96.5,-95.5,18,19); all_loc=rbind(all_loc,loc)
loc=c(-94.5,-91.9,17.5,18.6); all_loc=rbind(all_loc,loc)
loc=c(-93.0,-91.7,18.6,19.8); all_loc=rbind(all_loc,loc)
loc=c(-100+0.5,-99+0.3,19.2,19.8); all_loc=rbind(all_loc,loc)

locations=list()

emis=NULL
for(k in 1:7){
loc=all_loc[k,]	
ind1=(lon>=loc[1] & lon<=loc[2])
ind2=(lat>=loc[3] & lat<=loc[4])
ind=(Clusters[ind1,ind2])
ind[ind==0]=NA
ind=na.omit(as.numeric(ind))
locations[[k]]=ind
ap1=sum(WW[1,ind])
ap2=sum(WW[1,ind]*emis_by_sector[1,ind])
# ap2=sum(WW[1,ind]*ratio[ind,1])
emis=rbind(emis,c(ap1,ap2))
}

colnames(emis)=c("Prior","Posterior")
emis=t(emis)
colnames(emis)=c("R1","R2","R3","R4","R5","R6","R7")

Prior=WW[,1:1199]
Posterior=WW[,1:1199]*emis_by_sector[,1:1199]

ap=apply(Prior[7:10,],c(2),sum)#others
Prior2=rbind(Prior[1:6,],ap)
ap=apply(Posterior[7:10,],c(2),sum)
Posterior2=rbind(Posterior[1:6,],ap)
rownames(Prior2)=c("Oil/gas","Coal","Landfill","Livestock","Wastewater","Wetland","Others")
rownames(Posterior2)=c("Oil/gas","Coal","Landfill","Livestock","Wastewater","Wetland","Others")
figa_error=cal_error(1:1199)[c(1:6,10),]
figa_error[7,]=1.5*figa_error[7,]
figb_error=NULL
for(k in 1:7){
	figb_error = rbind(figb_error, cal_error(locations[[k]])[1,])
}

MAT=cbind(apply(Prior2,c(1),sum),apply(Posterior2,c(1),sum))
colnames(MAT)=c("Prior","Posterior")
MAT=t(MAT)

figa=MAT
figb=emis

all_figa=abind(all_figa, figa, along=3)
all_figb=abind(all_figb, figb, along=3)
}


pdf(file="/Users/lu/Documents/CH4_Mexico/Figures/Figure 4/Figure_4_uncertainty.pdf",width=6.5,height=5)

# dev.new(width=6.5, height=5)
par(mar=c(3,2,1,3))
par(mfrow=c(2,1))

par(mai=c(0.3, 0.6, 0.25, 0.2), mgp=c(1.2, 0.2, 0), tcl=-0.2, ps=10)
figa=apply(all_figa,c(1,2),mean)
barCenters=barplot(height= figa,main="", ylab="Emissions (Gg)", beside=TRUE,col=rep(c("#E69F00", "#56B4E9"),7),bty="n",ylim=c(0,2300),border=1)
means= figa[2,];standardErrors= apply(all_figa,c(1,2),sd)[2,]
arrows(barCenters[2,], means-standardErrors*2, barCenters[2,], means+standardErrors*2, lwd=1, angle=90, code=3,length=0.02)
box()
title("(a) Prior and Posterior emissions from different sectors",font.main=1,cex.main=1)

par(mai=c(0.3, 0.6, 0.25, 0.2), mgp=c(1.2, 0.2, 0), tcl=-0.2, ps=10)
figa=apply(all_figb,c(1,2),mean)
barCenters=barplot(height= figb,main="", ylab="Emissions (Gg)", beside=TRUE,col=rep(c("#E69F00", "#56B4E9"),7),bty="n",ylim=c(0,1300),border=1)
means= figb[2,];standardErrors= apply(all_figa,c(1,2),sd)[2,]
arrows(barCenters[2,], means-standardErrors*2, barCenters[2,], means+standardErrors*2, lwd=1, angle=90, code=3,length=0.02)

legend("topright", c("Prior","Posterior"), cex=1.0, bty="n",fill=c("#E69F00","#56B4E9"),border=1)
# legend(x=18,y=1050, c("  Aircraft"), cex=1.0, bty="n",pch=16,text.col=2,border=1,col=2)
box()
title("(b) Prior and posterior emissions in major oil/gas basins",font.main=1,cex.main=1)

# Flaring_emis=c(10.2, 15.2, 20.3)*3
# arrows(18.25, 30.6, 18.25, 60.9, lwd=1, angle=90, code=3,length=0.02,col=2)

# Flaring_emis=c(18, 90)
# arrows(15.25, 18, 15.25, 90, lwd=1.5, angle=90, code=3,length=0.02,col=4)

dev.off()