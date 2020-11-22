rm(list=ls())
library(fields); library(maps); library(ncdf4);library(abind)
library(geosphere)
setwd("~/Documents")
source('Function/get_geo.R')
source('Function/get_met.R')
source('Function/read_met.R')
source('Function/read_method.R')

#-- for onshore region only 
ss2=load("/Users/lu/Documents/CH4_Mexico/Data/Aircraft/low_resolution/overlaping_fraction.Rdata")
# plot.field(grid_fraction,lon,lat)
sum(overlaping_fraction[1:1199]*WW[1,1:1199])

#===== define functions =====
filename="/Users/lu/Documents/CH4_Mexico_archive/Data/Emissions_new/emis_Flexgrid/create_new_mask/East_mexico_clusters_v2.nc"
datafile=nc_open(filename)
lon=ncvar_get(datafile,varid="lon")
lat=ncvar_get(datafile,varid="lat")
Clusters=ncvar_get(datafile,varid="Clusters")
nc_close(datafile)

dirname="GC_format_offshoreX13%"

dirnames=c("GC_format_orig","GC_format_offshoreX13%","GC_format_offshoreX13%_mexicox50%","GC_format_offshoreX13%_onshorex100%")
all_gammas=c(0.1,0.25,0.5)

all_figs=list()
igamma=2;idirname=2

for(idirname in 1:length(dirnames)){
for(igamma in 1:length(all_gammas))	{
	
dirname=dirnames[idirname]
setwd(paste0("/Users/lu/Documents/CH4_Mexico/Data/", dirname,"/inversion"))

datafile=nc_open("inversion_result.nc")
all_part1 =ncvar_get(datafile,varid="all_part1")
all_part2 =ncvar_get(datafile,varid="all_part2")
ratio =ncvar_get(datafile,varid="ratio")
nc_close(datafile)

emis_error=rep(0.5^2,1207)
inv_Sa =diag(1/emis_error)
inv_gamma = all_gammas[igamma]
ratio=solve(inv_gamma*all_part1+inv_Sa)%*%(inv_gamma*all_part2)+1
temp=(ratio-1)^2/emis_error
S_posterior= solve(inv_gamma*all_part1+inv_Sa)
A=diag(1207)-S_posterior%*%inv_Sa
sum(diag(A))

ss=load("../Flexgrid/WW.Rdata")## WW is the emissions by sector
emis_by_sector=array(NA,c(10,1207))

cal_error=function(loc, cal_AK=F){
	WW2=array(0,dim(WW))
	WW2[,loc]=WW[,loc]
	S_red= WW2%*%S_posterior%*%t(WW2)
	error=cbind(sqrt(diag(WW2%*%diag(emis_error)%*%t(WW2))),sqrt(diag(S_red)))
	AK=NULL
	if(cal_AK){
	 W_star=t(WW2)%*%solve(WW2%*%t(WW2))
	 A_red=WW2%*%A%*%W_star
	 AK=diag(A_red)
	}
	colnames(error)=c("prior","posterior")
	return(list(error=error,AK=AK))
}

# 1 "Oil_gas"    "Coal"       "Landfill"   "Livestock"  "Wastewater"
# 6 "Wetland"    "Rice"       "Biomass"    "Termites"   "Others"    
nation_error=c(0.402, 0.135, 0.510, 0.158, 0.205, 1.40, 0.17, 1.06, 0.5, 0.5)/2

igrid=1195
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

sum(WW[1,1:1199]*ratio[1:1199])
sum(WW[1,1:1199]*emis_by_sector[1,1:1199])
WW2=WW*emis_by_sector
cbind(apply(WW2[,1:1199],c(1),sum))

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

figb=NULL
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
figb=rbind(figb,c(ap1,ap2))
}

#============================================================================
#============================================================================
colnames(figb)=c("Prior","Posterior")
figb=t(figb)
colnames(figb)=c("R1","R2","R3","R4","R5","R6","R7")

Prior=WW[,1:1199]
Posterior=WW[,1:1199]*emis_by_sector[,1:1199]

avg_prior=apply(Prior,c(1),sum)
avg_posterior=apply(Posterior,c(1),sum)
met=cal_error(1:1199,cal_AK=T)
error=met[["error"]]

figa=array(NA,c(2,7))
figa[1,1:6]= avg_prior[1:6]
figa[2,1:6]= avg_posterior[1:6]
figa[1,7]=sum(avg_prior[7:10])
figa[2,7]=sum(avg_posterior[7:10])
figa_error=array(NA,c(7,2))
figa_error[1:6,]= error[1:6,]
figa_error[7,1]=sum(error[7:10,1])
figa_error[7,2]=sum(error[7:10,2])
colnames(figa)=c("Oil/gas","Coal","Landfill","Livestock","Wastewater","Wetland","Others")
rownames(figa)=c("Prior","Posterior")
AK=abind(met[["AK"]][1:6], sum(met[["AK"]][7:10]))

figb_error=NULL
for(k in 1:7){
	met=cal_error(locations[[k]])
	figb_error = rbind(figb_error, met[["error"]][1,])
}

#--- anthropogenic
# 1 "Oil_gas"    "Coal"       "Landfill"   "Livestock"  "Wastewater"
# 6 "Wetland"    "Rice"       "Biomass"    "Termites"   "Others" 
new_WW=rbind(apply(WW[c(1,2,3,4,5,7,8,10),],c(2),sum))
new_WW[,1200:1207]=0
S_red= new_WW%*%S_posterior%*%t(new_WW)
anthro_emis=c(sum(Posterior[c(1,2,3,4,5,7,8,10),]), sqrt(diag(S_red)))
sum(new_WW)

new_WW=rbind(apply(WW[c(6,9),],c(2),sum))
new_WW[,1200:1207]=0
S_red= new_WW%*%S_posterior%*%t(new_WW)
nature_emis=c(sum(Posterior[c(6,9),]),sqrt(diag(S_red)))
sum(new_WW)

#--- onshore ---
onshore=c(sum(overlaping_fraction[1:1199]*Prior[1,]),sum(overlaping_fraction[1:1199]*Posterior[1,]))
#------

name=paste0(dirname,"__", inv_gamma)
all_figs[[name]]=list(dirname=dirname, gamma= inv_gamma, figa=figa, figb=figb, figa_error= figa_error, figb_error= figb_error, figa_AK= AK, onshore= onshore, anthro_emis= anthro_emis, nature_emis= nature_emis)
}
}
save(all_figs,file="~/Documents/CH4_Mexico/Figures/Figure 4/multiple_choices/all_figs.Rdata")