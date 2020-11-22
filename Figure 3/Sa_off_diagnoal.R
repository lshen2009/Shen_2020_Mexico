rm(list=ls())
library(fields); library(maps); library(ncdf4);library(abind)
setwd("~/Documents")
source('Function/get_geo.R')
source('Function/get_met.R')

PN.cols = colorRampPalette(c("#36348F","#2D6AAF","#209ED6","#9ED8DE", "white","#FFE37C","#F59431","#EF662F","#EB3C2F"), space="rgb")

#========
filename="/Users/lu/Documents/CH4_Mexico_archive/Data/Emissions_new/emis_Flexgrid/create_new_mask/East_mexico_clusters_v2.nc"
datafile=nc_open(filename)
lon=ncvar_get(datafile,varid="lon")
lat=ncvar_get(datafile,varid="lat")
Clusters=ncvar_get(datafile,varid="Clusters")
nc_close(datafile)

setwd("~/Documents/CH4_Mexico_archive/CPU/TROPOMI_and_GC_ND49_new/Step4_inversion_pixel/offshorex10percent")
datafile=nc_open("inversion_result.nc")
all_part1 =ncvar_get(datafile,varid="all_part1")
all_part2 =ncvar_get(datafile,varid="all_part2")
ratio =ncvar_get(datafile,varid="ratio")
nc_close(datafile)


#=== center lon's and lat's ===
center_loc=NULL
for(iclust in 1:1207){
	loc=which(Clusters == iclust,arr.ind=T)
	if(iclust<=1199){
	  center_loc = rbind(center_loc, c(lon[loc[1]], lat[loc[2]]))	
	}else{
	  center_loc = rbind(center_loc, c(mean(lon[loc[,1]]), mean(lat[loc[,2]]) ))
	}
}

Correlations=array(NA,c(1207,1207))
for(i in 1:1207){
	for(j in 1:1207){
		lat1=center_loc[i,2]
		lon1=center_loc[i,1]
		lat2=center_loc[j,2]
		lon2=center_loc[j,1]
		if(i==j){
  		 Correlations[i,j]=1
		}else{
		 ap2=dist.latlon(lat1,lon1,lat2,lon2)	
		 Correlations[i,j]= 0.97101235*exp(-0.03212729*ap2)		 
		}
	}
}
# plot.field(Correlations[1:100,1:100],1:100,1:100)

emis_error=array(NA,c(1207,1207))
for(i in 1:1207){
	for(j in 1:1207){
		emis_error[i,j]= Correlations[i,j]*0.5*0.5
	}
}

# emis_error=rep(0.5^2,1207);emis_error[1200:1207]=0.5^2
inv_Sa =solve(emis_error)

gamma=0.1
ratio=solve(gamma*all_part1+inv_Sa)%*%(gamma*all_part2)+1
S_posterior= solve(gamma*all_part1+inv_Sa)
A=diag(1207)-S_posterior%*%inv_Sa
sum(diag(A))

grid_ratio=array(1,c(length(lon),length(lat)))
grid_A=array(NA,c(length(lon),length(lat)))
for(k in 1:1199)grid_ratio[Clusters==k]=ratio[k]
for(k in 1:1199)grid_A[Clusters==k]=A[k,k]

# min(grid_ratio)
grid_ratio[grid_ratio<=0]=0
#dev.new(width=7,height=2.8)

setwd("~/Documents/CH4_Mexico/Figures/Figure 3")
pdf(file="inverse_test.pdf",width=6.5,height=2.5)
par(mar=c(3,2,1,3))
par(mfrow=c(1,2))
par(mai=c(0.3, 0.6, 0.25, 0.2), mgp=c(1.4, 0.4, 0), tcl=-0.2, ps=13)

plot.field(grid_ratio,lon,lat,type="sign",zlim=c(-0.5,2.5),legend.mar=5,mai=c(0.1,0.2,0.5,0.2),ps=10)
title("(a) Scaling factor (posterior/prior)",font.main=1,cex.main=0.9,line=+0.3)

plot.field(grid_A,lon,lat,type="def",zlim=c(0,0.4),legend.mar=7,col= rwb.colors(32)[16:32],mai=c(0.1,0.2,0.5,0.2),ps=10)
title("(b) Degree of freedoms (DOFS)",font.main=1,cex.main=0.9,line=+0.3)
ap=round(sum(grid_A,na.rm=T),digits=1)
text(x=-100,y=12,paste0("DOFS = ", ap),cex=0.9)

mtext("Posterior scaling factors and DOFS in eastern Mexico",side=3,outer=T,line=-1.5)

dev.off()