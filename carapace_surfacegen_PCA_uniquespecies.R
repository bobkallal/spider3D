###kallal_carapace_script
#R 3.6.1

library(Morpho)
library(rgl)
library(geomorph)
library(phytools)
library(abind)
library(Rvcg)
library(paleomorph)

setwd("/carapace/pts_onedip")
ntax = 39 

n.lm <- 62  # all point (2) and curve (60) semilandmarks
sp.dim <- 3  # dimensions (2D or 3D) spatial dimensionality of data (2-D or 3-D)

# importing the template (which is araneus marmoreus in my case)
{
  template.lm<-read.table(file="/carapace/template/araneus_carapace_resampled_v3b_template_pointscurvesonly.pts",skip=2,header=F,sep="") #template tps file of the template (landmark and curves but not surface points)
  template.lm<-template.lm[,2:4] #the coordinate data are in columns 2, 3, 4 (column 1 is the landmark name/label)
  template.mesh<-ply2mesh("/carapace/template/araneus_v3_ascii.ply") #template surface (converted from binary to ascii)
  patch<-read.table(file="/carapace/template/araneus_carapace_resampled_v3b_template_surfaceonly.pts",skip=2,header=F,sep="") #template tps file with only surface points
  patch<-patch[,2:4] #the coordinate data are in columns 2, 3, 4 (column 1 is the landmark name/label)
  skmeshlist<-dir("/carapace/ply_onedip",pattern=".ply") # list with link of mesh (link=directory?)
  }

# set working directory with .ply and .pts with all samples
ptslist<-dir(pattern='.pts')
ptsarray<-array(dim=c(n.lm,sp.dim,ntax)) #dim=c(number of landmarks and semilandmarks, number of dimensions, number of specimens)
for(i in 1:length(ptslist))
{
  ptsarray[,,i]<-as.matrix(read.table(file=ptslist[i],skip=2,header=F,sep="",row.names=1))
}

#maybe working...
dimnames(ptsarray)[3]<-list(
  substr(dir("/carapace/ply_onedip",pattern=".ply"),1,(nchar(dir("/carapace/ply_onedip",pattern=".ply"))-4)))###gives specimens names based on file names in ply folder
arraylm<-ptsarray

dim(ptsarray)
dim(arraylm)
skmeshlist

{
  fixed <- as.integer(c(1:2))  # fixed landmarks
  SC1<-as.integer(c(1,3:32,2))  # lateral curve
  SC2<-as.integer(c(1,33:62,2)) # median curve
  }
  
curve<-list(SC1,SC2) 

#list point of curves without landmark at each extremity
curvein<-curve
for(i in 1:length(curvein)){
  curvein[[i]]<-tail(head(curve[[i]],-1),-1)
}

# creation of a vector containing rowindices of sliding landmarks of curves
slidings<-c(tail(fixed,1):dim(arraylm)[1])
surface<-(c((dim(template.lm)[1]+1):(dim(template.lm)[1]+dim(patch)[1])))

length(surface)
dimnames(ptsarray)## this is specimen order
dim(ptsarray)

#check if pts and ply files are in sync
#checkLM(ptsarray,path="/carapace/ply/",suffix=".ply",render="s",alpha=1,begin=1)

#making the atlas/template
{
  atlas<-createAtlas(template.mesh,as.matrix(template.lm),as.matrix(patch), corrCurves=curvein, patchCurves=NULL, keep.fix=fixed)
  #save(atlas,file="./atlas.R")
  shade3d(atlas[["mesh"]],col=7)
  spheres3d(atlas[["landmarks"]][fixed,],col=1,radius=50)
  spheres3d(atlas[["landmarks"]][slidings,],col=3,radius=30)
  spheres3d(atlas[["patch"]],col=4,radius=20)
  #load("./atlas.R")
}

#see your atlas
plotAtlas(atlas)

#inflate = 250 best so far; 500 and 1000 too much
patchedtest250<-placePatch(atlas,ptsarray,path="/carapace/ply_onedip",inflate=250)
#checkLM(patchedtest250, path="/carapace/ply/", atlas=atlas)

#slide surface landmarks and curves
setwd("/carapace/ply_onedip")
slidewithcurves<-slider3d(patchedtest250,SMvector = fixed, outlines = curvein, deselect = TRUE, surp = surface, sur.path="/carapace/ply_onedip/", meshlist = skmeshlist, iterations=3,mc.cores=1)

specimens<-slidewithcurves$dataslide

#mirror left to right
add_col_or_row = function(x, n = 1, add_col = T, fill = 0)
{
    m1 = matrix(x, ncol = if(add_col) nrow(x) * ncol(x) else nrow(x), byrow = T)
    m2 = matrix(fill, nrow = if(add_col) dim(x)[3] else prod(dim(x)[-1]), 
                ncol = if(add_col) nrow(x) * n else n)
    array(t(cbind(m1, m2)), 
          c(nrow(x) + ((!add_col) * n), ncol(x) + (add_col * n), dim(x)[3]))
}

specimens2<-add_col_or_row(specimens,n=220,add_col=FALSE,fill=NA)
dimnames(specimens2)[3]<-dimnames(specimens)[3]
bilats<-cbind(c(3:32,63:252),c(253:472)) #define left and right pairs, creates dummy points 253-427 on right side
newarray<-mirrorfill(specimens2,l1=as.integer(c(1,2,33:62)),l2=bilats)
dimnames(newarray)[3]<-dimnames(specimens)[3]
car_onedip.gpa<-gpagen(newarray) #procrustes align
car_onedip.gpa$coords<-car_onedip.gpa$coords[-c(253:472),,]

#PCA and plot in morphospace
car.pca<-gm.prcomp(car.gpa$coords)
plot(car.pca,phylo=FALSE)

#extreme min and max of PC1
spheres3d(car.pca$shapes$shapes.comp1$min,radius=.001)
spheres3d(car.pca$shapes$shapes.comp1$max,radius=.001)

#borealis

devtools::install_github("aphanotus/borealis")
library(borealis)
gg.shape.space(pca)

tree<-read.tree("/carapace/bigtree.tre")


#input species, group names and abbreviations
species <- c("alaranea","aotearoa","aphonopelma","araneus","arkys","australomimetus","cheiracanthium","chilenodes","dictyna","diplocephalusF","diplocephalusM","dolichognatha","dysdera","gradungulidae","hogna","holarchaea","kukulcania","latrodectus","liphistius","mysmenidae","nesticus","nicodamidae","novanapis","oecobius","opopaea","pahoroides","pararchaea","patu","phonognatha","pimoa","platyoides","plectreurys","salticus","scytodes","stegodyphus","synaphris","synotaxus","telemidae","theridiosoma","uloborus")
car.gpa$gdf$species<-species
group <- c("ara","plp","myg","ara","ara","ara","rta","ara","rta","ara","ara","ara","syn","aus","rta","ara","fil","ara","mes","ara","ara","nic","ara","udo","syn","ara","ara","ara","ara","ara","rta","syn","rta","syn","ere","ara","ara","syn","ara","udo")
car.gpa$gdf$group<-group
abbreviations <- c("AL","AO","AP","AR","AK","AU","CR","CH","DC","DF","DM","DO","DY","GR","HG","HA","KU","LA","LI","MY","NE","NI","NO","OE","OP","PH","PA","PT","PG","PI","PS","PY","SA","SC","ST","SN","SX","TE","TH","UL")
car.gpa$gdf$abbreviations<-abbreviations

#morphospace plot
gg.shape.space(car.pca,group=car.gpa$gdf$group,group.title='group',include.legend = TRUE,convex.hulls=TRUE,pt.size = 2, labels=abbreviations, color=c("black","red","purple","royalblue","brown","green","gray","orange","cyan","maroon","magenta"))
