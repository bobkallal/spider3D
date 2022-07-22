###kallal_fang_script
#R 3.6.1

library(Morpho)
library(rgl)
library(geomorph)
library(phytools)
library(abind)
library(Rvcg)
library(paleomorph)

setwd("/fang/pts")
ntax = 40
n.lm <- 43 #(3 points, 20 x 2 curve pts)
sp.dim <- 3

{
  template.lm<-read.table(file="/fang/fang_template/araneus_fang_template_pointsandcurvesonlyResampled.pts",skip=2,header=F,sep="") #template tps file of the template (landmark and curves but not surface points)
  template.lm<-template.lm[,2:4] #the coordinate data are in columns 2, 3, 4 (column 1 is the landmark name/label)
  template.mesh<-ply2mesh("/fang/ply/ascii/araneus.ply") #template surface (converted from binary to ascii)
  patch<-read.table(file="/fang/fang_template/araneus_fang_template_patchResampled.pts",skip=2,header=F,sep="") #template tps file with only surface points
  patch<-patch[,2:4] #the coordinate data are in columns 2, 3, 4 (column 1 is the landmark name/label)
  skmeshlist<-dir("/fang/ply/ascii",pattern=".ply") # list with link of mesh (link=directory?)
  }
  
# set working directory with .ply and .pts with all samples
ptslist<-dir(pattern='.pts')
ptsarray<-array(dim=c(n.lm,sp.dim,ntax)) #dim=c(number of landmarks and semilandmarks, number of dimensions, number of specimens)
for(i in 1:length(ptslist))
{
  ptsarray[,,i]<-as.matrix(read.table(file=ptslist[i],skip=2,header=F,sep="",row.names=1))
}


dimnames(ptsarray)[3]<-list(
  substr(dir("/fang/ply/ascii",pattern=".ply"),1,(nchar(dir("/fang/ply/ascii",pattern=".ply"))-4)))###gives specimens names based on file names in ply folder
arraylm<-ptsarray

dim(ptsarray)
dim(arraylm)
skmeshlist

{
  fixed <- as.integer(c(1:3))  # fixed landmarks
  SC1<-as.integer(c(1,4:23,2))  # outer curve
  SC2<-as.integer(c(1,24:43,2)) # inner curve
  }
  curve<-list(SC1,SC2)
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
checkLM(ptsarray,path="/fang/ply/ascii/",suffix=".ply",render="s",alpha=1,begin=1)

{
  atlas<-createAtlas(template.mesh,as.matrix(template.lm),as.matrix(patch), corrCurves=curvein, patchCurves=NULL, keep.fix=fixed)
  #save(atlas,file="./atlas.R")
  shade3d(atlas[["mesh"]],col=7)
  spheres3d(atlas[["landmarks"]][fixed,],col=1,radius=50)
  spheres3d(atlas[["landmarks"]][slidings,],col=3,radius=30)
  spheres3d(atlas[["patch"]],col=4,radius=20)
  #load("./atlas.R")
}

patchedtest100<-placePatch(atlas,ptsarray,path="/fang/ply/ascii",inflate=100)
#checkLM(patchedtest50, path="/fang/ply/ascii/", atlas=atlas)

#slide surface landmarks and curves
setwd("/fang/ply/ascii")
slidewithcurves<-slider3d(patchedtest100,SMvector = fixed, outlines = curvein, deselect = TRUE, surp = surface, sur.path="/fang/ply/ascii/", meshlist = skmeshlist, iterations=3,mc.cores=1)

specimens<-slidewithcurves$dataslide

#PCA borealis
fang.gpa<-align.procrustes(specimens)
fang.pca <- gm.prcomp(fang.gpa$gdf$coords)

#PCA geomorph
fang.gpa<-gpagen(specimens)
fang.pca <- gm.prcomp(fang.gpa$coords)

species <- c("alaranea","aotearoa","aphonopelma","araneus","arkys","australomimetus","cheiracanthium","chilenodes","dictyna","diplocephalusF","diplocephalusM","dolichognatha","dysdera","gradungulidae","hogna","holarchaea","kukulcania","latrodectus","liphistius","mysmenidae","nesticus","nicodamidae","novanapis","oecobius","opopaea","pahoroides","pararchaea","patu","phonognatha","pimoa","platyoides","plectreurys","salticus","scytodes","stegodyphus","synaphris","synotaxus","telemidae","theridiosoma","uloborus")
fang.gpa$gdf$species<-species
group <- c("ara","plp","myg","ara","ara","ara","rta","ara","rta","ara","ara","ara","syn","aus","rta","ara","fil","ara","mes","ara","ara","nic","ara","udo","syn","ara","ara","ara","ara","ara","rta","syn","rta","syn","ere","ara","ara","syn","ara","udo")
fang.gpa$gdf$group<-group

gg.shape.space(fang.pca,labels=fang.gpa$gdf$species)
color <- c("ara","plp","myg","ara","ara","ara","rta","ara","rta","ara","ara","ara","syn","aus","rta","ara","fil","ara","mes","ara","ara","nic","ara","udo","syn","ara","ara","ara","ara","ara","rta","syn","rta","syn","ere","ara","ara","syn","ara","udo")

abbreviations <- c("AL","AO","AP","AR","AK","AU","CR","CH","DC","DF","DM","DO","DY","GR","HG","HA","KU","LA","LI","MY","NE","NI","NO","OE","OP","PH","PA","PT","PG","PI","PS","PY","SA","SC","ST","SN","SX","TE","TH","UL")
fang.gpa$gdf$abbreviations<-abbreviations
gg.shape.space(fang.pca,group=fang.gpa$gdf$group,group.title='group',include.legend = TRUE,convex.hulls=TRUE,pt.size = 2, labels=abbreviations, color=c("black","red","purple","royalblue","brown","green","gray","orange","cyan","maroon","magenta"))


