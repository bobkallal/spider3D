###kallal_paturon_script
#R 3.6.1

library(Morpho)
library(rgl)
library(geomorph)
library(phytools)
library(abind)
library(Rvcg)
library(paleomorph)
library(borealis)

setwd("/paturon/pts")
ntax = 40
n.lm <- 84 #(4 points, 20 x 4 curve pts)
sp.dim <- 3



{
  template.lm<-read.table(file="/paturon/template/araneus_paturon_template_pointscurves.pts",skip=2,header=F,sep="") #template tps file of the template (landmark and curves but not surface points)
  template.lm<-template.lm[,2:4] #the coordinate data are in columns 2, 3, 4 (column 1 is the landmark name/label)
  template.mesh<-ply2mesh("/paturon/ply/araneus.ply") #template surface (converted from binary to ascii)
  patch<-read.table(file="/paturon/template/araneus_paturon_template_surface.pts",skip=2,header=F,sep="") #template tps file with only surface points
  patch<-patch[,2:4] #the coordinate data are in columns 2, 3, 4 (column 1 is the landmark name/label)
  skmeshlist<-dir("/paturon/ply",pattern=".ply") # list with link of mesh (link=directory?)
  }
  
# set working directory with .ply and .pts with all samples
ptslist<-dir(pattern='.pts')
ptsarray<-array(dim=c(n.lm,sp.dim,ntax)) #dim=c(number of landmarks and semilandmarks, number of dimensions, number of specimens)
for(i in 1:length(ptslist))
{
  ptsarray[,,i]<-as.matrix(read.table(file=ptslist[i],skip=2,header=F,sep="",row.names=1))
}


dimnames(ptsarray)[3]<-list(
  substr(dir("/paturon/ply",pattern=".ply"),1,(nchar(dir("/paturon/ply",pattern=".ply"))-4)))###gives specimens names based on file names in ply folder
arraylm<-ptsarray

dim(ptsarray)
dim(arraylm)
skmeshlist

{
  fixed <- as.integer(c(1:4))  # fixed landmarks (good)
  SC1<-as.integer(c(1,5:24,2))  # top outer curve
  SC2<-as.integer(c(1,25:44,2)) # top inner curve
  SC3<-as.integer(c(3,45:64,4)) # bottom outer curve
  SC4<-as.integer(c(3,65:84,4)) # bottom inner curve
  }
  curve<-list(SC1,SC2,SC3,SC4)
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
#checkLM(ptsarray,path="/paturon/ply/",suffix=".ply",render="s",alpha=1,begin=1)

{
  atlas<-createAtlas(template.mesh,as.matrix(template.lm),as.matrix(patch), corrCurves=curvein, patchCurves=NULL, keep.fix=fixed)
  #save(atlas,file="./atlas.R")
  shade3d(atlas[["mesh"]],col=7)
  spheres3d(atlas[["landmarks"]][fixed,],col=1,radius=50)
  spheres3d(atlas[["landmarks"]][slidings,],col=3,radius=30)
  spheres3d(atlas[["patch"]],col=4,radius=20)
  #load("./atlas.R")
}

patchedtest100<-placePatch(atlas,ptsarray,path="/paturon/ply",inflate=100)
#checkLM(patchedtest100, path="/paturon/ply/", atlas=atlas)

#slide surface landmarks and curves
setwd("/paturon/ply")
slidewithcurves<-slider3d(patchedtest100,SMvector = fixed, outlines = curvein, deselect = TRUE, surp = surface, sur.path="/paturon/ply/", meshlist = skmeshlist, iterations=3,mc.cores=1)

specimens<-slidewithcurves$dataslide

#GPA (borealis)
pat.gpa<-align.procrustes(specimens)

#GPA (geomorph)
#pat.gpa2<-gpagen(specimens)

#PCA (geomorph & borealis)
pat.pca <- gm.prcomp(pat.gpa$gdf$coords)

#PCA (geomorph only)
#pat.pca2 <- gm.prcomp(pat.gpa2$coords)


species <- c("Alaranea","Aotearoa","Aphonopelma","Araneus","Arkys","Australomimetus","Cheiracanthium","Chilenodes","Dictyna","DiplocephalusF","DiplocephalusM","Dolichognatha","Dysdera","Gradungulidae","Hogna","Holarchaea","Kukulcania","Latrodectus","Liphistius","Mysmenidae","Nesticus","Nicodamidae","Novanapis","Oecobius","Opopaea","Pahoroides","Pararchaea","Patu","Phonognatha","Pimoa","Platyoides","Plectreurys","Salticus","Scytodes","Stegodyphus","Synaphris","Synotaxus","Telemidae","Theridiosoma","Uloborus")
pat.gpa$gdf$species<-species
group <- c("ara","plp","myg","ara","ara","ara","rta","ara","rta","ara","ara","ara","syn","aus","rta","ara","fil","ara","mes","ara","ara","nic","ara","udo","syn","ara","ara","ara","ara","ara","rta","syn","rta","syn","ere","ara","ara","syn","ara","udo")
pat.gpa$gdf$group<-group

gg.shape.space(pat.pca,labels=x.gpa$gdf$species)
color <- c("ara","plp","myg","ara","ara","ara","rta","ara","rta","ara","ara","ara","syn","aus","rta","ara","fil","ara","mes","ara","ara","nic","ara","udo","syn","ara","ara","ara","ara","ara","rta","syn","rta","syn","ere","ara","ara","syn","ara","udo")

 gg.shape.space(pat.pca,group=pat.gpa$gdf$group,group.title='group',include.legend = TRUE,convex.hulls=TRUE,pt.size = 2, labels=abbreviations, color=c("black","red","purple","royalblue","brown","green","gray","orange","cyan","maroon","magenta"))

