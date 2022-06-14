############ prepare the phylogenetic and functional layers############

rm(list=ls())   ##### clear environemnt 

# packages
library(maptools)    
library(ggplot2)
library(ape)
library(sf)
library(classInt)
library(rgdal) 
library(letsR)
library(rasterVis)
library(viridis)
library(raster)
library(tmap)
library(dplyr)
library(tidyr)

TC_shps <- rgdal::readOGR(dsn =  "./input_shp", layer = "TCSum")
TC_shps

print(bbox(TC_shps), digits=12)
set_ll_warn(FALSE) 
set_ll_TOL(0.2) 
st_crs(TC_shps)

colnames(TC_shps@data) <- "binomial"


## using letsR to rasterize the layer with 1 degree resolution
PAMAhull_TCshps <- lets.presab(TC_shps, xmn = -180, xmx = 180, ymn = -90, ymx = 90,resol = 1, 
                             remove.cells = TRUE, remove.sp = TRUE, show.matrix = FALSE,count = TRUE,   
                             crs = CRS("+proj=longlat +datum=WGS84"))
summary(PAMAhull_TCshps)  ### summary the presence-absence matrix
plot(PAMAhull_TCshps, main = "TC Species Richness")


tree.pam_new <- PAMAhull_TCshps

rownames(tree.pam_new[[1]]) <- cellFromXY(object=tree.pam_new$Rich,xy=tree.pam_new$Pre[,1:2])
##Detect the cell of outlier in richness. Remove species in the outlier cell, but not in the eight neighbour cells
#The cell of outlier and its eight neighbour cells
id.out <- which(as.vector(tree.pam_new$Rich>3500))
row.nr <- rowFromCell(tree.pam_new$Rich,id.out)
col.nr <- colFromCell(tree.pam_new$Rich,id.out)
id.out.nb <- cellFromRowColCombine(tree.pam_new$Rich,(row.nr-1):(row.nr+1),(col.nr-1):(col.nr+1))
id.out.nb <- id.out.nb[!id.out.nb %in% id.out]

#Remove species in the outlier cell, but not in the neighbour 8 cells
id.pre.out <- which(rownames(tree.pam_new$Pre) %in% id.out)
id.pre.out.nb <- which(rownames(tree.pam_new$Pre) %in% id.out.nb)
occ.out.nb <- colSums(tree.pam_new$Pre[id.pre.out.nb,-c(1:2)])
occ.out <- tree.pam_new$Pre[id.pre.out,-c(1:2)]
id.spe.out <- which(occ.out.nb==0 & occ.out==1)
tree.pam_new$Presence_and_Absence_Matrix[id.pre.out,id.spe.out+2] <- 0

#Update the richnes of outlier cell
rich.out.corrected <- sum(tree.pam_new$Pre[id.pre.out,-c(1:2)])
tree.pam_new[[2]][id.out] <- rich.out.corrected

##check the corrected data
plot(tree.pam_new)

#### functional traits #######

trait <- read.csv("./TC_species_Traits.csv",header = T,sep=";",stringsAsFactors = F) %>% # Supplementary Materials dataset 3
  as_tibble()
trait %>% str()

#--store shapefile names as a list
TC_list <- list.files("./ahull_range", pattern = "\\.tif$", all.files=TRUE, full.names=T) ## tree species' range maps 

## keep species only have range maps

TC_traits <- left_join(TC_list, trait, by = "species") 

# a function that cut numeric trait into intervals
# @ df: a dataframe that contains numeric trait
# @ n: how many groups should be cut, by default, n=3
cut_trait <- function(df,n=NULL){
  n <- ifelse(is.null(n),3,n)
  df.n <- mutate_if(df,is.numeric, function(x){
    x.quantil <- cut(x,breaks = c(quantile(x,probs = seq(0,1,by=1/n))),
                     labels = paste0("Quantil_",1:n),include.lowest = T) 
  })
  return(df.n)
}

trait.cat <- cut_trait(TC_traits,20) # 20 levels for each trait

# numbers of observations for each trait in each groups
sapply(trait.cat[,-1],table)

# wide data into long data
trait.cat1 <- trait.cat %>% 
  gather(trait,value,2:9)


# add another column 
trait.cat1$tmp <- 1
trait.cat.f <- reshape2::dcast(trait.cat1,
                               species ~ trait+value,
                               value.var = "tmp") 

trait.cat.f <- mutate_if(trait.cat.f,is.numeric,replace_na,0) # replace na with 0
trait.cat.f %>% str()
dim(TC_traits)
dim(trait.cat.f)

############################################################
#### creat a raster layer for common use
head(tree.pam_new$Presence_and_Absence_Matrix[, c(1:2)])
## generate species richness raster
tree.general <- tree.pam_new$Richness_Raster
plot(tree.general)

### multiply the species richness matrix using each trait level

sp.df <- as.matrix(tree.pam_new$Presence_and_Absence_Matrix[, -c(1:2)]) 

sp.trait.pb.abun.ls <- lapply(2:ncol(trait.cat.f),function(x){
  x.trait.matrix <- as.matrix(trait.cat.f[,x,drop=F])
  x.trait.abun <- sp.df %*% x.trait.matrix
})

names(sp.trait.pb.abun.ls) <- colnames(trait.cat.f)[2:ncol(trait.cat.f)]
for (i in 1:160) {
 # message(sprintf("output: %s",[i]))
  tree.pam_new_leaf <- tree.general
  id <- as.numeric(rownames(tree.pam_new$Pre))
  tree.pam_new_leaf[id] <- sp.trait.pb.abun.ls[[i]][,1]
 # plot(tree.pam_new_leaf)
  tree.pam_new_leaf<-  mask(tree.pam_new_leaf, MAINL)
  tree.pam_new_leaf[tree.pam_new_leaf < 1] <- NA
  writeRaster(x = tree.pam_new_leaf, 
              filename = paste0("./traits/", names(sp.trait.pb.abun.ls)[i]),
              format="GTiff", overwrite = TRUE) ## save the corrected raster
}

### ends here

#######   phylogeny  #########
# get the genus-level phylogeny using the genus-level tree of the whole data

library(ape)
#### species-level phylogeny
tree_range <- read.tree("./tree_species_phylogeny.tre")  # Supplementary Materials dataset 2
### I will try to obtain the genus level data and to see how many species in each genus are
tips_genus <- read.csv("./tips_genus.csv", header = T, sep = ';')
summary(tips_genus)
tree_genus <- read.tree("./tree_species_phylogeny_genus.tre")  #4157 genus

#But get the tips to keep
keep.spp <-levels(summary_genus$Var1)

keeptips <- tree_genus$tip.label[match(keep.spp, tree_genus$tip.label)]

remove_taxa = setdiff(tree_genus$tip.label, keep.spp)
tree_less_genus_phy <- drop.tip(tree_genus,remove_taxa)

pvr_tree_less_genus <- PVRdecomp(tree_less_genus_phy, scale = TRUE)  ## PCOA analysis

## obtain the percentage of each PCOA axis
pvr_tree_less_genus.per <- round(pvr_tree_less_genus@Eigen$values/sum(pvr_tree_less_genus@Eigen$values)*100,1)  ## obtain the contribution of each axis
# check how many are bigger than one
sum(pvr_tree_less_genus.per >= 1) ## explain 40.6% of the total variations

pvr_tree_less_genus_eigen15 <- as.data.frame(pvr_tree_less_genus@Eigen$vectors[,1:15])  

pvr_tree_genus_eigen <- cbind(pvr_tree_less_genus_eigen15, summary_genus)

species_phylogeny <- as.data.frame(tree_range$tip.label)

## combine the species with pcoa 
names(pvr_tree_genus_eigen)[names(pvr_tree_genus_eigen) == "Var1"] <- "genus"
TC_phylogeny_pcoa <- left_join(species_phylogeny, pvr_tree_genus_eigen, by = "genus")
TC_phylogeny_pcoa <-TC_phylogeny_pcoa[, -c(2,18)]


# a function that cut numeric trait into intervals
# @ df: a dataframe that contains numeric trait
# @ n: how many groups should be cut, by default, n=3
cut_trait <- function(df,n=NULL){
  n <- ifelse(is.null(n),3,n)
  df.n <- mutate_if(df,is.numeric, function(x){
    x.quantil <- cut(x,breaks = c(quantile(x,probs = seq(0,1,by=1/n))),
                     labels = paste0("Quantil_",1:n),include.lowest = T) 
  })
  return(df.n)
}

phylogeny.cat <- cut_trait(TC_phylogeny_pcoa,20) # 20 levels for each

# numbers of observations for each trait in each groups
sapply(phylogeny.cat[,-1],table)
library(dplyr)
library(tidyr)
# wide data into long data
phylogeny.cat1 <- phylogeny.cat %>% 
  gather(trait,value,2:16)

# add another column 
phylogeny.cat1$tmp <- 1
phylogeny.cat.f <- reshape2::dcast(phylogeny.cat1,
                                   species ~ trait+value,
                                   value.var = "tmp") #long data into wide

phylogeny.cat.f  <- mutate_if(phylogeny.cat.f,is.numeric,replace_na,0) # replace na with 0

phylogeny.cat.f  %>% str()

dim(phylogeny.cat.f)

### then multiply the species richness matrix using each trait level

sp.df <- as.matrix(tree.pam_new$Presence_and_Absence_Matrix[, -c(1:2)]) 

sp.phy.pb.abun.ls <- lapply(2:ncol(phylogeny.cat.f),function(x){
  x.phy.matrix <- as.matrix(phylogeny.cat.f[,x,drop=F])
  x.phy.abun <- sp.df %*% x.phy.matrix
})

names(sp.phy.pb.abun.ls) <- colnames(phylogeny.cat.f)[2:ncol(phylogeny.cat.f)]

for (i in 1:300) {
  # message(sprintf("output: %s",[i]))
  tree.pam_new_phy <- tree.general
  id <- as.numeric(rownames(tree.pam_new$Pre))
  tree.pam_new_phy[id] <- sp.phy.pb.abun.ls[[i]]
  # plot(tree.pam_new_leaf)
  tree.pam_new_phy<-  mask(tree.pam_new_phy, MAINL)
  tree.pam_new_phy[tree.pam_new_phy < 1] <- NA
  writeRaster(x = tree.pam_new_phy, 
              filename = paste0("./phylogeny_tif/", names(sp.phy.pb.abun.ls)[i]),
              format="GTiff", overwrite = TRUE) ## save the corrected raster
}

### ends here

## then import these rasters into zonation

################## prepare the three NGO frameworks #####################
##import the three basic layers
## G200 (global 200 ecoregions)
## biodiversity hotspots
## Last of the wild

G200 <- readOGR(dsn = "./global200ecoregions", layer = "g200_terr")
Ghotspot <- readOGR(dsn = "./hotspots_2016_1", layer = "hotspots_2016_1")
GLastW <- readOGR(dsn = "./wilderness_WGS", layer = "Pressure_free_lands_09_Proje")

## change the projection of the last of the world map
WGS84 <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
st_crs(GLastW)
GLastW1 <- spTransform(GLastW, WGS84)
                       
plot(G200)
plot(Ghotspot)
plot(GLastW1)

## rasterise the above three shp files

worldRaster <- raster(system.file("external/bioclim/current/bio3.grd", package = "biomod2"))
worldRaster[!is.na(worldRaster)] <- 0
worldRaster <- resample(worldRaster, trait_rank,method="bilinear")  ## original is 3 degree, now change to 1 degree
plot(worldRaster, axes = F, box = F, legend = F, main = "The world")

g200_raster <- shp2raster(shp = G200, mask.raster = worldRaster, label = "G200", transform = FALSE, value = 1)
plot(g200_raster)
g200_raster

Ghotspot_r <- shp2raster(shp = Ghotspot, mask.raster = worldRaster, label = "Ghotspot", transform = FALSE, value = 1)
plot(Ghotspot_r)
Ghotspot_r

GLastW_r <- shp2raster(shp = GLastW1, mask.raster = worldRaster, label = "GLastW", transform = FALSE, value = 1)
plot(GLastW_r)
GLastW_r


##### downland WDPA data ##############
library(wdpar)
library(dplyr)
library(ggmap)
library(sf)
#download global data
global_raw_data <- wdpa_fetch("global", wait = F)
plot(global_raw_data)
# clean the downloaded data

WDPA_web_levIV <- global_raw_data %>% 
  filter(
    MARINE == 0,
    IUCN_CAT %in% c("Ia", "Ib", "II", "III", "IV"), 
    STATUS %in% c("Designated", "Inscribed", "Established"),
    REP_AREA > 0,
    DESIG_ENG != "UNESCO-MAB Biosphere Reserve"
  )
# erase overlaps
WDPA_web_levIV <- st_erase_overlaps(WDPA_web_levIV)


st_write(WDPA_web_levIV, "./WDPA_web_levIV.shp")