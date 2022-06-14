rm(list=ls())   ##### clear environemnt 


library(raster)
library(sf)
library(rgeos)
library(rgdal)
library(plyr)
library(ggplot2)
library(ggalluvial)
library(sf)
library(cowplot)
library(reshape2)
library(viridis)
library(ggpubr)
library(egg)

## check the proporation of each species located within PAs
GLWDPA_lev4_1deg_Raster <- raster("GLWDPA_lev4_1deg_Raster.tif")

setwd("./raster_shp")
inFiles <- list.files(pattern="\\.tif$", all.files=TRUE, full.names=T)
nFiles <-  length(inFiles)

inputRaster <- raster()
intersect_raster <- raster()   

for (iCtr in 1 : nFiles){  
  message(sprintf("intersection file: %s",inFiles[iCtr]))
  inputRaster <- raster(inFiles[iCtr])
  abc[iCtr] <- cellStats(inputRaster, "sum") 
  projection(inputRaster) <- projection(GLWDPA_lev4_1deg_Raster)
  intersect_raster <- raster::mask(inputRaster, GLWDPA_lev4_1deg_Raster)
  nbc[iCtr]  = cellStats(intersect_raster, "sum")
  
}   

original_cellNO <- as.data.frame(abc)
protected_cellNO <- as.data.frame(nbc)
TC_cellNO <- cbind(inFiles,original_cellNO, protected_cellNO)
TC_cellNO$ratio <- "0"
TC_cellNO$ratio <- TC_cellNO$nbc/TC_cellNO$abc
str(TC_cellNO)
hist(TC_cellNO$ratio)

####  do the same for the top 17 and 50% priority 
sr_rank_original <- raster("./outputs/zonation.CAZ_E.rank.compressed.tif") # import zonation result
MAINL <- readOGR(dsn = "./COUNTRIES", layer = "GSHHS_i_L1_simple") #import mainland

sr_rank <- raster::mask(sr_rank_original, MAINL)
## SR
sr_rank_top17 <- sr_rank
sr_rank_top17[sr_rank_top17 < 0.83] <- NA
sr_rank_top17[sr_rank_top17 >= 0.83] <- 1
plot(sr_rank_top17)

sr_rank_top50 <- sr_rank
sr_rank_top50[sr_rank_top50 < 0.50] <- NA
sr_rank_top50[sr_rank_top50 >= 0.50] <- 1
plot(sr_rank_top50)

### obtain the overlap between each species and the top 17% priority area

for (iCtr in 1 : nFiles){
  message(sprintf("intersection file: %s",inFiles[iCtr]))
  inputRaster <- raster(inFiles[iCtr])
  abc17[iCtr]  <- cellStats(inputRaster, "sum") 
  projection(inputRaster) <- projection(sr_rank_top17)
  intersect_raster <- raster::mask(inputRaster, sr_rank_top17)
  nbc17[iCtr]  = cellStats(intersect_raster, "sum") 
}   

original_cellNO <- as.data.frame(abc17)
priority17_cellNO <- as.data.frame(nbc17)

priority_TC_cellNO <- cbind(inFiles,original_cellNO, priority17_cellNO)
priority_TC_cellNO$ratio <- "0"
priority_TC_cellNO$ratio <- priority_TC_cellNO$nbc17/priority_TC_cellNO$abc
str(priority_TC_cellNO)
hist(priority_TC_cellNO$ratio)
## combine the original and priority 17 and plot the density plot
TC_cellNO$group <- "protected"
priority_TC_cellNO$group <- "Priority17"
names(priority_TC_cellNO)[3] <- "nbc"
names(priority_TC_cellNO)[2] <- "abc"
TC_priority_cellNO <- rbind(TC_cellNO, priority_TC_cellNO)

ggplot(TC_priority_cellNO, aes(x=ratio, color=as.factor(group))) +
  geom_density()
# Add mean lines
mu <- ddply(TC_priority_cellNO, "group", summarise, grp.mean=mean(ratio))
head(mu)

ggplot(TC_priority_cellNO, aes(x=ratio, color=group)) +
  geom_density()+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=group),
             linetype="dashed") +
  scale_color_brewer(palette="Dark2")

## then top 50%
for (iCtr in 1 : nFiles){
  message(sprintf("intersection file: %s",inFiles[iCtr]))
  inputRaster <- raster(inFiles[iCtr])
  # abc50[iCtr]  <- cellStats(inputRaster, "sum") 
  projection(inputRaster) <- projection(sr_rank_top50)
  intersect_raster <- raster::mask(inputRaster, sr_rank_top50)
  nbc50[iCtr]   = cellStats(intersect_raster, "sum") 
} 

priority50_cellNO <- as.data.frame(nbc50)
priority50_TC_cellNO <- cbind(inFiles, original_cellNO, priority50_cellNO)
priority50_TC_cellNO$ratio <- "0"
priority50_TC_cellNO$ratio <- priority50_TC_cellNO$nbc50/priority50_TC_cellNO$abc
str(priority50_TC_cellNO)
hist(priority50_TC_cellNO$ratio)
## combine the original and priority 17 and plot the density plot
#TC_cellNO$group <- "protected"
priority50_TC_cellNO$group <- "Priority50"
names(priority50_TC_cellNO)[3] <- "nbc"
names(priority50_TC_cellNO)[2] <- "abc"
TC_priority1750_cellNO <- rbind(TC_priority_cellNO, priority50_TC_cellNO)
str(TC_priority1750_cellNO)

mu1 <- ddply(TC_priority1750_cellNO, "group", summarise, grp.mean=mean(ratio))
head(mu1)


ggplot(TC_priority1750_cellNO, aes(x=ratio, color=group)) +
  geom_density()+
  geom_vline(data=mu1, aes(xintercept=grp.mean, color=group),
             linetype="dashed") +
  scale_color_brewer(palette="Dark2") +
  theme_bw()

### ranking with all three dimension at the same time  #####
# richness ranking 
rank_3d_original <- raster("./outputs/do_3_aspects_together.CAZ_E.rank.compressed.tif") # import zonation result
rank_3d_mask <- raster::mask(rank_3d_original, MAINL)
plot(rank_3d_original)
plot(rank_3d_mask)

### change the conservation maps into binary

rank_3d_top17 <- rank_3d_mask
rank_3d_top17[rank_3d_top17 < 0.83] <- NA
rank_3d_top17[rank_3d_top17 >= 0.83] <- 1
plot(rank_3d_top17)

rank_3d_top50 <- rank_3d_mask
rank_3d_top50[rank_3d_top50 < 0.50] <- NA
rank_3d_top50[rank_3d_top50 >= 0.50] <- 1
plot(rank_3d_top50)

## obtain the percentage of range within the top 17% priority areas for each species
for (iCtr in 1 : nFiles){
  message(sprintf("intersection file: %s",inFiles[iCtr]))
  inputRaster <- raster(inFiles[iCtr])
  projection(inputRaster) <- projection(rank_3d_top17)
  intersect_raster <- raster::mask(inputRaster, rank_3d_top17)
  nbc_3d_17[iCtr]  = cellStats(intersect_raster, "sum") 
}   

for (iCtr in 1 : nFiles){
  message(sprintf("intersection file: %s",inFiles[iCtr]))
  inputRaster <- raster(inFiles[iCtr])
  projection(inputRaster) <- projection(rank_3d_top50)
  intersect_raster <- raster::mask(inputRaster, rank_3d_top50)
  nbc_3d_50[iCtr]  = cellStats(intersect_raster, "sum") 
}

priority17_3d_cellNO <- as.data.frame(nbc_3d_17)

priority17_3d_TC_cellNO <- cbind(inFiles, original_cellNO, priority17_3d_cellNO)
priority17_3d_TC_cellNO$ratio <- "0"
priority17_3d_TC_cellNO$ratio <- priority17_3d_TC_cellNO$nbc_3d_17/priority17_3d_TC_cellNO$abc
str(priority17_3d_TC_cellNO)
hist(priority17_3d_TC_cellNO$ratio)
## combine the original and priority 17 and plot the density plot
TC_cellNO$group <- "protected"
priority17_3d_TC_cellNO$group <- "Priority17_3d"
names(priority17_3d_TC_cellNO)[3] <- "nbc"
names(priority17_3d_TC_cellNO)[2] <- "abc"
TC_priority17_3d_cellNO <- rbind(TC_priority1750_cellNO, priority17_3d_TC_cellNO)
str(TC_priority17_3d_cellNO)

mu2 <- ddply(TC_priority17_3d_cellNO, "group", summarise, grp.mean=mean(ratio))
head(mu2)

ggplot(TC_priority17_3d_cellNO, aes(x=ratio, color=group)) +
  geom_density()+
  geom_vline(data=mu2, aes(xintercept=grp.mean, color=group),
             linetype="dashed") +
  scale_color_brewer(palette="Dark2") +
  theme_bw()

nbc_3d_50

priority50_3d_cellNO <- as.data.frame(nbc_3d_50)

priority50_3d_TC_cellNO <- cbind(inFiles, original_cellNO, priority50_3d_cellNO)
priority50_3d_TC_cellNO$ratio <- "0"
priority50_3d_TC_cellNO$ratio <- priority50_3d_TC_cellNO$nbc_3d_50/priority50_3d_TC_cellNO$abc
str(priority50_3d_TC_cellNO)
hist(priority50_3d_TC_cellNO$ratio)
## combine the original and priority 17 and plot the density plot
#TC_cellNO$group <- "protected"
priority50_3d_TC_cellNO$group <- "Priority50_3d"
names(priority50_3d_TC_cellNO)[3] <- "nbc"
names(priority50_3d_TC_cellNO)[2] <- "abc"
TC_priority17_50_3d_cellNO <- rbind(TC_priority17_3d_cellNO, priority50_3d_TC_cellNO)
str(TC_priority17_50_3d_cellNO)

mu3 <- ddply(TC_priority17_50_3d_cellNO, "group", summarise, grp.mean=mean(ratio))
head(mu3)
ddply(TC_priority17_50_3d_cellNO, "group", summarise, grp.median=median(ratio))

### combine them all
TC_ratio <- cbind( TC_cellNO, priority_TC_cellNO, priority50_TC_cellNO, priority17_3d_TC_cellNO,priority50_3d_TC_cellNO)

str(TC_ratio)

names(TC_ratio)[4] <- "PA_ratio"
names(TC_ratio)[9] <- "top17_ratio"
names(TC_ratio)[14] <- "top50_ratio"
names(TC_ratio)[19] <- "top17_3D_ratio"
names(TC_ratio)[24] <- "top50_3D_ratio"

names(TC_ratio)[7] <- "abc2"
names(TC_ratio)[8] <- "nbc_top17"
names(TC_ratio)[10] <- "group_top17"
names(TC_ratio)[12] <- "abc3"
names(TC_ratio)[13] <- "nbc_top50"
names(TC_ratio)[15] <- "group_top50"
names(TC_ratio)[17] <- "abc4"
names(TC_ratio)[18] <- "nbc_top17_3D"
names(TC_ratio)[20] <- "group_top17_3D"
names(TC_ratio)[22] <- "abc5"
names(TC_ratio)[23] <- "nbc_top50_3D"
names(TC_ratio)[25] <- "group_top50_3D"

ggplot(TC_priority17_50_3d_cellNO) +
  geom_histogram(aes(ratio, fill= group),position="fill", binwidth = 0.1) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,1)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  # geom_histogram(binwidth = 0.05) +
  scale_fill_brewer(palette="Set1") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #strip.background = element_rect(colour="white", fill="white"),
        text = element_text(size = 16),
        legend.position = "right") 


## prepare the bin
library(dplyr)
tags <- c("[0]", "(0-0.1]","(0.1-0.2]", "(0.2-0.3]", "(0.3-0.4]", "(0.4-0.5]", 
          "(0.5-0.6]","(0.6-0.7]", "(0.7-0.8]","(0.8-0.9]", "(0.9-1)", "[1]")

TC_ratio_PA_ratio <- TC_ratio %>% select(PA_ratio) #pick the variable 
TC_ratio_PA_ratio_group <- as_tibble(TC_ratio_PA_ratio) %>% 
  mutate(tag = case_when(
    PA_ratio <= 0 ~ tags[1],
    PA_ratio > 0 & PA_ratio <= 0.1 ~ tags[2],
    PA_ratio > 0.1 & PA_ratio <= 0.2 ~ tags[3],
    PA_ratio > 0.2 & PA_ratio <= 0.3 ~ tags[4],
    PA_ratio > 0.3 & PA_ratio <= 0.4 ~ tags[5],
    PA_ratio > 0.4 & PA_ratio <= 0.5 ~ tags[6],
    PA_ratio > 0.5 & PA_ratio <= 0.6 ~ tags[7],
    PA_ratio > 0.6 & PA_ratio <= 0.7 ~ tags[8],
    PA_ratio > 0.7 & PA_ratio <= 0.8 ~ tags[9],
    PA_ratio > 0.8 & PA_ratio <= 0.9 ~ tags[10],
    PA_ratio > 0.9 & PA_ratio < 1 ~ tags[11],
    PA_ratio >= 1 ~ tags[12]
  ))
summary(TC_ratio_PA_ratio_group)


TC_ratio_PA_ratio_group$tag <- factor(TC_ratio_PA_ratio_group$tag,
                                      levels = tags,
                                      ordered = FALSE)
summary(TC_ratio_PA_ratio_group$tag)

TC_ratio_top17_ratio <- TC_ratio %>% select(top17_ratio) 
TC_ratio_top17_ratio_group <- as_tibble(TC_ratio_top17_ratio) %>% 
  mutate(tag = case_when(
    top17_ratio <= 0 ~ tags[1],
    top17_ratio > 0 & top17_ratio <= 0.1 ~ tags[2],
    top17_ratio > 0.1 & top17_ratio <= 0.2 ~ tags[3],
    top17_ratio > 0.2 & top17_ratio <= 0.3 ~ tags[4],
    top17_ratio > 0.3 & top17_ratio <= 0.4 ~ tags[5],
    top17_ratio > 0.4 & top17_ratio <= 0.5 ~ tags[6],
    top17_ratio > 0.5 & top17_ratio <= 0.6 ~ tags[7],
    top17_ratio > 0.6 & top17_ratio <= 0.7 ~ tags[8],
    top17_ratio > 0.7 & top17_ratio <= 0.8 ~ tags[9],
    top17_ratio > 0.8 & top17_ratio <= 0.9 ~ tags[10],
    top17_ratio > 0.9 & top17_ratio < 1 ~ tags[11],
    top17_ratio >= 1 ~ tags[12]
  ))
summary(TC_ratio_top17_ratio_group)


TC_ratio_top17_ratio_group$tag <- factor(TC_ratio_top17_ratio_group$tag,
                                         levels = tags,
                                         ordered = FALSE)
summary(TC_ratio_top17_ratio_group$tag)

## top 50
TC_ratio_top50_ratio <- TC_ratio %>% select(top50_ratio) 
TC_ratio_top50_ratio_group <- as_tibble(TC_ratio_top50_ratio) %>% 
  mutate(tag = case_when(
    top50_ratio <= 0 ~ tags[1],
    top50_ratio > 0 & top50_ratio <= 0.1 ~ tags[2],
    top50_ratio > 0.1 & top50_ratio <= 0.2 ~ tags[3],
    top50_ratio > 0.2 & top50_ratio <= 0.3 ~ tags[4],
    top50_ratio > 0.3 & top50_ratio <= 0.4 ~ tags[5],
    top50_ratio > 0.4 & top50_ratio <= 0.5 ~ tags[6],
    top50_ratio > 0.5 & top50_ratio <= 0.6 ~ tags[7],
    top50_ratio > 0.6 & top50_ratio <= 0.7 ~ tags[8],
    top50_ratio > 0.7 & top50_ratio <= 0.8 ~ tags[9],
    top50_ratio > 0.8 & top50_ratio <= 0.9 ~ tags[10],
    top50_ratio > 0.9 & top50_ratio < 1 ~ tags[11],
    top50_ratio >= 1 ~ tags[12]
  ))
summary(TC_ratio_top50_ratio_group)


TC_ratio_top50_ratio_group$tag <- factor(TC_ratio_top50_ratio_group$tag,
                                         levels = tags,
                                         ordered = FALSE)
summary(TC_ratio_top50_ratio_group$tag)

## top 17_3D
TC_ratio_top17_3D_ratio <- TC_ratio %>% select(top17_3D_ratio) 
TC_ratio_top17_3D_ratio_group <- as_tibble(TC_ratio_top17_3D_ratio) %>% 
  mutate(tag = case_when(
    top17_3D_ratio <= 0 ~ tags[1],
    top17_3D_ratio > 0 & top17_3D_ratio <= 0.1 ~ tags[2],
    top17_3D_ratio > 0.1 & top17_3D_ratio <= 0.2 ~ tags[3],
    top17_3D_ratio > 0.2 & top17_3D_ratio <= 0.3 ~ tags[4],
    top17_3D_ratio > 0.3 & top17_3D_ratio <= 0.4 ~ tags[5],
    top17_3D_ratio > 0.4 & top17_3D_ratio <= 0.5 ~ tags[6],
    top17_3D_ratio > 0.5 & top17_3D_ratio <= 0.6 ~ tags[7],
    top17_3D_ratio > 0.6 & top17_3D_ratio <= 0.7 ~ tags[8],
    top17_3D_ratio > 0.7 & top17_3D_ratio <= 0.8 ~ tags[9],
    top17_3D_ratio > 0.8 & top17_3D_ratio <= 0.9 ~ tags[10],
    top17_3D_ratio > 0.9 & top17_3D_ratio < 1 ~ tags[11],
    top17_3D_ratio >= 1 ~ tags[12]
  ))
summary(TC_ratio_top17_3D_ratio_group)


TC_ratio_top17_3D_ratio_group$tag <- factor(TC_ratio_top17_3D_ratio_group$tag,
                                            levels = tags,
                                            ordered = FALSE)
summary(TC_ratio_top17_3D_ratio_group$tag)

## top50_3D
TC_ratio_top50_3D_ratio <- TC_ratio %>% select(top50_3D_ratio) 
TC_ratio_top50_3D_ratio_group <- as_tibble(TC_ratio_top50_3D_ratio) %>% 
  mutate(tag = case_when(
    top50_3D_ratio <= 0 ~ tags[1],
    top50_3D_ratio > 0 & top50_3D_ratio <= 0.1 ~ tags[2],
    top50_3D_ratio > 0.1 & top50_3D_ratio <= 0.2 ~ tags[3],
    top50_3D_ratio > 0.2 & top50_3D_ratio <= 0.3 ~ tags[4],
    top50_3D_ratio > 0.3 & top50_3D_ratio <= 0.4 ~ tags[5],
    top50_3D_ratio > 0.4 & top50_3D_ratio <= 0.5 ~ tags[6],
    top50_3D_ratio > 0.5 & top50_3D_ratio <= 0.6 ~ tags[7],
    top50_3D_ratio > 0.6 & top50_3D_ratio <= 0.7 ~ tags[8],
    top50_3D_ratio > 0.7 & top50_3D_ratio <= 0.8 ~ tags[9],
    top50_3D_ratio > 0.8 & top50_3D_ratio <= 0.9 ~ tags[10],
    top50_3D_ratio > 0.9 & top50_3D_ratio < 1 ~ tags[11],
    top50_3D_ratio >= 1 ~ tags[12]
  ))
summary(TC_ratio_top50_3D_ratio_group)


TC_ratio_top50_3D_ratio_group$tag <- factor(TC_ratio_top50_3D_ratio_group$tag,
                                            levels = tags,
                                            ordered = FALSE)
summary(TC_ratio_top50_3D_ratio_group$tag)

dat_pa <- as.data.frame(summary(TC_ratio_PA_ratio_group$tag))
dat_top17 <- as.data.frame(summary(TC_ratio_top17_ratio_group$tag))
dat_top50 <- as.data.frame(summary(TC_ratio_top50_ratio_group$tag))
dat_top17_3d <- as.data.frame(summary(TC_ratio_top17_3D_ratio_group$tag))
dat1 <- as.data.frame(summary(TC_ratio_top50_3D_ratio_group$tag))

dat_pa$group <- "PA"
dat_top17$group <- "Top 17%"
dat_top50$group <- "Top 50%"
dat_top17_3d$group <- "Top 30% (3D)"
dat1$group <- "Top 50% (3D)"

names(dat_pa)[1] <- "Fraction"
names(dat_top17)[1] <- "Fraction"
names(dat_top50)[1] <- "Fraction"
names(dat_top17_3d)[1] <- "Fraction"
names(dat1)[1] <- "Fraction"

dat_pa$tag <- rownames(dat_pa)
dat_top17$tag <- rownames(dat_top17)
dat_top50$tag <- rownames(dat_top50)
dat_top17_3d$tag <- rownames(dat_top17_3d)
dat1$tag <-rownames(dat1)

TC_ratio_summary <- rbind(dat_pa, dat_top17, dat_top50, dat_top17_3d, dat1)

TC_ratio_summary$percent <- TC_ratio_summary$Fraction/46752
write.csv(TC_ratio_summary, file = "D:/Wenyong/2_WDPA_3D/analysis/TC_gap_analysis_percentage.csv")

TC_ratio_summary$tag <- factor(TC_ratio_summary$tag, levels = c("[1]","(0.9-1)", "(0.8-0.9]", "(0.7-0.8]","(0.6-0.7]", 
                                                                "(0.5-0.6]","(0.4-0.5]", "(0.3-0.4]", "(0.2-0.3]",
                                                                "(0.1-0.2]","(0-0.1]", "[0]"), ordered = TRUE)

library(viridis)

p_percent <- ggplot(TC_ratio_summary, aes(fill=tag, y=Fraction, x= group)) + 
  geom_bar(position="fill", stat="identity") +
   scale_fill_viridis(discrete=TRUE, option = "D") +
    theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 16),
        legend.position = "right") + 
  coord_flip() +
  ylab("Percentage (%)") +
  xlab("Groups")


p2 <- ggplot(TC_priority17_50_3d_cellNO, aes(x=log10(abc), y=ratio, color=group)) + 
  geom_point(size=0.5) +
  geom_smooth(span = 0.3) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,3.5)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  scale_color_brewer(palette="Set1") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 16),
        legend.position = "right") + 
  ylab("Fraction of distribution covered (%)") +
  xlab("log10(species range size)")

## plot as one 
library(cowplot)
pdf("./analysis/Fig.gap_analysis_species.pdf", useDingbats=FALSE, width=12, height=12)
plot_grid(p_percent, p2,  nrow = 2,rel_heights =  c(1, 2))

dev.off()

#### phylogenetic and functional dimensions

#### then, I will do the same for the top 17 and 50% priority 
phylogenetic_rank_original <- raster("./analysis/do_TC_phylogeny.CAZ_E.rank.compressed.tif")  # import zonation result for phylogentic diversity

phy_crop <- raster::crop(phylogenetic_rank_original, extent(sr_rank))

pd_rank_top17 <- phy_crop
pd_rank_top17[pd_rank_top17 < 0.83] <- NA
pd_rank_top17[pd_rank_top17 >= 0.83] <- 1
plot(pd_rank_top17)

pd_rank_top50 <- phy_crop
pd_rank_top50[pd_rank_top50 < 0.50] <- NA
pd_rank_top50[pd_rank_top50 >= 0.50] <- 1
plot(pd_rank_top50)



### obtain the overlap between each species and the top 17% priority area
for (iCtr in 1 : nFiles){
  message(sprintf("intersection file: %s",inFiles[iCtr]))
  inputRaster <- raster(inFiles[iCtr])
  projection(inputRaster) <- projection(pd_rank_top17)
  intersect_raster <- raster::mask(inputRaster, pd_rank_top17)
  nbc_pd17[iCtr]  = cellStats(intersect_raster, "sum") 
}   

priority17_pd_cellNO <- as.data.frame(nbc_pd17)

priority_TC_17pd_cellNO <- cbind(inFiles,original_cellNO, priority17_pd_cellNO)
priority_TC_17pd_cellNO$ratio <- "0"
priority_TC_17pd_cellNO$ratio <- priority_TC_17pd_cellNO$nbc_pd17/priority_TC_17pd_cellNO$abc
str(priority_TC_17pd_cellNO)
hist(priority_TC_17pd_cellNO$ratio)
## combine the original and priority 17 and plot the density plot
TC_cellNO$group <- "protected"
priority_TC_17pd_cellNO$group <- "Priority17_pd"
names(priority_TC_17pd_cellNO)[3] <- "nbc"
names(priority_TC_17pd_cellNO)[2] <- "abc"
priority_TC_PA_17pd_cellNO <- rbind(TC_cellNO, priority_TC_17pd_cellNO)

# Add mean lines
mu_priority_TC_17pd <- ddply(priority_TC_PA_17pd_cellNO, "group", summarise, grp.mean=mean(ratio))
head(mu_priority_TC_17pd)

## then top 50%
for (iCtr in 1 : nFiles){
  message(sprintf("intersection file: %s",inFiles[iCtr]))
  inputRaster <- raster(inFiles[iCtr])
  projection(inputRaster) <- projection(pd_rank_top50)
  intersect_raster <- raster::mask(inputRaster, pd_rank_top50)
  nbc_pd50[iCtr]   = cellStats(intersect_raster, "sum") 
} 

priority50_pd_cellNO <- as.data.frame(nbc_pd50)

priority50_pd_TC_cellNO <- cbind(inFiles, original_cellNO, priority50_pd_cellNO)
priority50_pd_TC_cellNO$ratio <- "0"
priority50_pd_TC_cellNO$ratio <- priority50_pd_TC_cellNO$nbc_pd50/priority50_pd_TC_cellNO$abc
str(priority50_pd_TC_cellNO)
hist(priority50_pd_TC_cellNO$ratio)
## combine the original and priority 17 and plot the density plot
#TC_cellNO$group <- "protected"
priority50_pd_TC_cellNO$group <- "Priority50_pd"
names(priority50_pd_TC_cellNO)[3] <- "nbc"
names(priority50_pd_TC_cellNO)[2] <- "abc"
TC_priority1750_pd_cellNO <- rbind(priority_TC_17pd_cellNO, priority50_pd_TC_cellNO)
str(TC_priority1750_pd_cellNO)

mu_priority_TC_1750pd <- ddply(TC_priority1750_pd_cellNO, "group", summarise, grp.mean=mean(ratio))
head(mu_priority_TC_1750pd)


## functional diversity
functional_rank_original <- raster("./analysis/do_TC_20traits_10Wrap.CAZ_E.rank.compressed.tif")  ## import zonation result for functional diversity
trait_rank <- raster::mask(functional_rank_original, MAINL)
trait_crop <- raster::crop(trait_rank, extent(sr_rank))
## trait
fd_rank_top17 <- trait_crop
fd_rank_top17[fd_rank_top17 < 0.83] <- NA
fd_rank_top17[fd_rank_top17 >= 0.83] <- 1
plot(fd_rank_top17)

fd_rank_top50 <- trait_crop
fd_rank_top50[fd_rank_top50 < 0.50] <- NA
fd_rank_top50[fd_rank_top50 >= 0.50] <- 1
plot(fd_rank_top50)



### obtain the overlap between each species and the top 17% priority area
for (iCtr in 1 : nFiles){
  message(sprintf("intersection file: %s",inFiles[iCtr]))
  inputRaster <- raster(inFiles[iCtr])
  projection(inputRaster) <- projection(fd_rank_top17)
  intersect_raster <- raster::mask(inputRaster, fd_rank_top17)
  nbc_fd17[iCtr]  = cellStats(intersect_raster, "sum") 
}   

priority17_fd_cellNO <- as.data.frame(nbc_fd17)

priority_TC_17fd_cellNO <- cbind(inFiles,original_cellNO, priority17_fd_cellNO)
priority_TC_17fd_cellNO$ratio <- "0"
priority_TC_17fd_cellNO$ratio <- priority_TC_17fd_cellNO$nbc_fd17/priority_TC_17fd_cellNO$abc
str(priority_TC_17fd_cellNO)
hist(priority_TC_17fd_cellNO$ratio)
## combine the original and priority 17 and plot the density plot
TC_cellNO$group <- "protected"
priority_TC_17fd_cellNO$group <- "Priority17_fd"
names(priority_TC_17fd_cellNO)[3] <- "nbc"
names(priority_TC_17fd_cellNO)[2] <- "abc"
priority_TC_PA_17fd_cellNO <- rbind(TC_cellNO, priority_TC_17fd_cellNO)
mu_priority_TC_17fd <- ddply(priority_TC_PA_17fd_cellNO, "group", summarise, grp.mean=mean(ratio))
head(mu_priority_TC_17fd)

## then top 50%
for (iCtr in 1 : nFiles){
  message(sprintf("intersection file: %s",inFiles[iCtr]))
  inputRaster <- raster(inFiles[iCtr])
  # abc50[iCtr]  <- cellStats(inputRaster, "sum") 
  projection(inputRaster) <- projection(fd_rank_top50)
  intersect_raster <- raster::mask(inputRaster, fd_rank_top50)
  nbc_fd50[iCtr]   = cellStats(intersect_raster, "sum") 
} 

priority50_fd_cellNO <- as.data.frame(nbc_fd50)

priority50_fd_TC_cellNO <- cbind(inFiles, original_cellNO, priority50_fd_cellNO)
priority50_fd_TC_cellNO$ratio <- "0"
priority50_fd_TC_cellNO$ratio <- priority50_fd_TC_cellNO$nbc_fd50/priority50_fd_TC_cellNO$abc
str(priority50_fd_TC_cellNO)
hist(priority50_fd_TC_cellNO$ratio)
## combine the original and priority 17 and plot the density plot
priority50_fd_TC_cellNO$group <- "Priority50_fd"
names(priority50_fd_TC_cellNO)[3] <- "nbc"
names(priority50_fd_TC_cellNO)[2] <- "abc"
TC_priority1750_fd_cellNO <- rbind(priority_TC_17fd_cellNO, priority50_fd_TC_cellNO)
str(TC_priority1750_fd_cellNO)

mu_priority_TC_1750fd <- ddply(TC_priority1750_fd_cellNO, "group", summarise, grp.mean=mean(ratio))
head(mu_priority_TC_1750fd)


### combine them 
TC_ratio_3d_single <- cbind(TC_cellNO, priority_TC_cellNO, priority50_TC_cellNO, priority_TC_17pd_cellNO,
                            priority50_pd_TC_cellNO,
                            priority_TC_17fd_cellNO, priority50_fd_TC_cellNO)

str(TC_ratio_3d_single)


names(TC_ratio_3d_single)[4] <- "PA_ratio"
names(TC_ratio_3d_single)[9] <- "top17_sr_ratio"
names(TC_ratio_3d_single)[14] <- "top50_sr_ratio"
names(TC_ratio_3d_single)[19] <- "top17_pd_ratio"
names(TC_ratio_3d_single)[24] <- "top50_pd_ratio"
names(TC_ratio_3d_single)[29] <- "top17_fd_ratio"
names(TC_ratio_3d_single)[34] <- "top50_fd_ratio"
TC_ratio_3d_single <- TC_ratio_3d_single[,-c(6,7,11,12,16,17,21,22,26,27,31,32)]
# prepare the bin
library(dplyr)
tags <- c("[0]", "(0-0.1]","(0.1-0.2]", "(0.2-0.3]", "(0.3-0.4]", "(0.4-0.5]", 
          "(0.5-0.6]","(0.6-0.7]", "(0.7-0.8]","(0.8-0.9]", "(0.9-1)", "[1]")

TC_ratio_PA_ratio <- TC_ratio_3d_single %>% select(PA_ratio) #pick the variable 
TC_ratio_PA_ratio_group <- as_tibble(TC_ratio_PA_ratio) %>% 
  mutate(tag = case_when(
    PA_ratio <= 0 ~ tags[1],
    PA_ratio > 0 & PA_ratio <= 0.1 ~ tags[2],
    PA_ratio > 0.1 & PA_ratio <= 0.2 ~ tags[3],
    PA_ratio > 0.2 & PA_ratio <= 0.3 ~ tags[4],
    PA_ratio > 0.3 & PA_ratio <= 0.4 ~ tags[5],
    PA_ratio > 0.4 & PA_ratio <= 0.5 ~ tags[6],
    PA_ratio > 0.5 & PA_ratio <= 0.6 ~ tags[7],
    PA_ratio > 0.6 & PA_ratio <= 0.7 ~ tags[8],
    PA_ratio > 0.7 & PA_ratio <= 0.8 ~ tags[9],
    PA_ratio > 0.8 & PA_ratio <= 0.9 ~ tags[10],
    PA_ratio > 0.9 & PA_ratio < 1 ~ tags[11],
    PA_ratio >= 1 ~ tags[12]
  ))
summary(TC_ratio_PA_ratio_group)


TC_ratio_PA_ratio_group$tag <- factor(TC_ratio_PA_ratio_group$tag,
                                      levels = tags,
                                      ordered = FALSE)
summary(TC_ratio_PA_ratio_group$tag)

TC_ratio_top17_ratio <- TC_ratio_3d_single %>% select(top17_sr_ratio) #pick the variable 
TC_ratio_top17_ratio_group <- as_tibble(TC_ratio_top17_ratio) %>% 
  mutate(tag = case_when(
    top17_sr_ratio <= 0 ~ tags[1],
    top17_sr_ratio > 0 & top17_sr_ratio <= 0.1 ~ tags[2],
    top17_sr_ratio > 0.1 & top17_sr_ratio <= 0.2 ~ tags[3],
    top17_sr_ratio > 0.2 & top17_sr_ratio <= 0.3 ~ tags[4],
    top17_sr_ratio > 0.3 & top17_sr_ratio <= 0.4 ~ tags[5],
    top17_sr_ratio > 0.4 & top17_sr_ratio <= 0.5 ~ tags[6],
    top17_sr_ratio > 0.5 & top17_sr_ratio <= 0.6 ~ tags[7],
    top17_sr_ratio > 0.6 & top17_sr_ratio <= 0.7 ~ tags[8],
    top17_sr_ratio > 0.7 & top17_sr_ratio <= 0.8 ~ tags[9],
    top17_sr_ratio > 0.8 & top17_sr_ratio <= 0.9 ~ tags[10],
    top17_sr_ratio > 0.9 & top17_sr_ratio < 1 ~ tags[11],
    top17_sr_ratio >= 1 ~ tags[12]
  ))
summary(TC_ratio_top17_ratio_group)


TC_ratio_top17_ratio_group$tag <- factor(TC_ratio_top17_ratio_group$tag,
                                         levels = tags,
                                         ordered = FALSE)
summary(TC_ratio_top17_ratio_group$tag)

## top 50
TC_ratio_top50_ratio <- TC_ratio_3d_single %>% select(top50_sr_ratio) #pick the variable 
TC_ratio_top50_ratio_group <- as_tibble(TC_ratio_top50_ratio) %>% 
  mutate(tag = case_when(
    top50_sr_ratio <= 0 ~ tags[1],
    top50_sr_ratio > 0 & top50_sr_ratio <= 0.1 ~ tags[2],
    top50_sr_ratio > 0.1 & top50_sr_ratio <= 0.2 ~ tags[3],
    top50_sr_ratio > 0.2 & top50_sr_ratio <= 0.3 ~ tags[4],
    top50_sr_ratio > 0.3 & top50_sr_ratio <= 0.4 ~ tags[5],
    top50_sr_ratio > 0.4 & top50_sr_ratio <= 0.5 ~ tags[6],
    top50_sr_ratio > 0.5 & top50_sr_ratio <= 0.6 ~ tags[7],
    top50_sr_ratio > 0.6 & top50_sr_ratio <= 0.7 ~ tags[8],
    top50_sr_ratio > 0.7 & top50_sr_ratio <= 0.8 ~ tags[9],
    top50_sr_ratio > 0.8 & top50_sr_ratio <= 0.9 ~ tags[10],
    top50_sr_ratio > 0.9 & top50_sr_ratio < 1 ~ tags[11],
    top50_sr_ratio >= 1 ~ tags[12]
  ))
summary(TC_ratio_top50_ratio_group)


TC_ratio_top50_ratio_group$tag <- factor(TC_ratio_top50_ratio_group$tag,
                                         levels = tags,
                                         ordered = FALSE)
summary(TC_ratio_top50_ratio_group$tag)

## top 17_pd
TC_ratio_top17_pd_ratio <- TC_ratio_3d_single %>% select(top17_pd_ratio) #pick the variable 
TC_ratio_top17_pd_ratio_group <- as_tibble(TC_ratio_top17_pd_ratio) %>% 
  mutate(tag = case_when(
    top17_pd_ratio <= 0 ~ tags[1],
    top17_pd_ratio > 0 & top17_pd_ratio <= 0.1 ~ tags[2],
    top17_pd_ratio > 0.1 & top17_pd_ratio <= 0.2 ~ tags[3],
    top17_pd_ratio > 0.2 & top17_pd_ratio <= 0.3 ~ tags[4],
    top17_pd_ratio > 0.3 & top17_pd_ratio <= 0.4 ~ tags[5],
    top17_pd_ratio > 0.4 & top17_pd_ratio <= 0.5 ~ tags[6],
    top17_pd_ratio > 0.5 & top17_pd_ratio <= 0.6 ~ tags[7],
    top17_pd_ratio > 0.6 & top17_pd_ratio <= 0.7 ~ tags[8],
    top17_pd_ratio > 0.7 & top17_pd_ratio <= 0.8 ~ tags[9],
    top17_pd_ratio > 0.8 & top17_pd_ratio <= 0.9 ~ tags[10],
    top17_pd_ratio > 0.9 & top17_pd_ratio < 1 ~ tags[11],
    top17_pd_ratio >= 1 ~ tags[12]
  ))
summary(TC_ratio_top17_pd_ratio_group)


TC_ratio_top17_pd_ratio_group$tag <- factor(TC_ratio_top17_pd_ratio_group$tag,
                                            levels = tags,
                                            ordered = FALSE)
summary(TC_ratio_top17_pd_ratio_group$tag)

## top50_pd
TC_ratio_top50_pd_ratio <- TC_ratio_3d_single %>% select(top50_pd_ratio) #pick the variable 
TC_ratio_top50_pd_ratio_group <- as_tibble(TC_ratio_top50_pd_ratio) %>% 
  mutate(tag = case_when(
    top50_pd_ratio <= 0 ~ tags[1],
    top50_pd_ratio > 0 & top50_pd_ratio <= 0.1 ~ tags[2],
    top50_pd_ratio > 0.1 & top50_pd_ratio <= 0.2 ~ tags[3],
    top50_pd_ratio > 0.2 & top50_pd_ratio <= 0.3 ~ tags[4],
    top50_pd_ratio > 0.3 & top50_pd_ratio <= 0.4 ~ tags[5],
    top50_pd_ratio > 0.4 & top50_pd_ratio <= 0.5 ~ tags[6],
    top50_pd_ratio > 0.5 & top50_pd_ratio <= 0.6 ~ tags[7],
    top50_pd_ratio > 0.6 & top50_pd_ratio <= 0.7 ~ tags[8],
    top50_pd_ratio > 0.7 & top50_pd_ratio <= 0.8 ~ tags[9],
    top50_pd_ratio > 0.8 & top50_pd_ratio <= 0.9 ~ tags[10],
    top50_pd_ratio > 0.9 & top50_pd_ratio < 1 ~ tags[11],
    top50_pd_ratio >= 1 ~ tags[12]
  ))
summary(TC_ratio_top50_pd_ratio_group)


TC_ratio_top50_pd_ratio_group$tag <- factor(TC_ratio_top50_pd_ratio_group$tag,
                                            levels = tags,
                                            ordered = FALSE)
summary(TC_ratio_top50_pd_ratio_group$tag)


## top 17_fd
TC_ratio_top17_fd_ratio <- TC_ratio_3d_single %>% select(top17_fd_ratio) #pick the variable 
TC_ratio_top17_fd_ratio_group <- as_tibble(TC_ratio_top17_fd_ratio) %>% 
  mutate(tag = case_when(
    top17_fd_ratio <= 0 ~ tags[1],
    top17_fd_ratio > 0 & top17_fd_ratio <= 0.1 ~ tags[2],
    top17_fd_ratio > 0.1 & top17_fd_ratio <= 0.2 ~ tags[3],
    top17_fd_ratio > 0.2 & top17_fd_ratio <= 0.3 ~ tags[4],
    top17_fd_ratio > 0.3 & top17_fd_ratio <= 0.4 ~ tags[5],
    top17_fd_ratio > 0.4 & top17_fd_ratio <= 0.5 ~ tags[6],
    top17_fd_ratio > 0.5 & top17_fd_ratio <= 0.6 ~ tags[7],
    top17_fd_ratio > 0.6 & top17_fd_ratio <= 0.7 ~ tags[8],
    top17_fd_ratio > 0.7 & top17_fd_ratio <= 0.8 ~ tags[9],
    top17_fd_ratio > 0.8 & top17_fd_ratio <= 0.9 ~ tags[10],
    top17_fd_ratio > 0.9 & top17_fd_ratio < 1 ~ tags[11],
    top17_fd_ratio >= 1 ~ tags[12]
  ))
summary(TC_ratio_top17_fd_ratio_group)


TC_ratio_top17_fd_ratio_group$tag <- factor(TC_ratio_top17_fd_ratio_group$tag,
                                            levels = tags,
                                            ordered = FALSE)
summary(TC_ratio_top17_fd_ratio_group$tag)

## top50_fd
TC_ratio_top50_fd_ratio <- TC_ratio_3d_single %>% select(top50_fd_ratio) #pick the variable 
TC_ratio_top50_fd_ratio_group <- as_tibble(TC_ratio_top50_fd_ratio) %>% 
  mutate(tag = case_when(
    top50_fd_ratio <= 0 ~ tags[1],
    top50_fd_ratio > 0 & top50_fd_ratio <= 0.1 ~ tags[2],
    top50_fd_ratio > 0.1 & top50_fd_ratio <= 0.2 ~ tags[3],
    top50_fd_ratio > 0.2 & top50_fd_ratio <= 0.3 ~ tags[4],
    top50_fd_ratio > 0.3 & top50_fd_ratio <= 0.4 ~ tags[5],
    top50_fd_ratio > 0.4 & top50_fd_ratio <= 0.5 ~ tags[6],
    top50_fd_ratio > 0.5 & top50_fd_ratio <= 0.6 ~ tags[7],
    top50_fd_ratio > 0.6 & top50_fd_ratio <= 0.7 ~ tags[8],
    top50_fd_ratio > 0.7 & top50_fd_ratio <= 0.8 ~ tags[9],
    top50_fd_ratio > 0.8 & top50_fd_ratio <= 0.9 ~ tags[10],
    top50_fd_ratio > 0.9 & top50_fd_ratio < 1 ~ tags[11],
    top50_fd_ratio >= 1 ~ tags[12]
  ))
summary(TC_ratio_top50_fd_ratio_group)


TC_ratio_top50_fd_ratio_group$tag <- factor(TC_ratio_top50_fd_ratio_group$tag,
                                            levels = tags,
                                            ordered = FALSE)
summary(TC_ratio_top50_fd_ratio_group$tag)



dat_pa <- as.data.frame(summary(TC_ratio_PA_ratio_group$tag))
dat_top17 <- as.data.frame(summary(TC_ratio_top17_ratio_group$tag))
dat_top50 <- as.data.frame(summary(TC_ratio_top50_ratio_group$tag))
dat_top17_pd <- as.data.frame(summary(TC_ratio_top17_pd_ratio_group$tag))
dat_top50_pd <- as.data.frame(summary(TC_ratio_top50_pd_ratio_group$tag))
dat_top17_fd <- as.data.frame(summary(TC_ratio_top17_fd_ratio_group$tag))
dat_top50_fd <- as.data.frame(summary(TC_ratio_top50_fd_ratio_group$tag))


dat_pa$group <- "PA"
dat_top17$group <- "Top 17% (SR)"
dat_top50$group <- "Top 50% (SR)"
dat_top17_pd$group <- "Top 17% (PD)"
dat_top50_pd$group <- "Top 50% (PD)"
dat_top17_fd$group <- "Top 17% (FD)"
dat_top50_fd$group <- "Top 50% (FD)"

names(dat_pa)[1] <- "Fraction"
names(dat_top17)[1] <- "Fraction"
names(dat_top50)[1] <- "Fraction"
names(dat_top17_pd)[1] <- "Fraction"
names(dat_top50_pd)[1] <- "Fraction"
names(dat_top17_fd)[1] <- "Fraction"
names(dat_top50_fd)[1] <- "Fraction"

dat_pa$tag <- rownames(dat_pa)
dat_top17$tag <- rownames(dat_top17)
dat_top50$tag <- rownames(dat_top50)
dat_top17_pd$tag <- rownames(dat_top17_pd)
dat_top50_pd$tag <- rownames(dat_top50_pd)
dat_top17_fd$tag <- rownames(dat_top17_fd)
dat_top50_fd$tag <- rownames(dat_top50_fd)

TC_ratio_3_single_summary <- rbind(dat_pa, dat_top17, dat_top50, dat_top17_pd, dat_top50_pd,
                                   dat_top17_fd, dat_top50_fd)

TC_ratio_3_single_summary$percent <- TC_ratio_3_single_summary$Fraction/46752
write.csv(TC_ratio_3_single_summary, file = "D:/Wenyong/2_WDPA_3D/analysis/TC_gap_analysis_percentage.csv")

TC_ratio_3_single_summary$tag <- factor(TC_ratio_3_single_summary$tag, levels = c("[1]","(0.9-1)", "(0.8-0.9]", "(0.7-0.8]","(0.6-0.7]", 
                                                                                  "(0.5-0.6]","(0.4-0.5]", "(0.3-0.4]", "(0.2-0.3]",
                                                                                  "(0.1-0.2]","(0-0.1]", "[0]"), ordered = TRUE)


TC_ratio_3_single_summary$group <- factor(TC_ratio_3_single_summary$group, levels = c("PA","Top 17% (SR)", "Top 17% (PD)","Top 17% (FD)","Top 50% (SR)", "Top 50% (PD)", 
                                                                                      "Top 50% (FD)"), ordered = TRUE)


p_single_percent <- ggplot(TC_ratio_3_single_summary, aes(fill=tag, y=Fraction, x= group)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_viridis(discrete=TRUE, option = "D") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
          text = element_text(size = 16),
        legend.position = "right") + 
  ylab("Percentage (%)") +
  xlab("Groups")

## prepare the data for area and protection relationship
TC_ratio_3_singleDim <- rbind(TC_priority1750_cellNO, TC_priority1750_pd_cellNO, TC_priority1750_fd_cellNO)

TC_ratio_3_singleDim$group <- factor(TC_ratio_3_singleDim$group, levels = c("protected","Priority17", "Priority17_pd","Priority17_fd",
                                                                            "Priority50", "Priority50_pd", "Priority50_fd"), ordered = TRUE)

cbp <- c("#fb5607", "#4D9221", "#7FBC41", "#B8E186",   "#4575B4", "#74ADD1", "#ABD9E9")


p_area_protected_percent <-  ggplot(TC_ratio_3_singleDim, aes(x=log10(abc), y=ratio, color=group)) + 
  geom_point(size=0.5) +
  geom_smooth(span = 0.3) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,3.5)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  scale_color_manual(values = cbp) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 16),
        legend.position = "right") + 
  ylab("Percentage of range protected (%)") +
  xlab("log10(species range size)")
# scale_color_viridis(discrete=TRUE, option= "C")

## plot as one 

pdf("./analysis/Fig.gap_analysis_species_single_dimension.pdf", useDingbats=FALSE, width=12, height=12)
plot_grid(p_single_percent, p_area_protected_percent,  nrow = 2,rel_heights =  c(1, 2))

dev.off()

## prepare the four quantiles
str(TC_cellNO2)

TC_cellNO2 <- TC_cellNO %>% mutate(inFiles = gsub(".tif","", inFiles))

names(TC_cellNO2)[1] <- "files"

library(reshape2)
current_hmi_values_sp_all_wide <- dcast(current_hmi_values_sp_all_long, files ~ group, value.var = "value")
TC_cellNO2_hmi_df <- left_join(TC_cellNO2, current_hmi_values_sp_all_wide, by = "files")
TC_cellNO2_hmi_df[is.na(TC_cellNO2_hmi_df)] <- 0
library(Hmisc)
TC_cell2_him_df_cut <- table(cut2(TC_cellNO2_hmi_df$nbc, g = 4))  ## the break points are 2, 6, 17, i.e., 

## prepare Q1 data
TC_cellNO2_hmi_df_Q1 <- TC_cellNO2_hmi_df %>%  filter (nbc < 2) %>% droplevels()
TC_cellNO2_hmi_df_Q2 <- TC_cellNO2_hmi_df %>%  filter (nbc >= 2 & nbc < 6) %>% droplevels()
TC_cellNO2_hmi_df_Q3 <- TC_cellNO2_hmi_df %>%  filter (nbc >=6 & nbc <17) %>% droplevels()
TC_cellNO2_hmi_df_Q4 <- TC_cellNO2_hmi_df %>%  filter (nbc >=17) %>% droplevels()

TC_cellNO2_hmi_df_Q234 <- TC_cellNO2_hmi_df %>%  filter (nbc >=2 ) %>% droplevels()

##obtain the quantiles

TC_cellNO2_hmi_df %>% 
  summarise(mean  =mean(ratio))

TC_cellNO2_hmi_df %>% 
  summarise(mean  =mean(HMc))

TC_cellNO2_hmi_df %>% 
  summarise(mean  =mean(HMc_insidePA))

TC_cellNO2_hmi_df %>% 
  summarise(mean  =mean(HMc_outsidePA))

TC_cellNO2_hmi_df %>% 
  summarise(x  =quantile(ratio, c(0.25, 0.5, 0.75)), q = c(0.25, 0.5, 0.75))

TC_cellNO2_hmi_df %>% 
  summarise(x  =quantile(HMc, c(0.25, 0.5, 0.75)), q = c(0.25, 0.5, 0.75))

TC_cellNO2_hmi_df %>% 
  summarise(x  =quantile(HMc_insidePA, c(0.25, 0.5, 0.75)), q = c(0.25, 0.5, 0.75))

TC_cellNO2_hmi_df %>% 
  summarise(x  =quantile(HMc_outsidePA, c(0.25, 0.5, 0.75)), q = c(0.25, 0.5, 0.75))

count_global<- count(TC_cellNO2_hmi_df, vars = "ratio")

count(TC_cellNO2_hmi_df_Q1_clean, vars = "Protected")
count_Q124<- count(TC_cellNO2_hmi_df_Q234_clean, vars = "Protected")

write.csv(TC_cellNO2_hmi_df, file = "D:/Wenyong/2_WDPA_3D/analysis/TC_cellNO2_protected_ratio&hmi_df.csv")
## wide to long
names(TC_cellNO2_hmi_df_Q1)[4] <- "Protected"
TC_cellNO2_hmi_df_Q1_clean <- TC_cellNO2_hmi_df_Q1[, -c(2:3, 5)]
TC_cellNO2_hmi_df_Q1_long <- melt(TC_cellNO2_hmi_df_Q1_clean, id.vars = "files",
                                  variable.name="group")

names(TC_cellNO2_hmi_df_Q2)[4] <- "Protected"
TC_cellNO2_hmi_df_Q2_clean <- TC_cellNO2_hmi_df_Q2[, -c(2:3, 5)]
TC_cellNO2_hmi_df_Q2_long <- melt(TC_cellNO2_hmi_df_Q2_clean, id.vars = "files",
                                  variable.name="group")

names(TC_cellNO2_hmi_df_Q3)[4] <- "Protected"
TC_cellNO2_hmi_df_Q3_clean <- TC_cellNO2_hmi_df_Q3[, -c(2:3, 5)]
TC_cellNO2_hmi_df_Q3_long <- melt(TC_cellNO2_hmi_df_Q3_clean, id.vars = "files",
                                  variable.name="group")

names(TC_cellNO2_hmi_df_Q4)[4] <- "Protected"
TC_cellNO2_hmi_df_Q4_clean <- TC_cellNO2_hmi_df_Q4[, -c(2:3, 5)]
TC_cellNO2_hmi_df_Q4_long <- melt(TC_cellNO2_hmi_df_Q4_clean, id.vars = "files",
                                  variable.name="group")

names(TC_cellNO2_hmi_df_Q234)[4] <- "Protected"
TC_cellNO2_hmi_df_Q234_clean <- TC_cellNO2_hmi_df_Q234[, -c(2:3, 5)]
TC_cellNO2_hmi_df_Q234_long <- melt(TC_cellNO2_hmi_df_Q234_clean, id.vars = "files",
                                    variable.name="group")


hmi_q1_mean <- ddply(TC_cellNO2_hmi_df_Q1_long, "group", summarise, grp.mean=mean(value))
head(hmi_q1_mean)
hmi_q1_median <- ddply(TC_cellNO2_hmi_df_Q1_long, "group", summarise, grp.median=median(value))
head(hmi_q1_median)
##obtain the quantile for the first quantile data
TC_cellNO2_hmi_df_Q1_clean %>% 
   summarise(x  =quantile(Protected, c(0.25, 0.5, 0.75)), q = c(0.25, 0.5, 0.75))

TC_cellNO2_hmi_df_Q1_clean %>% 
  summarise(x  =quantile(HMc, c(0.25, 0.5, 0.75)), q = c(0.25, 0.5, 0.75))

TC_cellNO2_hmi_df_Q1_clean %>% 
  summarise(x  =quantile(HMc_insidePA, c(0.25, 0.5, 0.75)), q = c(0.25, 0.5, 0.75))

TC_cellNO2_hmi_df_Q1_clean %>% 
  summarise(x  =quantile(HMc_outsidePA, c(0.25, 0.5, 0.75)), q = c(0.25, 0.5, 0.75))


hmi_q2_mean <- ddply(TC_cellNO2_hmi_df_Q2_long, "group", summarise, grp.mean=mean(value))
head(hmi_q4_mean)
hmi_q3_mean <- ddply(TC_cellNO2_hmi_df_Q3_long, "group", summarise, grp.mean=mean(value))
head(hmi_q4_mean)
hmi_q4_mean <- ddply(TC_cellNO2_hmi_df_Q4_long, "group", summarise, grp.mean=mean(value))
head(hmi_q4_mean)
hmi_q234_mean <- ddply(TC_cellNO2_hmi_df_Q234_long, "group", summarise, grp.mean=mean(value))
head(hmi_q234_mean)

##obtain the quantile for the first quantile data
TC_cellNO2_hmi_df_Q234_clean %>% 
  summarise(x  =quantile(Protected, c(0.25, 0.5, 0.75)), q = c(0.25, 0.5, 0.75))

TC_cellNO2_hmi_df_Q234_clean %>% 
  summarise(x  =quantile(HMc, c(0.25, 0.5, 0.75)), q = c(0.25, 0.5, 0.75))

TC_cellNO2_hmi_df_Q234_clean %>% 
  summarise(x  =quantile(HMc_insidePA, c(0.25, 0.5, 0.75)), q = c(0.25, 0.5, 0.75))

TC_cellNO2_hmi_df_Q234_clean %>% 
  summarise(x  =quantile(HMc_outsidePA, c(0.25, 0.5, 0.75)), q = c(0.25, 0.5, 0.75))


Plot_Q1 <- ggplot(TC_cellNO2_hmi_df_Q1_long, aes(group, value, fill=group))  +
  geom_half_violin(side=2,  scale = "width") + 
  geom_boxplot(width=.1, position=position_nudge(x=-.2), outlier.shape = NA) +
  scale_fill_brewer(palette = "Spectral") +
  stat_summary(fun=mean, geom="crossbar", width = 0.09,size = 0.3, position=position_nudge(x=-.2),color="white") +
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin =0 , ymax =0.1,
           fill = "#6699CC",  alpha = 0.4) +
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin =0.1 , ymax =0.4,
           fill = "#999966",  alpha = 0.4) +
  scale_y_continuous(expand = c(0, 0),limits = c(0, 1))+
  theme_bw() +
  raincloud_theme + theme(legend.position = "none")


Plot_Q234 <- ggplot(TC_cellNO2_hmi_df_Q234_long, aes(group, value, fill=group))  +
  geom_half_violin(side=2,  scale = "width") + 
  geom_boxplot(width=.1, position=position_nudge(x=-.2), outlier.shape = NA) +
  scale_fill_brewer(palette = "Spectral") +
  stat_summary(fun=mean, geom="crossbar", width = 0.09,size = 0.3, position=position_nudge(x=-.2),color="white") +
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin =0 , ymax =0.1,
           fill = "#6699CC",  alpha = 0.4) +
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin =0.1 , ymax =0.4,
           fill = "#999966",  alpha = 0.4) +
  scale_y_continuous(expand = c(0, 0),limits = c(0, 1))+
  theme_bw() +
  raincloud_theme + theme(legend.position = "none")


pdf("./Fig.1.pdf", 
    useDingbats=FALSE, width=6, height=10)
plot_grid(p_global, Plot_Q1, Plot_Q234, ncol = 1,
          labels = "auto", label_size = 16, align = "hv" )
dev.off()


##### obtain the HMc for each species using the priority regions###########

# combine it with the protected precentage by priority areas
## protected precentage
TC_ratio_3d_single <- cbind(TC_cellNO, priority_TC_cellNO, priority50_TC_cellNO, priority_TC_17pd_cellNO,priority50_pd_TC_cellNO,
                            priority_TC_17fd_cellNO, priority50_fd_TC_cellNO)

str(TC_ratio_3d_single)
protected_percentage_prioity <- rbind(priority_TC_cellNO,priority50_TC_cellNO, priority_TC_17pd_cellNO,priority50_pd_TC_cellNO,
                                      priority_TC_17fd_cellNO, priority50_fd_TC_cellNO)
protected_percentage_prioity$group <- as.factor(protected_percentage_prioity$group)
summary(protected_percentage_prioity)
protected_percentage_prioity$group <- factor(protected_percentage_prioity$group, 
                                             levels = c("Priority17","Priority17_fd", "Priority17_pd", "Priority50",
                                                        "Priority50_fd","Priority50_pd"), ordered = TRUE)

p1_protected <- ggplot(protected_percentage_prioity, aes(group, ratio, fill=group))  +
  geom_half_violin(side=2,  scale = "area") + 
  geom_boxplot(width=.1, position=position_nudge(x=-.2), outlier.shape = NA) +
  stat_summary(fun=mean, geom="crossbar", width = 0.09,size = 0.3, position=position_nudge(x=-.2),color="white") +
  scale_fill_brewer(palette = "Spectral") +
  geom_hline(yintercept=0.4978920, linetype="dotted", color = "#2874C5", size=1)+ 
  coord_flip() + 
  theme_bw()+  raincloud_theme +theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      text = element_text(size = 16),
                                      legend.position = "none")
## inaide
hmi_inside_values_all <- rbind(hmi_inside_sr17_values_long_sp, hmi_inside_sr50_values_long_sp,
                               hmi_inside_pd17_values_long_sp, hmi_inside_pd50_values_long_sp,
                               hmi_inside_fd17_values_long_sp, hmi_inside_fd50_values_long_sp)
hmi_inside_values_all[is.na(hmi_inside_values_all)] <- 0
hmi_inside_values_all$group <- as.factor(hmi_inside_values_all$group)
summary(hmi_inside_values_all)

hmi_inside_values_all$group <- factor(hmi_inside_values_all$group, 
                                      levels = c("HMc_SR_in17", "HMc_FD_in17", "HMc_PD_in17","HMc_SR_in50","HMc_FD_in50",
                                                 "HMc_PD_in50"), ordered = TRUE)

p2_hmi_inside <- ggplot(hmi_inside_values_all, aes(group, value, fill=group))  +
  geom_half_violin(side=2, scale = "area") + 
  geom_boxplot(width=.1, position=position_nudge(x=-.2), outlier.shape = NA) +
  stat_summary(fun=mean, geom="crossbar",  width = 0.09,size = 0.3, position=position_nudge(x=-.2),color="white") +
  scale_fill_brewer(palette = "Spectral") +
  geom_hline(yintercept=0.115, linetype="dotted", color = "#2874C5", size=1)+ 
  coord_flip() +  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin =0 , ymax =0.1,
                           fill = "#6699CC",  alpha = 0.4) +
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin =0.1 , ymax =0.4,
           fill = "#999966",  alpha = 0.4) +
  scale_y_continuous(expand = c(0, 0),limits = c(0, 1))+
   theme_bw() +
  raincloud_theme +theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         axis.title.y=element_blank(),
                         axis.text.y=element_blank(),
                         text = element_text(size = 16),
                         legend.position = "none")

## outside

hmi_outside_values_all <- rbind(hmi_outside_sr17_values_long_sp, hmi_outside_sr50_values_long_sp,
                                hmi_outside_pd17_values_long_sp, hmi_outside_pd50_values_long_sp,
                                hmi_outside_fd17_values_long_sp, hmi_outside_fd50_values_long_sp)
hmi_outside_values_all[is.na(hmi_outside_values_all)] <- 0
hmi_outside_values_all$group <- as.factor(hmi_outside_values_all$group)
summary(hmi_outside_values_all)

hmi_outside_values_all$group <- factor(hmi_outside_values_all$group, 
                                       levels = c("HMc_SR_out17","HMc_FD_out17", "HMc_PD_out17",
                                                  "HMc_SR_out50", "HMc_FD_out50","HMc_PD_out50"), ordered = TRUE)

p3_hmi_outside <-ggplot(hmi_outside_values_all, aes(group, value, fill=group))  +
  geom_half_violin(side=2, scale = "area") + 
  geom_boxplot(width=.1, position=position_nudge(x=-.2), outlier.shape = NA) +
  stat_summary(fun=mean, geom="crossbar", width = 0.09,size = 0.3, position=position_nudge(x=-.2),color="white") +
  scale_fill_brewer(palette = "Spectral") +
  geom_hline(yintercept=0.257, linetype="dotted", color = "#2874C5", size=1)+  
  coord_flip() +  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin =0 , ymax =0.1,
                           fill = "#6699CC",  alpha = 0.4) +
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin =0.1 , ymax =0.4,
           fill = "#999966",  alpha = 0.4) +
  scale_y_continuous(expand = c(0, 0),limits = c(0, 1))+
  theme_bw() +
  raincloud_theme +theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         axis.title.y=element_blank(),
                         axis.text.y=element_blank(),
                         text = element_text(size = 16),
                         legend.position = "none")

#library(cowplot)
pdf("./analysis/Fig.pretected_pressure_priority_areas_13Oct2020.pdf", 
    useDingbats=FALSE, width=13, height=6)
plot_grid(p1_protected, p2_hmi_inside, p3_hmi_outside,  nrow = 1,
          labels = "auto", label_size = 16, align = "hv" )
dev.off()


raincloud_theme <- theme(
  text = element_text(size = 10),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text = element_text(size = 14),
  axis.text.x = element_text(angle = 0, vjust = 0.5),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 16),
  legend.position = c(0.9, 0.9),
  plot.title = element_text(lineheight = .8, face = "bold", size = 16),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
  axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"))

#### obtain the mean values for each group
mean_current <- ddply(TC_current_PA_pressure, "group", summarise, grp.mean=mean(value))
head(mean_current)

mean_priority_proportion <- ddply(protected_percentage_prioity, "group", summarise, grp.mean=mean(ratio))
head(mean_priority_proportion)
mean_hmi_priority_inside <- ddply(hmi_inside_values_all, "group", summarise, grp.mean=mean(value))
head(mean_hmi_priority_inside)
mean_hmi_priority_outside <- ddply(hmi_outside_values_all, "group", summarise, grp.mean=mean(value))
head(mean_hmi_priority_outside)



### use Sankey plot (alluvial diagram) to present the three (current, top 17% and top 50%) time/scanrios

TC_cellNO2 <- TC_cellNO %>% mutate(inFiles = gsub("./","", inFiles))
TC_cellNO2 <- TC_cellNO2 %>% mutate(inFiles = gsub(".tif","", inFiles))
priority_TC_cellNO1 <- priority_TC_cellNO %>% mutate(inFiles = gsub("./","", inFiles))
priority_TC_cellNO1 <- priority_TC_cellNO1 %>% mutate(inFiles = gsub(".tif","",inFiles))
priority50_TC_cellNO1 <- priority50_TC_cellNO %>% mutate(inFiles = gsub("./","", inFiles))
priority50_TC_cellNO1 <- priority50_TC_cellNO1 %>% mutate(inFiles = gsub(".tif","",inFiles))

TC_ratio_SR_single <- rbind(TC_cellNO2, priority_TC_cellNO1, priority50_TC_cellNO1)

TC_ratio_SR_single_ratio <- TC_ratio_SR_single %>% select(ratio) #pick the variable 
TC_ratio_SR_single_ratio_level <- as_tibble(TC_ratio_SR_single_ratio) %>% 
  mutate(tag = case_when(
    TC_ratio_SR_single_ratio <= 0 ~ tags[1],
    TC_ratio_SR_single_ratio > 0 & TC_ratio_SR_single_ratio <= 0.1 ~ tags[2],
    TC_ratio_SR_single_ratio > 0.1 & TC_ratio_SR_single_ratio <= 0.2 ~ tags[3],
    TC_ratio_SR_single_ratio > 0.2 & TC_ratio_SR_single_ratio <= 0.3 ~ tags[4],
    TC_ratio_SR_single_ratio > 0.3 & TC_ratio_SR_single_ratio <= 0.4 ~ tags[5],
    TC_ratio_SR_single_ratio > 0.4 & TC_ratio_SR_single_ratio <= 0.5 ~ tags[6],
    TC_ratio_SR_single_ratio > 0.5 & TC_ratio_SR_single_ratio <= 0.6 ~ tags[7],
    TC_ratio_SR_single_ratio > 0.6 & TC_ratio_SR_single_ratio <= 0.7 ~ tags[8],
    TC_ratio_SR_single_ratio > 0.7 & TC_ratio_SR_single_ratio <= 0.8 ~ tags[9],
    TC_ratio_SR_single_ratio > 0.8 & TC_ratio_SR_single_ratio <= 0.9 ~ tags[10],
    TC_ratio_SR_single_ratio > 0.9 & TC_ratio_SR_single_ratio < 1 ~ tags[11],
    TC_ratio_SR_single_ratio >= 1 ~ tags[12]
  ))
summary(TC_ratio_SR_single_ratio_level)


TC_ratio_SR_single_ratio_level$tag <- factor(TC_ratio_SR_single_ratio_level$tag,
                                             levels = tags,
                                             ordered = FALSE)
summary(TC_ratio_SR_single_ratio_level$tag)

TC_ratio_SR_single_final <- cbind(TC_ratio_SR_single, TC_ratio_SR_single_ratio_level[,2])


TC_ratio_SR_single_final$group <- factor(TC_ratio_SR_single_final$group, 
                                         levels = c("protected","Priority17", "Priority50"), ordered = TRUE)
TC_ratio_SR_single_final$tag <- factor(TC_ratio_SR_single_final$tag, 
                                       levels = c("[1]","(0.9-1)", "(0.8-0.9]", 
                                                  "(0.7-0.8]","(0.6-0.7]",
                                                  "(0.5-0.6]","(0.4-0.5]", "(0.3-0.4]", "(0.2-0.3]",
                                                  "(0.1-0.2]","(0-0.1]", "[0]"), ordered = TRUE)

p_protected_proportion <-ggplot(TC_ratio_SR_single_final,
                                aes(x = group, stratum = tag, alluvium = inFiles,
                                    fill = tag, label = tag)) +
  geom_flow() +
  geom_stratum(alpha = .9) +
  theme_bw() +
  raincloud_theme +
  scale_fill_viridis(discrete=TRUE, option = "D") +
  geom_text(stat = "stratum", size = 3) +
  theme(legend.position = "right") 

##  plot the same figure for HMI inside
current_hmi_values_sp_inside <- current_hmi_values_sp_all_long %>% filter(group == "HMc_insidePA") %>% droplevels()
hmi_inside_values_sr <- rbind(current_hmi_values_sp_inside, hmi_inside_sr17_values_long_sp, hmi_inside_sr50_values_long_sp)
hmi_inside_values_sr[is.na(hmi_inside_values_sr)] <- 0

summary(hmi_inside_values_sr)

hmi_inside_sr_values <- hmi_inside_values_sr %>% select(value) #pick the variable 
tags1 <- c("[0-0.1]","(0.1-0.4]", "(0.4-1]")

hmi_inside_sr_values_level <- as_tibble(hmi_inside_sr_values) %>% 
  mutate(tag = case_when(
    hmi_inside_sr_values >= 0 & hmi_inside_sr_values <= 0.1 ~ tags1[1],
    hmi_inside_sr_values > 0.1 & hmi_inside_sr_values <= 0.4 ~ tags1[2],
    hmi_inside_sr_values > 0.2 & hmi_inside_sr_values <= 1 ~ tags1[3]
  ))
summary(hmi_inside_sr_values_level)


hmi_inside_sr_values_level$tag <- factor(hmi_inside_sr_values_level$tag,
                                         levels = tags1,
                                         ordered = FALSE)
summary(hmi_inside_sr_values_level$tag)

hmi_inside_values_sr_final <- cbind(hmi_inside_values_sr, hmi_inside_sr_values_level[,2])


hmi_inside_values_sr_final$group <- factor(hmi_inside_values_sr_final$group, 
                                           levels = c("HMc_insidePA","HMc_SR_in17", "HMc_SR_in50"), ordered = TRUE)
hmi_inside_values_sr_final$tag <- factor(hmi_inside_values_sr_final$tag, 
                                         levels = c("(0.4-1]",  "(0.1-0.4]", "[0-0.1]"), ordered = TRUE)
library(ggalluvial)
p_hmi_inside <- ggplot(hmi_inside_values_sr_final,
                       aes(x = group, stratum = tag, alluvium = files,
                           fill = tag, label = tag)) +
  geom_flow() +
  geom_stratum(alpha = 0.8) +
  theme_bw() +
  raincloud_theme +
  scale_fill_viridis(discrete=TRUE, option = "C") +
  geom_text(stat = "stratum", size = 3) +
  theme(legend.position = "right") 


pdf("./analysis/Fig.SR pretected & pressure_priority_areas.pdf", 
    useDingbats=FALSE, width=20, height=8)

plot_grid(p_protected_proportion, p_hmi_inside,  nrow = 1,
          labels = "auto", label_size = 16, align = "hv" )
dev.off()


## obtain the summary of the presure and protected proportions

hmi_inside_values_sr_final

TC_ratio_SR_single_final

ddply(TC_ratio_SR_single_final, c("group", "tag"), summarise, freq = length(group))
ddply(hmi_inside_values_sr_final, c("group", "tag"), summarise, freq = length(group))
head(mean_priority_proportion)
## plot the sankey charts for the other two 
# priority_TC_17fd_cellNO, 
# priority50_fd_TC_cellNO
## FD
TC_cellNO2 <- TC_cellNO %>% mutate(inFiles = gsub("./","", inFiles))
TC_cellNO2 <- TC_cellNO2 %>% mutate(inFiles = gsub(".tif","", inFiles))  ## PAs
priority_TC_17fd_cellNO1 <- priority_TC_17fd_cellNO %>% mutate(inFiles = gsub("./","", inFiles))
priority_TC_17fd_cellNO1 <- priority_TC_17fd_cellNO1 %>% mutate(inFiles = gsub(".tif","",inFiles))
priority50_fd_TC_cellNO1 <- priority50_fd_TC_cellNO %>% mutate(inFiles = gsub("./","", inFiles))
priority50_fd_TC_cellNO1 <- priority50_fd_TC_cellNO1 %>% mutate(inFiles = gsub(".tif","",inFiles))

TC_ratio_fd_single <- rbind(TC_cellNO2, priority_TC_17fd_cellNO1, priority50_fd_TC_cellNO1)

TC_ratio_fd_single_ratio <- TC_ratio_fd_single %>% select(ratio) #pick the variable 
TC_ratio_fd_single_ratio_level <- as_tibble(TC_ratio_fd_single_ratio) %>% 
  mutate(tag = case_when(
    TC_ratio_fd_single_ratio <= 0 ~ tags[1],
    TC_ratio_fd_single_ratio > 0 & TC_ratio_fd_single_ratio <= 0.1 ~ tags[2],
    TC_ratio_fd_single_ratio > 0.1 & TC_ratio_fd_single_ratio <= 0.2 ~ tags[3],
    TC_ratio_fd_single_ratio > 0.2 & TC_ratio_fd_single_ratio <= 0.3 ~ tags[4],
    TC_ratio_fd_single_ratio > 0.3 & TC_ratio_fd_single_ratio <= 0.4 ~ tags[5],
    TC_ratio_fd_single_ratio > 0.4 & TC_ratio_fd_single_ratio <= 0.5 ~ tags[6],
    TC_ratio_fd_single_ratio > 0.5 & TC_ratio_fd_single_ratio <= 0.6 ~ tags[7],
    TC_ratio_fd_single_ratio > 0.6 & TC_ratio_fd_single_ratio <= 0.7 ~ tags[8],
    TC_ratio_fd_single_ratio > 0.7 & TC_ratio_fd_single_ratio <= 0.8 ~ tags[9],
    TC_ratio_fd_single_ratio > 0.8 & TC_ratio_fd_single_ratio <= 0.9 ~ tags[10],
    TC_ratio_fd_single_ratio > 0.9 & TC_ratio_fd_single_ratio < 1 ~ tags[11],
    TC_ratio_fd_single_ratio >= 1 ~ tags[12]
  ))
str(TC_ratio_fd_single_ratio_level)


TC_ratio_fd_single_ratio_level$tag <- factor(TC_ratio_fd_single_ratio_level$tag,
                                             levels = tags,
                                             ordered = FALSE)
summary(TC_ratio_fd_single_ratio_level$tag)

TC_ratio_fd_single_final <- cbind(TC_ratio_fd_single, TC_ratio_fd_single_ratio_level[,2])


TC_ratio_fd_single_final$group <- factor(TC_ratio_fd_single_final$group, 
                                         levels = c("protected","Priority17_fd", "Priority50_fd"), ordered = TRUE)
TC_ratio_fd_single_final$tag <- factor(TC_ratio_fd_single_final$tag, 
                                       levels = c("[1]","(0.9-1)", "(0.8-0.9]", 
                                                  "(0.7-0.8]","(0.6-0.7]",
                                                  "(0.5-0.6]","(0.4-0.5]", "(0.3-0.4]", "(0.2-0.3]",
                                                  "(0.1-0.2]","(0-0.1]", "[0]"), ordered = TRUE)
library(ggalluvial)
p_protected_fd_proportion <-ggplot(TC_ratio_fd_single_final,
                                   aes(x = group, stratum = tag, alluvium = inFiles,
                                       #weight = abc,
                                       fill = tag, label = tag)) +
  geom_flow() +
  geom_stratum(alpha = .9) +
  stat_summary(aes(y = ratio),fun = mean, 
               geom = "crossbar",
               width = 0.2, colour = "red", 
               position = position_dodge(width = .2)
  ) +
  theme_bw() +
  raincloud_theme +
  scale_fill_viridis(discrete=TRUE, option = "D") +
  geom_text(stat = "stratum", size = 3) +
  theme(legend.position = "right") 

### PD

# priority_TC_17pd_cellNO,
# priority50_pd_TC_cellNO,

priority_TC_17pd_cellNO1 <- priority_TC_17pd_cellNO %>% mutate(inFiles = gsub("./","", inFiles))
priority_TC_17pd_cellNO1 <- priority_TC_17pd_cellNO1 %>% mutate(inFiles = gsub(".tif","",inFiles))
priority50_pd_TC_cellNO1 <- priority50_pd_TC_cellNO %>% mutate(inFiles = gsub("./","", inFiles))
priority50_pd_TC_cellNO1 <- priority50_pd_TC_cellNO1 %>% mutate(inFiles = gsub(".tif","",inFiles))

TC_ratio_pd_single <- rbind(TC_cellNO2, priority_TC_17pd_cellNO1, priority50_pd_TC_cellNO1)

TC_ratio_pd_single_ratio <- TC_ratio_pd_single %>% select(ratio) #pick the variable 
TC_ratio_pd_single_ratio_level <- as_tibble(TC_ratio_pd_single_ratio) %>% 
  mutate(tag = case_when(
    TC_ratio_pd_single_ratio <= 0 ~ tags[1],
    TC_ratio_pd_single_ratio > 0 & TC_ratio_pd_single_ratio <= 0.1 ~ tags[2],
    TC_ratio_pd_single_ratio > 0.1 & TC_ratio_pd_single_ratio <= 0.2 ~ tags[3],
    TC_ratio_pd_single_ratio > 0.2 & TC_ratio_pd_single_ratio <= 0.3 ~ tags[4],
    TC_ratio_pd_single_ratio > 0.3 & TC_ratio_pd_single_ratio <= 0.4 ~ tags[5],
    TC_ratio_pd_single_ratio > 0.4 & TC_ratio_pd_single_ratio <= 0.5 ~ tags[6],
    TC_ratio_pd_single_ratio > 0.5 & TC_ratio_pd_single_ratio <= 0.6 ~ tags[7],
    TC_ratio_pd_single_ratio > 0.6 & TC_ratio_pd_single_ratio <= 0.7 ~ tags[8],
    TC_ratio_pd_single_ratio > 0.7 & TC_ratio_pd_single_ratio <= 0.8 ~ tags[9],
    TC_ratio_pd_single_ratio > 0.8 & TC_ratio_pd_single_ratio <= 0.9 ~ tags[10],
    TC_ratio_pd_single_ratio > 0.9 & TC_ratio_pd_single_ratio < 1 ~ tags[11],
    TC_ratio_pd_single_ratio >= 1 ~ tags[12]
  ))
str(TC_ratio_pd_single_ratio_level)


TC_ratio_pd_single_ratio_level$tag <- factor(TC_ratio_pd_single_ratio_level$tag,
                                             levels = tags,
                                             ordered = FALSE)
summary(TC_ratio_pd_single_ratio_level$tag)

TC_ratio_pd_single_final <- cbind(TC_ratio_pd_single, TC_ratio_pd_single_ratio_level[,2])


TC_ratio_pd_single_final$group <- factor(TC_ratio_pd_single_final$group, 
                                         levels = c("protected","Priority17_pd", "Priority50_pd"), ordered = TRUE)
TC_ratio_pd_single_final$tag <- factor(TC_ratio_pd_single_final$tag, 
                                       levels = c("[1]","(0.9-1)", "(0.8-0.9]", 
                                                  "(0.7-0.8]","(0.6-0.7]",
                                                  "(0.5-0.6]","(0.4-0.5]", "(0.3-0.4]", "(0.2-0.3]",
                                                  "(0.1-0.2]","(0-0.1]", "[0]"), ordered = TRUE)

p_protected_pd_proportion <-ggplot(TC_ratio_pd_single_final,
                                   aes(x = group, stratum = tag, alluvium = inFiles,
                                       #weight = abc,
                                       fill = tag, label = tag)) +
  geom_flow() +
  geom_stratum(alpha = .9) +
  theme_bw() +
  raincloud_theme +
  scale_fill_viridis(discrete=TRUE, option = "D") +
  geom_text(stat = "stratum", size = 3) +
  theme(legend.position = "right") 


##  plot the same figure for HMI inside
#fd
current_hmi_values_sp_inside <- current_hmi_values_sp_all_long %>% filter(group == "HMc_insidePA") %>% droplevels()
hmi_inside_values_fd <- rbind(current_hmi_values_sp_inside, hmi_inside_fd17_values_long_sp, hmi_inside_fd50_values_long_sp)
hmi_inside_values_fd[is.na(hmi_inside_values_fd)] <- 0

summary(hmi_inside_values_fd)

hmi_inside_fd_values <- hmi_inside_values_fd %>% select(value) #pick the variable 
tags1 <- c("[0-0.1]","(0.1-0.4]", "(0.4-1]")

hmi_inside_fd_values_level <- as_tibble(hmi_inside_fd_values) %>% 
  mutate(tag = case_when(
    hmi_inside_fd_values >= 0 & hmi_inside_fd_values <= 0.1 ~ tags1[1],
    hmi_inside_fd_values > 0.1 & hmi_inside_fd_values <= 0.4 ~ tags1[2],
    hmi_inside_fd_values > 0.2 & hmi_inside_fd_values <= 1 ~ tags1[3]
  ))
summary(hmi_inside_fd_values_level)


hmi_inside_fd_values_level$tag <- factor(hmi_inside_fd_values_level$tag,
                                         levels = tags1,
                                         ordered = FALSE)
summary(hmi_inside_fd_values_level$tag)

hmi_inside_values_fd_final <- cbind(hmi_inside_values_fd, hmi_inside_fd_values_level[,2])


hmi_inside_values_fd_final$group <- factor(hmi_inside_values_fd_final$group, 
                                           levels = c("HMc_insidePA","HMc_FD_in17", "HMc_FD_in50"), ordered = TRUE)
hmi_inside_values_fd_final$tag <- factor(hmi_inside_values_fd_final$tag, 
                                         levels = c("(0.4-1]",  "(0.1-0.4]", "[0-0.1]"), ordered = TRUE)
library(ggalluvial)
p_hmi_fd_inside <- ggplot(hmi_inside_values_fd_final,
                          aes(x = group, stratum = tag, alluvium = files,
                              #weight = abc,
                              fill = tag, label = tag)) +
  geom_flow() +
  geom_stratum(alpha = 0.8) +
  theme_bw() +
  raincloud_theme +
  scale_fill_viridis(discrete=TRUE, option = "C") +
  geom_text(stat = "stratum", size = 3) +
  theme(legend.position = "right") 


##PD
current_hmi_values_sp_inside <- current_hmi_values_sp_all_long %>% filter(group == "HMc_insidePA") %>% droplevels()
hmi_inside_values_pd <- rbind(current_hmi_values_sp_inside, hmi_inside_pd17_values_long_sp, hmi_inside_pd50_values_long_sp)
hmi_inside_values_pd[is.na(hmi_inside_values_pd)] <- 0

summary(hmi_inside_values_pd)

hmi_inside_pd_values <- hmi_inside_values_pd %>% select(value) #pick the variable 
tags1 <- c("[0-0.1]","(0.1-0.4]", "(0.4-1]")

hmi_inside_pd_values_level <- as_tibble(hmi_inside_pd_values) %>% 
  mutate(tag = case_when(
    hmi_inside_pd_values >= 0 & hmi_inside_pd_values <= 0.1 ~ tags1[1],
    hmi_inside_pd_values > 0.1 & hmi_inside_pd_values <= 0.4 ~ tags1[2],
    hmi_inside_pd_values > 0.2 & hmi_inside_pd_values <= 1 ~ tags1[3]
  ))
summary(hmi_inside_pd_values_level)


hmi_inside_pd_values_level$tag <- factor(hmi_inside_pd_values_level$tag,
                                         levels = tags1,
                                         ordered = FALSE)
summary(hmi_inside_pd_values_level$tag)

hmi_inside_values_pd_final <- cbind(hmi_inside_values_pd, hmi_inside_pd_values_level[,2])


hmi_inside_values_pd_final$group <- factor(hmi_inside_values_pd_final$group, 
                                           levels = c("HMc_insidePA","HMc_PD_in17", "HMc_PD_in50"), ordered = TRUE)
hmi_inside_values_pd_final$tag <- factor(hmi_inside_values_pd_final$tag, 
                                         levels = c("(0.4-1]",  "(0.1-0.4]", "[0-0.1]"), ordered = TRUE)
library(ggalluvial)
p_hmi_pd_inside <- ggplot(hmi_inside_values_pd_final,
                          aes(x = group, stratum = tag, alluvium = files,
                              #weight = abc,
                              fill = tag, label = tag)) +
  geom_flow() +
  geom_stratum(alpha = 0.8) +
  theme_bw() +
  raincloud_theme +
  scale_fill_viridis(discrete=TRUE, option = "C") +
  geom_text(stat = "stratum", size = 3) +
  theme(legend.position = "right") 


pdf("./analysis/Fig.FD_PD pretected & pressure_priority_areas.pdf", 
    useDingbats=FALSE, width=20, height=16)

plot_grid(p_protected_fd_proportion, p_hmi_fd_inside,  
          p_protected_pd_proportion, p_hmi_pd_inside,nrow = 2,
          labels = "auto", label_size = 16, align = "hv" )
dev.off()

## obtain the means
ddply(hmi_inside_values_sr_final, "group", summarise, mean = mean(value),median   = median(value))
ddply(hmi_inside_values_fd_final, "group", summarise, mean = mean(value),median   = median(value))
ddply(hmi_inside_values_pd_final, "group", summarise, mean = mean(value),median   = median(value))

ddply(TC_ratio_SR_single_final, "group", summarise, mean = mean(ratio),median   = median(ratio))
ddply(TC_ratio_fd_single_final, "group", summarise, mean = mean(ratio),median   = median(ratio))
ddply(TC_ratio_pd_single_final, "group", summarise, mean = mean(ratio),median   = median(ratio))


ddply(TC_ratio_SR_single_final, c("group", "tag"), summarise, freq = length(group))
ddply(hmi_inside_values_sr_final, c("group", "tag"), summarise, freq = length(group))

ddply(TC_ratio_fd_single_final, c("group", "tag"), summarise, freq = length(group))
ddply(hmi_inside_values_fd_final, c("group", "tag"), summarise, freq = length(group))

ddply(TC_ratio_pd_single_final, c("group", "tag"), summarise, freq = length(group))
ddply(hmi_inside_values_pd_final, c("group", "tag"), summarise, freq = length(group))


## plot the sankey plot using the jointing analysis, i.e. Fig. 3

TC_cellNO2 <- TC_cellNO %>% mutate(inFiles = gsub("./","", inFiles))
TC_cellNO2 <- TC_cellNO2 %>% mutate(inFiles = gsub(".tif","", inFiles))
priority17_3d_TC_cellNO1 <- priority17_3d_TC_cellNO %>% mutate(inFiles = gsub("./","", inFiles))
priority17_3d_TC_cellNO1 <- priority17_3d_TC_cellNO1 %>% mutate(inFiles = gsub(".tif","",inFiles))
priority50_3d_TC_cellNO1 <- priority50_3d_TC_cellNO %>% mutate(inFiles = gsub("./","", inFiles))
priority50_3d_TC_cellNO1 <- priority50_3d_TC_cellNO1 %>% mutate(inFiles = gsub(".tif","",inFiles))

TC_ratio_3d <- rbind(TC_cellNO2, priority17_3d_TC_cellNO1, priority50_3d_TC_cellNO1)

TC_ratio_3d_ratio <- TC_ratio_3d %>% select(ratio) #pick the variable 
TC_ratio_3d_ratio_level <- as_tibble(TC_ratio_3d_ratio) %>% 
  mutate(tag = case_when(
    TC_ratio_3d_ratio <= 0 ~ tags[1],
    TC_ratio_3d_ratio > 0 & TC_ratio_3d_ratio <= 0.1 ~ tags[2],
    TC_ratio_3d_ratio > 0.1 & TC_ratio_3d_ratio <= 0.2 ~ tags[3],
    TC_ratio_3d_ratio > 0.2 & TC_ratio_3d_ratio <= 0.3 ~ tags[4],
    TC_ratio_3d_ratio > 0.3 & TC_ratio_3d_ratio <= 0.4 ~ tags[5],
    TC_ratio_3d_ratio > 0.4 & TC_ratio_3d_ratio <= 0.5 ~ tags[6],
    TC_ratio_3d_ratio > 0.5 & TC_ratio_3d_ratio <= 0.6 ~ tags[7],
    TC_ratio_3d_ratio > 0.6 & TC_ratio_3d_ratio <= 0.7 ~ tags[8],
    TC_ratio_3d_ratio > 0.7 & TC_ratio_3d_ratio <= 0.8 ~ tags[9],
    TC_ratio_3d_ratio > 0.8 & TC_ratio_3d_ratio <= 0.9 ~ tags[10],
    TC_ratio_3d_ratio > 0.9 & TC_ratio_3d_ratio < 1 ~ tags[11],
    TC_ratio_3d_ratio >= 1 ~ tags[12]
  ))
summary(TC_ratio_3d_ratio_level)


TC_ratio_3d_ratio_level$tag <- factor(TC_ratio_3d_ratio_level$tag,
                                             levels = tags,
                                             ordered = FALSE)
summary(TC_ratio_3d_ratio_level$tag)

TC_ratio_3d_final <- cbind(TC_ratio_3d, TC_ratio_3d_ratio_level[,2])


TC_ratio_3d_final$group <- factor(TC_ratio_3d_final$group, 
                                         levels = c("protected","Priority17_3d", "Priority50_3d"), ordered = TRUE)
TC_ratio_3d_final$tag <- factor(TC_ratio_3d_final$tag, 
                                       levels = c("[1]","(0.9-1)", "(0.8-0.9]", 
                                                  "(0.7-0.8]","(0.6-0.7]",
                                                  "(0.5-0.6]","(0.4-0.5]", "(0.3-0.4]", "(0.2-0.3]",
                                                  "(0.1-0.2]","(0-0.1]", "[0]"), ordered = TRUE)
library(ggalluvial)
p_protected_3d_proportion <-ggplot(TC_ratio_3d_final,
                                aes(x = group, stratum = tag, alluvium = inFiles,
                                    #weight = abc,
                                    fill = tag, label = tag)) +
  geom_flow() +
  geom_stratum(alpha = .9) +
  theme_bw() +
  raincloud_theme +
  scale_fill_viridis(discrete=TRUE, option = "D") +
  geom_text(stat = "stratum", size = 3) +
  theme(legend.position = "right") 



gHM_original <- raster("./gHM/gHM.tif")
gHM_original ## 1000 m * 1000 m

## upscale the gHM to 5 km
gHM_5km <- aggregate(gHM_original, fact=5, fun=mean)

## change the projections

D3_rank_top17_high <- disaggregate(rank_3d_top17, 22)
D3_rank_top17_prj <- projectRaster(D3_rank_top17_high, crs=projection(gHM_5km))
D3_rank_top17_prj <- resample(D3_rank_top17_prj, gHM_5km)
plot(D3_rank_top17_prj)

D3_rank_top50_high <- disaggregate(rank_3d_top50, 22)
D3_rank_top50_prj <- projectRaster(D3_rank_top50_high, crs=projection(gHM_5km))
D3_rank_top50_prj <- resample(D3_rank_top50_prj, gHM_5km)
plot(D3_rank_top50_prj)



par(mfrow=c(1,2))
## prepare the inside and outside layers
gHM_3d_17_5km_inside <- mask(gHM_5km, D3_rank_top17_prj)
plot(gHM_3d_17_5km_inside)

gHM_3d_50_5km_inside <- mask(gHM_5km, D3_rank_top50_prj)
plot(gHM_3d_50_5km_inside)
TC_shp <- list.files(pattern = "\\.shp$", all.files=TRUE, full.names=T)
nFiles <-  length(TC_shp)
for (i in 1:nFiles
){
  print(i)
  inputshp <- readOGR(TC_shp[i])
  inputshp_prj <- spTransform(inputshp, projection(gHM_original)) 
  hmi_inside_3d_17_values[[i]] <- raster::extract(gHM_3d_17_5km_inside, inputshp_prj, method = "simple",
                                         fun = mean, na.rm = T,  df = F) 
}
## combine sp names and mean gHM values
hmi_inside_3d_17_values_long <- melt(hmi_inside_3d_17_values, id.vars = "") 
## combine sp names and mean gHM values
files<-sub(".shp","",TC_shp) files<-sub("./","",files)
hmi_inside_3d_17_values_long_sp <- cbind(files, hmi_inside_3d_17_values_long)
hist(hmi_inside_3d_17_values_long_sp$value)

### SR top 50% Inside

for (i in 1:nFiles
){
  print(i)
  inputshp <- readOGR(TC_shp[i])
  inputshp_prj <- spTransform(inputshp, projection(gHM_original)) 
  hmi_inside_3d_50_values[[i]] <- raster::extract(gHM_3d_50_5km_inside, inputshp_prj, method = "simple",
                                         fun = mean, na.rm = T,  df = F)
}
## combine sp names and mean gHM values
library(reshape2)
hmi_inside_3d_50_values_long <- melt(hmi_inside_3d_50_values, id.vars = "") 

## combine sp names and mean gHM values
hmi_inside_3d_50_values_long_sp <- cbind(files, hmi_inside_3d_50_values_long)
hist(hmi_inside_3d_50_values_long_sp$value)

##  plot the same figure for HMI inside
current_hmi_values_sp_inside <- current_hmi_values_sp_all_long %>% filter(group == "HMc_insidePA") %>% droplevels()
hmi_inside_3d_17_values_long_sp$group  <- "Priority 17% 3D"
hmi_inside_3d_50_values_long_sp$group  <- "Priority 50% 3D"
hmi_inside_values_3d <- rbind(current_hmi_values_sp_inside, hmi_inside_3d_17_values_long_sp, hmi_inside_3d_50_values_long_sp)
hmi_inside_values_3d[is.na(hmi_inside_values_3d)] <- 0

summary(hmi_inside_values_3d)

hmi_inside_3d_values <- hmi_inside_values_3d %>% select(value) #pick the variable 
tags1 <- c("[0-0.1]","(0.1-0.4]", "(0.4-1]")

hmi_inside_3d_values_level <- as_tibble(hmi_inside_3d_values) %>% 
  mutate(tag = case_when(
    hmi_inside_3d_values >= 0 & hmi_inside_3d_values <= 0.1 ~ tags1[1],
    hmi_inside_3d_values > 0.1 & hmi_inside_3d_values <= 0.4 ~ tags1[2],
    hmi_inside_3d_values > 0.2 & hmi_inside_3d_values <= 1 ~ tags1[3]
  ))
summary(hmi_inside_3d_values_level)


hmi_inside_3d_values_level$tag <- factor(hmi_inside_3d_values_level$tag,
                                         levels = tags1,
                                         ordered = FALSE)
summary(hmi_inside_3d_values_level$tag)

hmi_inside_values_3d_final <- cbind(hmi_inside_values_3d, hmi_inside_3d_values_level[,2])


hmi_inside_values_3d_final$group <- factor(hmi_inside_values_3d_final$group, 
                                           levels = c("HMc_insidePA","Priority 17% 3D", "Priority 50% 3D"), ordered = TRUE)
hmi_inside_values_3d_final$tag <- factor(hmi_inside_values_3d_final$tag, 
                                         levels = c("(0.4-1]",  "(0.1-0.4]", "[0-0.1]"), ordered = TRUE)
library(ggalluvial)
p_hmi_3d_inside <- ggplot(hmi_inside_values_3d_final,
                       aes(x = group, stratum = tag, alluvium = files,
                           #weight = abc,
                           fill = tag, label = tag)) +
  geom_flow() +
  geom_stratum(alpha = 0.8) +
  theme_bw() +
  raincloud_theme +
  scale_fill_viridis(discrete=TRUE, option = "C") +
  geom_text(stat = "stratum", size = 3) +
  theme(legend.position = "right") 

pdf("./analysis/Fig.3 pretected & pressure_priority_areas.pdf", 
    useDingbats=FALSE, width=20, height=8)

plot_grid(p_protected_3d_proportion, p_hmi_3d_inside,  nrow = 1,
          labels = "auto", label_size = 16, align = "hv" )
dev.off()


## obtain the summary of the presure and protected proportions

hmi_inside_values_3d_final

TC_ratio_SR_single_final

ddply(TC_ratio_3d_final, c("group", "tag"), summarise, freq = length(group))
ddply(hmi_inside_values_3d_final, c("group", "tag"), summarise, freq = length(group))

ddply(hmi_inside_values_3d_final, "group", summarise, mean = mean(value),median   = median(value))

ddply(TC_ratio_3d_final, "group", summarise, mean = mean(ratio),median   = median(ratio))

## get the cell summary for Figure 5

cellStats(g200_1deg_Raster, "sum") # 8987
cellStats(ghotspot_1deg_raster, "sum") # 11054
cellStats(GlastW_1deg_Raster, "sum")  # 6412
cellStats(sum_top17_lev_binary, "sum") # 3100

cellStats(GLWDPA_lev4_1deg_Raster, "sum") # 6945
cellStats(GLWDPA_lev6_1deg_Raster, "sum") # 8403

## top17
freq(overlay(g200_1deg_Raster, rank_3d_top17,fun=sum))  ## 1689
freq(overlay(ghotspot_1deg_raster, rank_3d_top17,fun=sum))  ## 1200
freq(overlay(GlastW_1deg_Raster, rank_3d_top17,fun=sum))  ## 253

freq(overlay(GLWDPA_lev4_1deg_Raster, rank_3d_top17,fun=sum))  ## 1012
freq(overlay(GLWDPA_lev4_1deg_Raster, g200_1deg_Raster,fun=sum))  ## 3436
freq(overlay(GLWDPA_lev4_1deg_Raster, ghotspot_1deg_raster,fun=sum))  ## 1765
freq(overlay(GLWDPA_lev4_1deg_Raster, GlastW_1deg_Raster,fun=sum))  ## 2162
freq(overlay(GLWDPA_lev4_1deg_Raster, g200_1deg_Raster,rank_3d_top17,fun=sum))  ## 897
freq(overlay(GLWDPA_lev4_1deg_Raster, ghotspot_1deg_raster,rank_3d_top17,fun=sum))  ## 687
freq(overlay(GLWDPA_lev4_1deg_Raster, GlastW_1deg_Raster,rank_3d_top17,fun=sum))  ## 144

## top 50

freq(overlay(g200_1deg_Raster, rank_3d_top50,fun=sum))  ## 4212
freq(overlay(ghotspot_1deg_raster, rank_3d_top50,fun=sum))  ## 2596
freq(overlay(GlastW_1deg_Raster, rank_3d_top50,fun=sum))  ## 875

freq(overlay(GLWDPA_lev4_1deg_Raster, rank_3d_top50,fun=sum))  ## 2635
freq(overlay(GLWDPA_lev4_1deg_Raster, g200_1deg_Raster,fun=sum))  ## 3436
freq(overlay(GLWDPA_lev4_1deg_Raster, ghotspot_1deg_raster,fun=sum))  ## 1765
freq(overlay(GLWDPA_lev4_1deg_Raster, GlastW_1deg_Raster,fun=sum))  ## 2162
# freq(overlay(g200_1deg_Raster, GlastW_1deg_Raster,fun=sum))

freq(overlay(GLWDPA_lev4_1deg_Raster, g200_1deg_Raster,rank_3d_top50,fun=sum))  ## 1926
freq(overlay(GLWDPA_lev4_1deg_Raster, ghotspot_1deg_raster,rank_3d_top50,fun=sum))  ## 1294
freq(overlay(GLWDPA_lev4_1deg_Raster, GlastW_1deg_Raster,rank_3d_top50,fun=sum))  ## 410

## transform them into percentage, and then plot them
df17 <- read.table(text = "
category	type	group	Percent
17prec	G200	outside	8.70
17prec	G200	WDPA	5.82
17prec	G200	washared	45.39
17prec	G200	template	40.08
17prec	hotspots	outside	22.82
17prec	hotspots	WDPA	16.45
17prec	hotspots	washared	34.77
17prec	hotspots	template	25.96
17prec	LastW	outside	43.27
17prec	LastW	WDPA	43.93
17prec	LastW	washared	7.29
17prec	LastW	template	5.52
 ", header = T)
 


p_df17 <- ggplot(df17, aes(x = type, y = Percent, fill = group)) + 
  geom_bar(position="stack", stat="identity") +
  geom_text(aes(label =Percent), color = "white", size = 5, position = position_stack(vjust = 0.5)) +
  scale_fill_viridis(discrete = T) +
  ggtitle("17%") +
  theme_article() +   xlab("")

df50 <- read.table(text = "category	type	group	Percent
50prec	G200	outside	16.58
50prec	G200	WDPA	12.02
50prec	G200	washared	32.65
50prec	G200	template	38.75
50prec	hotspots	outside	33.26
50prec	hotspots	WDPA	22.73
50prec	hotspots	washared	21.94
50prec	hotspots	template	22.07
50prec	LastW	outside	47.45
50prec	LastW	WDPA	37.72
50prec	LastW	washared	6.95
50prec	LastW	template	7.88


                                 ", header = T)
p_df50 <- ggplot(df50, aes(x = type, y = Percent, fill = group)) + 
  geom_bar(position="stack", stat="identity") +
  geom_text(aes(label =Percent), color = "white",size = 5, position = position_stack(vjust = 0.5)) +
  scale_fill_viridis(discrete = T) +
  ggtitle("50%") +
  theme_article() +
  xlab("")



pdf("Fig.5.pdf",  useDingbats=FALSE, width=8, height=8)

ggarrange(p_df17, p_df50, labels = c("A", "B"), ncol = 2, nrow = 1)

dev.off()

##### compare the two top17 or top50 layers

## top 17 comparison
## to obtain the inside range for the less 50% areas
plot(rank_3d_top17)
plot(sum_top17_lev_binary)
top17_inside <- mask(rank_3d_top17, sum_top17_lev_binary)  ### the shared part
plot(top17_inside)


## to obtian the outside range for the less 50% areas
top17_outside <- mask(rank_3d_top17, sum_top17_lev_binary, inverse = T)  ## only rank_3d_top17
plot(top17_outside)


top17_outside2 <- mask(sum_top17_lev_binary,rank_3d_top17,  inverse = T)  ## only sum_top17_lev_binary
plot(top17_outside2)

### top 50
plot(rank_3d_top50)
plot(sum_top50_lev_binary)
top50_inside <- mask(rank_3d_top50, sum_top50_lev_binary)  ### the shared part
plot(top50_inside)


## to obtian the outside range for the less 50% areas
top50_outside <- mask(rank_3d_top50, sum_top50_lev_binary, inverse = T)  ## only rank_3d_top17
plot(top50_outside)



top50_outside2 <- mask(sum_top50_lev_binary,rank_3d_top50,  inverse = T)  ## only sum_top17_lev_binary
plot(top50_outside2)

### I changed to the same percantage as the combined layers

### change the conservation maps into binary

rank_3d_top26 <- rank_3d_mask
rank_3d_top26[rank_3d_top26 < 0.737] <- NA
rank_3d_top26[rank_3d_top26 >= 0.737] <- 1
plot(rank_3d_top26)
g4 <- 1 * !is.na(rank_3d_top26)

plot(rank_3d_top26)
plot(sum_top17_lev_binary)
top26_inside <- mask(rank_3d_top26, sum_top17_lev_binary)  ### the shared part
plot(top26_inside)


## to obtian the outside range for the less 50% areas
top26_outside <- mask(rank_3d_top26, sum_top17_lev_binary, inverse = T)  ## only rank_3d_top17
plot(top26_outside)



top26_outside2 <- mask(sum_top17_lev_binary,rank_3d_top26,  inverse = T)  ## only sum_top17_lev_binary
plot(top26_outside2)
### top 61.5%
rank_3d_top61 <- rank_3d_mask
rank_3d_top61[rank_3d_top61 < 0.385] <- NA
rank_3d_top61[rank_3d_top61 >= 0.385] <- 1
plot(rank_3d_top61)
g4 <- 1 * !is.na(rank_3d_top61)

plot(rank_3d_top61)
plot(sum_top50_lev_binary)
top62_inside <- mask(rank_3d_top61, sum_top50_lev_binary)  ### the shared part
plot(top62_inside)


## to obtian the outside range for the less 50% areas
top62_outside <- mask(rank_3d_top61, sum_top50_lev_binary, inverse = T)  ## only rank_3d_top17
plot(top62_outside)



top62_outside2 <- mask(sum_top50_lev_binary,rank_3d_top61,  inverse = T)  ## only sum_top17_lev_binary
plot(top62_outside2)


## prepare and plot the proportition with the existing PAs and three frameworks

cellStats(rank_3d_top26, "sum") # 3076
cellStats(rank_3d_top61, "sum") # 7244

# writeRaster(sum_top50_lev, "sum_top50_lev.tiff", overwrite = TRUE)

cellStats(g200_1deg_Raster, "sum") # 8987
cellStats(ghotspot_1deg_raster, "sum") # 11054
cellStats(GlastW_1deg_Raster, "sum")  # 6412
cellStats(sum_top17_lev_binary, "sum") # 3100

cellStats(GLWDPA_lev4_1deg_Raster, "sum") # 6945

## top17
freq(overlay(g200_1deg_Raster, rank_3d_top26,fun=sum))  ## 2530
freq(overlay(ghotspot_1deg_raster, rank_3d_top26,fun=sum))  ## 1692
freq(overlay(GlastW_1deg_Raster, rank_3d_top26,fun=sum))  ## 413

freq(overlay(GLWDPA_lev4_1deg_Raster, rank_3d_top26,fun=sum))  ## 1486


freq(overlay(GLWDPA_lev4_1deg_Raster, g200_1deg_Raster,rank_3d_top26,fun=sum))  ## 1258
freq(overlay(GLWDPA_lev4_1deg_Raster, ghotspot_1deg_raster,rank_3d_top26,fun=sum))  ##905
freq(overlay(GLWDPA_lev4_1deg_Raster, GlastW_1deg_Raster,rank_3d_top26,fun=sum))  ##225


## top 50

freq(overlay(g200_1deg_Raster, rank_3d_top61,fun=sum))  ## 4726
freq(overlay(ghotspot_1deg_raster, rank_3d_top61,fun=sum))  ## 2780
freq(overlay(GlastW_1deg_Raster, rank_3d_top61,fun=sum))  ## 1071

freq(overlay(GLWDPA_lev4_1deg_Raster, rank_3d_top61,fun=sum))  # 3240


freq(overlay(GLWDPA_lev4_1deg_Raster, g200_1deg_Raster,rank_3d_top61,fun=sum))  ## 2155
freq(overlay(GLWDPA_lev4_1deg_Raster, ghotspot_1deg_raster,rank_3d_top61,fun=sum))  ## 1374
freq(overlay(GLWDPA_lev4_1deg_Raster, GlastW_1deg_Raster,rank_3d_top61,fun=sum))  ## 526

## transform them into percentage, and then plot them


df26 <- read.table(text = "
category	type	group	Percent
17prec	G200	outside	10.34
17prec	G200	WDPA	7.41
17prec	G200	washared	40.90
17prec	G200	template	41.35
17prec	hotspots	outside	26.11
17prec	hotspots	WDPA	18.89
17prec	hotspots	washared	29.42
17prec	hotspots	template	25.59
17prec	LastW	outside	45.58
17prec	LastW	WDPA	40.99
17prec	LastW	washared	7.31
17prec	LastW	template	6.11



                                 ", header = T)
library(reshape2)

library(ggplot2)
library(viridis)
library(ggpubr)
library(egg)

p_df26 <- ggplot(df26, aes(x = type, y = Percent, fill = group)) + 
  geom_bar(position="stack", stat="identity") +
  geom_text(aes(label =Percent), color = "white", size = 5, position = position_stack(vjust = 0.5)) +
  scale_fill_viridis(discrete = T) +
  ggtitle("17%") +
  theme_article() +  ## library(egg)
  xlab("")

df62 <- read.table(text = "category	type	group	Percent
50prec	G200	outside	19.78
50prec	G200	WDPA	14.98
50prec	G200	washared	29.75
50prec	G200	template	35.49
50prec	hotspots	outside	35.86
50prec	hotspots	WDPA	25.76
50prec	hotspots	washared	18.97
50prec	hotspots	template	19.41
50prec	LastW	outside	47.75
50prec	LastW	WDPA	37.47
50prec	LastW	washared	7.26
50prec	LastW	template	7.52



                                 ", header = T)
p_df62 <- ggplot(df62, aes(x = type, y = Percent, fill = group)) + 
  geom_bar(position="stack", stat="identity") +
  geom_text(aes(label =Percent), color = "white",size = 5, position = position_stack(vjust = 0.5)) +
  scale_fill_viridis(discrete = T) +
  ggtitle("50%") +
  theme_article() +
  xlab("")



pdf("Fig.S14.pdf",  useDingbats=FALSE, width=8, height=12)

ggarrange( p_df26, p_df62,  ncol = 2, nrow = 2)

dev.off()

