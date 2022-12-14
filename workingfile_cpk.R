library(dplyr)
library(tidyverse)
library(ape)
library(Momocs)
library(phytools)
library(RRphylo)
library(phytools)

setwd("~/Downloads")

##making the dataset
birdlist <- read.table(file = "datalist.txt", header = TRUE, sep = ",") %>% 
  as_tibble()
birdlist <- birdlist[,-5]

#note you were trying to compute mean of the common name. removed.
meanbirdlist <- birdlist %>% 
  group_by(scientific_name) %>% 
  summarise_at(vars(barb_distance, barbule_distance),mean) %>% 
  mutate(ratio=barbule_distance/barb_distance)

#no need for this?
#meanbirdlist$common_name <- birdlist$common_name

# commonnamelist <- data.frame(birdlist$scientific_name, birdlist$common_name)
# commonnamelist <- commonnamelist[!duplicated(commonnamelist),]
# commonnamelist <- rename(commonnamelist, scientific_name = birdlist.scientific_name, common_name = birdlist.common_name)
# 
# finalbirdlist <- merge (meanbirdlist, commonnamelist, by ="scientific_name" )
# finalbirdlist <- finalbirdlist[,-2]
# finalbirdlist <-relocate(finalbirdlist, common_name, .after = scientific_name)
# finalbirdlist$ratio <-(finalbirdlist$barbule_distance/finalbirdlist$barb_distance)
# finalbirdlist <- rename(finalbirdlist, species = scientific_name)
# 
# rownames(finalbirdlist) <- finalbirdlist$species

##tree

#load the many trees from the MCMC analysis
bird.trees <- read.tree("AllBirdsEricson1.tre")

#trim all the trees to the taxa included
trimmed.trees <- lapply(bird.trees,function(x) keep.tip(x, meanbirdlist$scientific_name))

#make a consensus tree
con.tree <- consensus.edges(trimmed.trees)
#must be rooted to work in RRphylo
con.tree <- root(con.tree,"Sayornis_phoebe")

plot(con.tree,cex=0.5)

##evolutionary rate
ratio <- meanbirdlist$ratio
names(ratio) <-meanbirdlist$scientific_name

birdRR <- RRphylo(tree=con.tree, ratio)

