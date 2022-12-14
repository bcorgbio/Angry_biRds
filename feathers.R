library(tidyverse)
library(Momocs)
library(vroom)
library(ape)
library(mgsub)
library(phytools)
library(RRphylo)

f <- list.files("Feathers",pattern=".txt|.csv",full.names = TRUE)

out.df <- vroom::vroom(f, id = "filename")
#make list
outs.l <- sapply(f,function(x) out.df %>% filter(filename==x) %>% select(X,Y) %>% as.matrix)
summary(outs.l)

birdlist <- tibble(xy.file=unique(basename(out.df$filename)), species=unique(mgsub(basename(out.df$filename), c("XY_",".txt","_"), c("",""," "))))

outs.l %>% 
  Out() %>% 
  coo_flipx()

##

#make a large df with vroom
out.df <- vroom::vroom(f, id = "filename") %>% 
  mutate(species=gsub("XY_.+_\\..+","\\1",basename(filename))) %>% 
  na.omit()

#make list
outs.l <- sapply(f,function(x) out.df %>% filter(filename==x) %>% select(X,Y) %>% as.matrix)

#extract wing info
species <- gsub("XY_.+_\\..+","\\1",basename(names(outs.l)))


outs <-  outs.l %>% 
  Out(fac=list(species=species)) %>% 
  coo_flipx()
#pipes the list of matrices into Our() function, flip, and visualize with stack

outs.min <- outs %>% 
  coo_nb() %>% 
  min()

outs %>% 
  coo_interpolate(outs.min) %>% 
  coo_slide(id=1) %>% 
  coo_align()  %>%
  fgProcrustes() %>% 
  stack()


outs.pca <- outs %>%
  coo_interpolate(outs.min) %>%
  coo_align()  %>%
  coo_slide(id=1) %>% 
  fgProcrustes() %>% 
  efourier(norm=FALSE) %>% 
  PCA()

outs.pca %>% 
  plot_PCA(title = "Primary Feathers")

outs.pca <-  tibble(xy.file=basename(rownames(outs.pca$x)),PC1=outs.pca$x[,1],PC2=outs.pca$x[,2])%>% 
  left_join(birdlist)


#comparative analysis

head(outs.pca$x,1)


#plotting
bird.tree <- ape::read.tree("bird_tree.tre") %>% ladderize()
bird.tree$tip.label <- gsub("_"," ",bird.tree$tip.label)
bird.tree <- keep.tip(bird.tree, birdlist$species)


plot(main = "Birds Phylogenetics Tree", bird.tree,cex=0.1)

#PC1s

outs.pc1 <- outs.pca %>% 
  filter(species%in% bird.tree$tip.label) %>% 
  group_by(species) %>% 
  summarize(PC1=mean(PC1)) %>% 
  pull

names(outs.pc1) <-  outs.pca%>% 
  filter(species%in% bird.tree$tip.label) %>% 
  group_by(species) %>% 
  summarize(PC1=mean(PC1)) %>% 
  pull(species)


#PC2s
outs.pc2 <- outs.pca %>% 
  filter(species%in% bird.tree$tip.label) %>% 
  group_by(species) %>% 
  summarize(PC2=mean(PC2)) %>% 
  pull(PC2)

names(outs.pc2) <-  outs.pca%>% 
  filter(species%in% bird.tree$tip.label) %>% 
  group_by(species) %>%
  summarize(PC2=mean(PC2)) %>% 
  pull(species)

outs.PC1.BM<-brownie.lite(bird.tree,outs.pc1*10)

outs.PC2.BM<-brownie.lite(bird.tree,outs.pc2*10)

#shifts in evolutionary rate
outs.PC1.RR <- RRphylo(tree=bird.tree,y=outs.pc1)


