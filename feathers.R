library(tidyverse)
library(Momocs)
library(vroom)
library(ape)
library(mgsub)
library(phytools)
library(RRphylo)
library(ggtree)
library(wesanderson)

f <- list.files("Feathers",pattern=".txt|.csv",full.names = TRUE)

out.df <- vroom::vroom(f, id = "filename")
#make list
outs.l <- sapply(f,function(x) out.df %>% filter(filename==x) %>% select(X,Y) %>% as.matrix)
summary(outs.l)

birdlist <- tibble(xy.file=unique(basename(out.df$filename)), species=unique(mgsub(basename(out.df$filename), c("XY_",".txt"), c("",""))))

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
con.tree <- read.tree("contree.tre")
con.tree <- root(con.tree,"Sayornis_phoebe")
plot(main = "Birds Phylogenetics Tree", con.tree,cex=0.2)

view(con.tree$tip.label)
#PC1s

outs.pc1 <- outs.pca %>% 
  filter(species%in% con.tree$tip.label) %>% 
  group_by(species) %>% 
  summarize(PC1=mean(PC1)) %>% 
  pull

names(outs.pc1) <-  outs.pca%>% 
  filter(species%in% con.tree$tip.label) %>% 
  group_by(species) %>% 
  summarize(PC1=mean(PC1)) %>% 
  pull(species)


#PC2s
outs.pc2 <- outs.pca %>% 
  filter(species%in% con.tree$tip.label) %>% 
  group_by(species) %>% 
  summarize(PC2=mean(PC2)) %>% 
  pull(PC2)

names(outs.pc2) <-  outs.pca%>% 
  filter(species%in% con.tree$tip.label) %>% 
  group_by(species) %>%
  summarize(PC2=mean(PC2)) %>% 
  pull(species)

outs.PC1.BM<-brownie.lite(con.tree,outs.pc1*10)

outs.PC2.BM<-brownie.lite(con.tree,outs.pc2*10)

#shifts in evolutionary rate
outs.PC1.RR <- RRphylo(tree=con.tree,y=outs.pc1)
outs.PC2.RR <- RRphylo(tree=con.tree,y=outs.pc2)

outs.PC1.SS<- search.shift(RR=outs.PC1.RR,status.type="clade")
outs.PC1.SS$single.clades

outs.PC2.SS<- search.shift(RR=outs.PC2.RR,status.type="clade")
outs.PC2.SS$single.clades

plot(con.tree, cex=0.2)
nodelabels(node = as.numeric(rownames(outs.PC1.SS$single.clades)),text = rownames(outs.PC1.SS$single.clades))

plot(con.tree, cex=0.2)
nodelabels(node = as.numeric(rownames(outs.PC2.SS$single.clades)),text = rownames(outs.PC2.SS$single.clades))

outs.PC1.plot <- plotShift(RR=outs.PC1.RR,SS=outs.PC1.SS)
outs.PC1.plot$plotClades()

outs.PC2.plot <- plotShift(RR=outs.PC2.RR,SS=outs.PC2.SS)
outs.PC2.plot$plotClades()

