library(tidyverse)
library(Momocs)
library(vroom)
library(ape)
library(mgsub)

f <- list.files("Feathers",pattern=".txt|.csv",full.names = TRUE)

out.df <- vroom::vroom(f, id = "filename")
#make list
outs.l <- sapply(f,function(x) out.df %>% filter(filename==x) %>% select(X,Y) %>% as.matrix)
summary(outs.l)

birdlist <- unique(mgsub(basename(out.df$filename), c("XY_",".txt"), c("","")))
view(birdlist)

filenames <-tibble(xy.file=out.df$filename, species=mgsub(basename(out.df$filename), c("XY_",".txt"), c("","")))

view(birdlist)

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
  left_join(filenames)

#comparative analysis

head(outs.pca$x,1)


#plotting
bird.tree <- ape::read.tree("bird_tree.tre") %>% ladderize()
bird.tree <- keep.tip(bird.tree, birdlist)
bird.tree$tip.label <- gsub("_"," ",bird.tree$tip.label)

##bird.tree$tip.label <- gsub("_"," ",bird.tree$tip.label)
plot(main = "Birds Phylogenetics Tree", bird.tree,cex=0.1)


basename(names(outs))[1:5]
head(outs.pca$x,1)

files_species <- out.df$filename



