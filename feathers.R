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
bird.trees <- read.tree("AllBirdsEricson1.tre")
trimmed.trees <- lapply(bird.trees,function(x) keep.tip(x, birdlist$species))
con.tree <- consensus.edges(trimmed.trees)
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

outs.PC1.SS<- search.shift(RR=outs.PC1.RR,status.type="clade")
outs.PC1.SS$single.clades

plot(con.tree, cex=0.2)
nodelabels(node = as.numeric(rownames(outs.PC1.SS$single.clades)),text = rownames(outs.PC1.SS$single.clades))

outs.PC1.plot <- plotShift(RR=outs.PC1.RR,SS=outs.PC1.SS)

outs.PC1.plot$plotClades()

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")



##custom function
plot_SS <- function(tre=NULL,SS=NULL,tax=NULL){
  
  
  nodes <- as.numeric(rownames(SS$single.clades))
  
  pal <- wes_palette("Zissou1",n=length(nodes))
  sp <- list()
  for(i in nodes){
    sp.i <- extract.clade(tre,i)$tip.label
    
    #print(head(tax))
    sub.names <- lapply(tax,function(x) x[x%in%sp.i]) 
    
    in.clades <- lapply(sub.names,function(x) length(x)>0) 
    all.of.clade <- lapply(sub.names,function(x) all(sapply(sp.i,function(z) z%in%x))) 
    
    high.clade <- names(sub.names)[last(which(all.of.clade==T))]
    all.clades <- names(sub.names)[which(in.clades==T)]
    crown <- ""
    if(high.clade!=last(names(sub.names))) crown <- "crown-"
    
    sub.clades <- NULL
    if(length(grepl("oidea",all.clades))>0) sub.clades <- all.clades[grepl("oidea",all.clades)]
    
    high.clade2 <- paste0(crown,high.clade,": ",paste0(sub.clades,collapse = "+"))
    sp[[paste0(i)]] <- tibble(n=i,species=sp.i,clade=high.clade2)
    
  }
  
  
  d <- do.call(rbind,sp)%>% 
    rename(label=species) 
  
  d2<- d %>% rename(clade_name=clade) 
  
  p <- ggtree(tre)+ scale_y_reverse()
  
  p$data <- p$data %>% left_join(d) %>% left_join(tibble(node=nodes,SS$single.clades) %>% mutate(shift=ifelse(rate.difference>0,"+","-")))
  
  p <-  p+geom_tiplab(aes(col=clade),geom="text",size=1.2)+
    geom_cladelab(data=d2,mapping=aes(node=n,col=clade_name,label=clade_name),offset=1,size=1.5)+
    geom_hilight(data=d2,mapping = aes(node = n,fill=clade_name),alpha = 0.01)+
    scale_fill_manual(values = pal)+
    scale_color_manual(values = pal)+
    theme(legend.position = "none")+geom_nodepoint(mapping=aes(subset = shift =="-"), size=5, shape=25,fill='blue',color='blue',alpha=0.7)+
    geom_nodepoint(mapping=aes(subset = shift =="+"), size=5, shape=24, fill='red',color='red',alpha=0.7)
  p <- p+xlim(NA,6)
  res <- tibble(n=nodes,SS$single.clades) %>% left_join(d %>% select(n,clade) %>% unique)
  
  return(list(plot=p,res=res))
  
}

outs.PC1.res <- plot_SS(con.tree,outs.PC1.SS,tax = birdlist$species)
