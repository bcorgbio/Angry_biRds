library(ape)
library(tidyverse)

bird.tree <- ape::read.tree("bird_tree.tre") %>% ladderize()
bird.tree$tip.label <- gsub("_"," ",bird.tree$tip.label)
plot(main = "Birds Phylogenetics Tree", bird.tree,cex=0.1)

