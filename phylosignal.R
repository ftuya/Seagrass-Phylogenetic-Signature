###############Phylogenetic Signal (PS) of seagrasses. Code by F. Tuya (October 2022)#############
library(tidyverse)
library(ape)
library(adephylo)
library(phylobase)
library(phylosignal)
library(vegan)
library(ggord)
library(ggvegan)

################Initial PCA on reproduction traits#################
reproPCA=read.table("./data/PCA_repro.txt", sep="\t", header=TRUE) %>% glimpse()
pca=prcomp(reproPCA[, 2:5], scale=TRUE)
summary(pca)
biplot(pca, scale=0) 
str(pca)
pca$x
pca2=cbind(reproPCA, pca$x[, 1:2])
head(pca2)
write.csv2(pca2, "pca.csv")

plotpca=ggord(pca)
plotpca + xlim(-2,2.8) + ylim(-1.8,1.6) +
  theme(text = element_text(size = 20)) +
  xlab("PC1") + ylab("PC2")
ggsave("PCA.tiff")

################Analysis of PS####################################
tree=read.nexus("./data/tree_seagrass.nex") #read tree in nexus format
write.nexus(tree) # check names of phylogeny; name in the phylogeny should be equal to those of the trait dataframe
write.csv2(tree$tip.label, "tips.csv")
traits=read.table("./data/phylo_traits.txt", sep="\t", header=TRUE) %>% glimpse() #read traits dataframe
rownames(traits)=traits[, 1]
traits[, 1]=NULL
p4d=phylo4d(tree, traits) #pool both data (tree and traits)

#Locating and visualizing the PS with LIPA
seagrass.lipa=lipaMoran(p4d)
seagrass.lipa.p4d=lipaMoran(p4d, as.p4d = TRUE)
barplot.phylo4d(p4d, bar.col=(seagrass.lipa$p.value < 0.05) + 1, center = FALSE , scale = FALSE)
barplot.phylo4d(seagrass.lipa.p4d, bar.col = (carni.lipa$p.value < 0.05) + 1, center = FALSE, scale = FALSE)

#other visualizations of PS
barplot.phylo4d(p4d, center = FALSE , scale = FALSE) #Fig. S5
barplot.phylo4d(p4d, tree.type = "phylo", tree.ladderize = TRUE) #Fig. S5
gridplot(p4d)
gridplot(p4d, tree.type = "fan", tip.cex = 0.6, show.trait = FALSE)

#Getting and exporting the metrics of PS
phylosignal=phyloSignal(p4d) 
write.csv2(phylosignal$stat, "stat.csv")
write.csv2(phylosignal$pvalue, "p.csv")

#Assessing the behavior of these methods with this phylogeny along a Brownian-Motion influence gradient
phylosim=phyloSim(tree = tree, method = "all", nsim = 100, reps = 99) #Fig. S3
plot(phylosim, stacked.methods = FALSE, quantiles = c(0.05, 0.95))

#Phylogenetic correlogram for particular traits 
Mindepth.cg=phyloCorrelogram(p4d, trait = "Min.depth")
plot(Mindepth.cg)

Maxdepth.cg=phyloCorrelogram(p4d, trait = "Max.depth")
plot(Maxdepth.cg)

Generation.cg=phyloCorrelogram(p4d, trait = "Generation")
plot(Generation.cg)

leafmaxlenght.cg=phyloCorrelogram(p4d, trait = "Leaf.max.lenght")
plot(leafmaxlenght.cg)

leafmaxwidth.cg=phyloCorrelogram(p4d, trait = "Leaf.max.width")
plot(leafmaxwidth.cg)

repropc1.cg=phyloCorrelogram(p4d, trait = "repro.pc1")
plot(repropc1.cg)

repropc2.cg=phyloCorrelogram(p4d, trait = "repro.pc2")
plot(repropc2.cg)

par(mfrow=c(2,4))
plot(Mindepth.cg)
plot(Maxdepth.cg)
plot(Generation.cg)
plot(leafmaxlenght.cg)
plot(leafmaxwidth.cg)
plot(repropc1.cg)
plot(repropc2.cg)

################ Figures Appendix ####################################
traits_leaf=read.table("./data/leaf_traits.txt", sep="\t", header=TRUE) %>% glimpse() #read leaf traits dataframe
library(ggstatsplot)
library(inspectdf)
traits_leaf$Leaf=as.factor(traits_leaf$Leaf)

traits_leaf %>% # Fig. S1 leaf traits
  mutate(Leaf = fct_reorder(Leaf, desc(Leaf.max.lenght))) %>%
  ggbetweenstats(x = Leaf,y = Leaf.max.lenght) +
  labs(y="Leaf maximum lenght", x=" ")
ggsave("Fig.S1.tiff")

traits_leaf %>%
  mutate(Leaf = fct_reorder(Leaf, desc(Leaf.max.width))) %>%
  ggbetweenstats(x = Leaf,y = Leaf.max.lenght) +
  labs(y="Leaf maximum widht", x=" ") 

traits %>% # Fig. S2
  inspect_cor() %>% show_plot() 
ggsave("Fig.S2.tiff")
