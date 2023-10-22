################################################################################
######### Tuya et al. Phylogenetic signal - Evolutionary models Seagrass #######

## Clean working directory

rm(list=ls())

## Load libraries

library(dplyr)
library(tidyr)
library(ggplot2)
library(picante)
library(ape)

### Set plotting defaults ----

Theme1 <-
  theme( # use theme_get() to see available options
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill="white"),
    legend.key = element_blank(), # switch off the rectangle around symbols in the legend
    legend.text = element_text(size=10),
    legend.title = element_text(size=14, face="bold"),
    legend.position = "right",
    legend.direction="vertical",
    text=element_text(size=14),
    strip.text.y = element_text(size = 14,angle = 0),
    axis.title.x=element_text(vjust=0.3, size=12),
    axis.title.y=element_text(vjust=0.6, angle=90, size=12),
    axis.text.x=element_text(size=10,angle = 0),
    axis.text.y=element_text(size=10),
    axis.line.x=element_line(colour="black", size=0.5,linetype='solid'),
    axis.line.y=element_line(colour="black", size=0.5,linetype='solid'),
    plot.title = element_text(size=14, face="bold"),
    strip.background = element_blank())

### Function to calculate standard errors ----

se <- function(x) sd(x) / sqrt(length(x))
se.min <- function(x) (mean(x)) - se(x)
se.max <- function(x) (mean(x)) + se(x)

# Set work directory----

work.dir=("C:/Users/22373243/Dropbox/Projects/Analysis/Analysis_Seagrass_Phylosignal") ### Desktop

# Set sub directories----

data.dir=paste(work.dir,"Data",sep="/")
plots.dir=paste(work.dir,"Plots",sep="/")
model.out=paste(work.dir,"ModelOut",sep="/")

## Bring in trait data

setwd(data.dir)
dir()

traits<-read.delim("phylo_traits.txt")%>%
  glimpse()

## Explore NAs in the database

sum(is.na(traits))/prod(dim(traits))*100
apply(traits,2,function(col)sum(is.na(col))/length(col))*100 ## 1.53% of biomass values missing

## Code traits

glimpse(traits)

traits$Min.depth<-as.numeric(traits$Min.depth)
traits$Max.depth<-as.numeric(traits$Max.depth)
traits$Generation<-as.numeric(traits$Generation)
traits$Leaf.max.lenght<-as.numeric(traits$Leaf.max.lenght)
traits$Leaf.max.width<-as.numeric(traits$Leaf.max.width)
traits$repro.pc1<-as.numeric(traits$repro.pc1)
traits$repro.pc2<-as.numeric(traits$repro.pc2)

## Add species names as rownames

rownames(traits) <- traits[,1] #Assigning row names from 1st column 
traits[,1] <- NULL #Removing the first column

## Bring in Phylogenetic tree

setwd(data.dir)
dir()

phy=read.nexus("tree.nex")

## Check species in trait and phylo match

x<-as.vector(rownames(traits))
y<-as.vector(phy$tip.label)

setdiff(x,y)
setdiff(y,x)

## Plot tree for exploration

plot(phy)
axisPhylo()
nodelabels()


################################################################################
##### What is the most appropiate evoluionary model for Seagrass Traits #####33

library(phytools)
library(geiger)
library(pmc)
library(reshape)

## Remove polytomies in the tree

phy_mod<-multi2di(phy,0.03) #remove polytomies
is.binary(phy_mod)
plot(phy_mod)

## Convert to an ultrametric tree

phy_mod<-force.ultrametric(phy_mod)
is.ultrametric(phy_mod)
plot(phy_mod)
is.binary(phy_mod)

#### (A) Univariate method - no trait covariance accounted for ----

## Store individual traits as vectors

names(traits)

MaxLength=traits[,"Leaf.max.lenght"]
names(MaxLength)=row.names(traits)

Mindepth=traits[,"Min.depth"]
names(Mindepth)=row.names(traits)

Maxdepth=traits[,"Max.depth"]
names(Maxdepth)=row.names(traits)

Generation=traits[,"Generation"]
names(Generation)=row.names(traits)

MaxWidth=traits[,"Leaf.max.width"]
names(MaxWidth)=row.names(traits)

repro1=traits[,"repro.pc1"]
names(repro1)=row.names(traits)

repro2=traits[, "repro.pc2"]
names(repro2)=row.names(traits)

## !Important! - trait data might need to be transformed to logarithmic scale to reduced variance between species

## (A) Max Length ----

brown_MaxLength=fitContinuous(phy_mod,MaxLength)

eb_MaxLength=fitContinuous(phy_mod,MaxLength,model = "EB")

ou_MaxLength=fitContinuous(phy_mod,MaxLength,model = "OU")

modelMaxLength=matrix(4,3,2,dimnames = list(c("Brownian Motion","Early Burst","Ornstein-Uhlenbeck"),c("log likelihood","AICc")))
modelMaxLength[,1]=c(brown_MaxLength$opt$lnL,eb_MaxLength$opt$lnL,ou_MaxLength$opt$lnL)
modelMaxLength[,2]=c(brown_MaxLength$opt$aicc,eb_MaxLength$opt$aicc,ou_MaxLength$opt$aicc)

## Plotting modes of trait evolution along the phylogeny

phenogram(phy_mod,MaxLength,fsize = 0.6,spread.costs=c(1,0),main = "Leaf Maximum Length")

obj1=contMap(phy_mod,MaxLength,plot=FALSE)
plot(obj1,legend=0.7*max(nodeHeights(phy_mod)),
     fsize=c(0.7,0.9))

## (B) Min Depth ----

brown_Mindepth=fitContinuous(phy_mod,Mindepth)

eb_Mindepth=fitContinuous(phy_mod,Mindepth,model = "EB")

ou_Mindepth=fitContinuous(phy_mod,Mindepth,model = "OU")

modelMindepth=matrix(4,3,2,dimnames = list(c("Brownian Motion","Early Burst","Ornstein-Uhlenbeck"),c("log likelihood","AICc")))
modelMindepth[,1]=c(brown_Mindepth$opt$lnL,eb_Mindepth$opt$lnL,ou_Mindepth$opt$lnL)
modelMindepth[,2]=c(brown_Mindepth$opt$aicc,eb_Mindepth$opt$aicc,ou_Mindepth$opt$aicc)

## Plotting modes of trait evolution along the phylogeny

phenogram(phy_mod,Mindepth,fsize = 0.6,spread.costs=c(1,0),main = "Minimum Depth")

obj1=contMap(phy_mod,Mindepth,plot=FALSE)
plot(obj1,legend=0.7*max(nodeHeights(phy_mod)),
     fsize=c(0.7,0.9))

## (C) Max Depth ----

brown_Maxdepth=fitContinuous(phy_mod,Maxdepth)

eb_Maxdepth=fitContinuous(phy_mod,Maxdepth,model = "EB")

ou_Maxdepth=fitContinuous(phy_mod,Maxdepth,model = "OU")

modelMaxdepth=matrix(4,3,2,dimnames = list(c("Brownian Motion","Early Burst","Ornstein-Uhlenbeck"),c("log likelihood","AICc")))
modelMaxdepth[,1]=c(brown_Maxdepth$opt$lnL,eb_Maxdepth$opt$lnL,ou_Maxdepth$opt$lnL)
modelMaxdepth[,2]=c(brown_Maxdepth$opt$aicc,eb_Maxdepth$opt$aicc,ou_Maxdepth$opt$aicc)

## Plotting modes of trait evolution along the phylogeny

phenogram(phy_mod,Maxdepth,fsize = 0.6,spread.costs=c(1,0),main = "Minimum Depth")

obj1=contMap(phy_mod,Maxdepth,plot=FALSE)
plot(obj1,legend=0.7*max(nodeHeights(phy_mod)),
     fsize=c(0.7,0.9))

## (D) Generation ----

brown_Generation=fitContinuous(phy_mod,Generation)

eb_Generation=fitContinuous(phy_mod,Generation,model = "EB")

ou_Generation=fitContinuous(phy_mod,Generation,model = "OU")

modelGeneration=matrix(4,3,2,dimnames = list(c("Brownian Motion","Early Burst","Ornstein-Uhlenbeck"),c("log likelihood","AICc")))
modelGeneration[,1]=c(brown_Generation$opt$lnL,eb_Generation$opt$lnL,ou_Generation$opt$lnL)
modelGeneration[,2]=c(brown_Generation$opt$aicc,eb_Generation$opt$aicc,ou_Generation$opt$aicc)

## Plotting modes of trait evolution along the phylogeny

phenogram(phy_mod,Generation,fsize = 0.6,spread.costs=c(1,0),main = "Minimum Depth")

obj1=contMap(phy_mod,Generation,plot=FALSE)
plot(obj1,legend=0.7*max(nodeHeights(phy_mod)),
     fsize=c(0.7,0.9))

## (E) Leaf Maximum Width ----

brown_MaxWidth=fitContinuous(phy_mod,MaxWidth)

eb_MaxWidth=fitContinuous(phy_mod,MaxWidth,model = "EB")

ou_MaxWidth=fitContinuous(phy_mod,MaxWidth,model = "OU")

modelMaxWidth=matrix(4,3,2,dimnames = list(c("Brownian Motion","Early Burst","Ornstein-Uhlenbeck"),c("log likelihood","AICc")))
modelMaxWidth[,1]=c(brown_MaxWidth$opt$lnL,eb_MaxWidth$opt$lnL,ou_MaxWidth$opt$lnL)
modelMaxWidth[,2]=c(brown_MaxWidth$opt$aicc,eb_MaxWidth$opt$aicc,ou_MaxWidth$opt$aicc)

## Plotting modes of trait evolution along the phylogeny

phenogram(phy_mod,MaxWidth,fsize = 0.6,spread.costs=c(1,0),main = "Minimum Depth")

obj1=contMap(phy_mod,MaxWidth,plot=FALSE)
plot(obj1,legend=0.7*max(nodeHeights(phy_mod)),
     fsize=c(0.7,0.9))

## (F) Repro PC1 ----

brown_repro1=fitContinuous(phy_mod,repro1)

eb_repro1=fitContinuous(phy_mod,repro1,model = "EB")

ou_repro1=fitContinuous(phy_mod,repro1,model = "OU")

modelrepro1=matrix(4,3,2,dimnames = list(c("Brownian Motion","Early Burst","Ornstein-Uhlenbeck"),c("log likelihood","AICc")))
modelrepro1[,1]=c(brown_repro1$opt$lnL,eb_repro1$opt$lnL,ou_repro1$opt$lnL)
modelrepro1[,2]=c(brown_repro1$opt$aicc,eb_repro1$opt$aicc,ou_repro1$opt$aicc)

## Plotting modes of trait evolution along the phylogeny

phenogram(phy_mod,repro1,fsize = 0.6,spread.costs=c(1,0),main = "Minimum Depth")

obj1=contMap(phy_mod,repro1,plot=FALSE)
plot(obj1,legend=0.7*max(nodeHeights(phy_mod)),
     fsize=c(0.7,0.9))

## (G) Repro PC2 ----

brown_repro2=fitContinuous(phy_mod,repro2)

eb_repro2=fitContinuous(phy_mod,repro2,model = "EB")

ou_repro2=fitContinuous(phy_mod,repro2,model = "OU")

modelrepro2=matrix(4,3,2,dimnames = list(c("Brownian Motion","Early Burst","Ornstein-Uhlenbeck"),c("log likelihood","AICc")))
modelrepro2[,1]=c(brown_repro2$opt$lnL,eb_repro2$opt$lnL,ou_repro2$opt$lnL)
modelrepro2[,2]=c(brown_repro2$opt$aicc,eb_repro2$opt$aicc,ou_repro2$opt$aicc)

## Plotting modes of trait evolution along the phylogeny

phenogram(phy_mod,repro2,fsize = 0.6,spread.costs=c(1,0),main = "Minimum Depth")

obj1=contMap(phy_mod,repro2,plot=FALSE)
plot(obj1,legend=0.7*max(nodeHeights(phy_mod)),
     fsize=c(0.7,0.9))

#### (B) Multivariate method - trait covariance accounted for ----

library(mvMORPH)
browseVignettes(package = "mvMORPH")

mvBM(phy_mod,traits)
mvOU(phy_mod,traits)
mvEB(phy,traits)

