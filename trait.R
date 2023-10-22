######################Seagrass trait data####################
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(vegan)
library(gawdis)
library(reshape2)
library(inspectdf)
library(vegan)
library(ggforce)
library(ggrepel)
library(dendextend)
library(circlize)
library(ape)
library(mgcv)
library(ggord)

#read data
traits=read.table("./data/traits.txt", sep="\t", header=TRUE) %>% glimpse()
ordersp=read.table("./data/order.txt", sep="\t", header=TRUE) %>% glimpse()
distances=read.table("./data/distancias_final.txt", sep="\t", header=TRUE) %>% glimpse()

#visual inspection of multi-trait data
traits %>% select(-Species) %>% inspect_cat() %>% show_plot()
ggsave("Fig.S1.tiff")
inspect_cor(traits) %>% show_plot()
ggsave("Fig.S2.tiff")
inspect_num(traits) %>% show_plot()

traits %>% group_by(Family) %>% count()

#Create three types of multi-trait dissimilarities
trait1 = gowdis(traits) %>% round(3)
class(trait1) #let's coerce the dist object into a data.frame
df.trait1=melt(as.matrix(trait1), varnames = c("row", "col")) %>% rename (value1 = value)
write.csv2(df.trait1, "traitdis.csv")

trait2=gawdis(traits[, 2:11], w.type="user") %>% round(3)
df.trait2=melt(as.matrix(trait2), varnames = c("row", "col")) %>% rename (value2 = value)

trait3 = gawdis(traits[, 2:11], w.type="equal", groups =c(1,1, 2, 3,3,3, 4, 5,5,5))
df.trait3=melt(as.matrix(trait3), varnames = c("row", "col")) %>% round(3) %>% rename (value3 = value)

#combine the 3 data.frames into one and check correlation between the 3 multi-trait dissimilarities
dfcorr=cbind(df.trait1, df.trait2, df.trait3) %>% select(value1, value2, value3) %>%cor()

#PCoa and customization via ggplot
pcoa=cmdscale (trait1, eig = TRUE)
ordiplot(pcoa, display = 'sites', type = 'text')
str(pcoa)
pcoa$points #extract scores to plot
scores=as.data.frame(pcoa$points) 
scores2=cbind(scores, traits) %>%rename(Dim1 = V1, Dim2 = V2)%>%
  select(Species, Family, Dim1,Dim2)

p=ggplot(scores2, aes(x = Dim1, y = Dim2), size = 6) + geom_point(aes(colour = factor(Family), size = 4, shape = factor(Family))) 
cols=c("Cymodoceaceae" = "black", "Hydrocharitaceae" = "magenta", "Zosteraceae" = "green", "Posidoniaceae" = "blue")
p + theme_bw(base_size = 14) +
  theme(legend.title=element_blank()) + 
  theme(legend.position = 'bottom') +
  labs(x = "PCoA I (22.3%)", y = "PCoA II (9.6%)") +
  scale_shape_manual(values = c(0, 1, 2, 6)) +
  scale_color_manual(values = cols)
ggsave("PCOA.tiff")

pcoa$eig #extract eigenvalues to assess variance explained
eig = as.data.frame(pcoa$eig)
eig %>% sum (pcoa$eig)
eig %>% mutate (porc = pcoa$eig/12.41)

#Cluster of multi-trait dissimilarities
clus=hclust(trait1) #agglomerative clustering using complete linkage
plot(clus)
str(clus)
order=as.data.frame(clus$order)
write.csv(order, "order.csv")

#let's match names and customize the dendogram
dend=as.dendrogram(hclust(trait1))
labels(dend)=ordersp$Species2
labels(dend)
plot(dend)
labels_colors(dend)=ordersp$Familycode
labels_colors(dend)
dend=set(dend, "labels_cex", 1)
plot(dend, horiz = TRUE)
plot(as.phylo(dend), cex = 0.8, label.offset = 0.01, type = "fan", tip.color = ordersp$Familycode)
ggsave("cluster.tiff")

#Assesing and plotting phylogenetic signal (overall and for each family)
distances = distances %>% filter(distrais != 0, disfil != 0) 
distances %>% select(distrais, disfil) %>% cor(method = c("pearson"))

distances %>% filter(Family1 == "Cymodoceaceae", Family2 == "Cymodoceaceae") %>% 
select(distrais, disfil) %>% cor(method = c("pearson"))

distances %>% filter(Family1 == "Hydrocharitaceae", Family2 == "Hydrocharitaceae") %>% 
  select(distrais, disfil) %>% cor(method = c("pearson"))

distances %>% filter(Family1 == "Zosteraceae", Family2 == "Zosteraceae") %>% 
  select(distrais, disfil) %>% cor(method = c("pearson")) 

distances %>% filter(Family1 == "Posidoniaceae", Family2 == "Posidoniaceae") %>% 
select(distrais, disfil) %>% cor(method = c("pearson")) 

ggplot(distances, aes(x = disfil, y = distrais)) +
  geom_point(size = 3) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) + 
  theme_bw(base_size = 18) +
  xlab("Phylogenetic distances") + ylab("Multi-trait distances") +
  geom_smooth()
ggsave("relation.tiff")

distances %>% filter(Family1 == "Zosteraceae", Family2 == "Zosteraceae") %>%
  ggplot(aes(x = disfil, y = distrais)) +
  geom_point(size = 3) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) + 
  theme_bw(base_size = 18) +
  xlab("Phylogenetic distances") + ylab("Multi-trait distances") +
  geom_smooth()

distances %>% filter(Family1 == "Posidoniaceae", Family2 == "Posidoniaceae") %>%
  ggplot(aes(x = disfil, y = distrais)) +
  geom_point(size = 3) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) + 
  theme_bw(base_size = 18) +
  xlab("Phylogenetic distances") + ylab("Multi-trait distances") +
  geom_smooth()

distances %>% filter(Family1 == "Cymodoceaceae", Family2 == "Cymodoceaceae") %>%
  ggplot(aes(x = disfil, y = distrais)) +
  geom_point(size = 3) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) + 
  theme_bw(base_size = 18) +
  xlab("Phylogenetic distances") + ylab("Multi-trait distances") +
  geom_smooth()

distances %>% filter(Family1 == "Hydrocharitaceae", Family2 == "Hydrocharitaceae") %>%
  ggplot(aes(x = disfil, y = distrais)) +
  geom_point(size = 3) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) + 
  theme_bw(base_size = 18) +
  xlab("Phylogenetic distances") + ylab("Multi-trait distances") +
  geom_smooth()
ggsave("multi-traitdistancesfam.tiff")

#Assesing and ploting multi-trait distances among families
zosteraTD = distances %>% group_by(Family1, Family2) %>% 
  summarise(meanTD = mean(distrais), SETD = sd(distrais)) %>%
  filter(Family1 == "Zosteraceae", Family2 == "Zosteraceae")

HydroTD = distances %>% group_by(Family1, Family2) %>% 
  summarise(meanTD = mean(distrais), SETD = sd(distrais)) %>%
  filter(Family1 == "Hydrocharitaceae", Family2 == "Hydrocharitaceae")  

PosiTD = distances %>% group_by(Family1, Family2) %>% 
  summarise(meanTD = mean(distrais), SETD = sd(distrais)) %>%
  filter(Family1 == "Posidoniaceae", Family2 == "Posidoniaceae")   
  
CymoTD = distances %>% group_by(Family1, Family2) %>% 
  summarise(meanTD = mean(distrais), SETD = sd(distrais)) %>%
  filter(Family1 == "Cymodoceaceae", Family2 == "Cymodoceaceae")   

TD_fam = rbind(CymoTD, PosiTD, HydroTD, zosteraTD)

ggplot(TD_fam, aes(x=Family1, y=meanTD, fill=Family1)) + 
  geom_bar(stat="identity") + geom_errorbar(aes(ymin = meanTD - SETD, ymax = meanTD + SETD)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) + 
  xlab(" ") + ylab("Multi-trait distances") +
  theme_bw(base_size = 16) 

#Assesing and ploting phylogenetic distances among families
zosteraPD = distances %>% group_by(Family1, Family2) %>% 
  summarise(meanPD = mean(disfil), SEPD = sd(disfil)) %>%
  filter(Family1 == "Zosteraceae", Family2 == "Zosteraceae")

HydroPD = distances %>% group_by(Family1, Family2) %>% 
  summarise(meanPD = mean(disfil), SEPD = sd(disfil)) %>%
  filter(Family1 == "Hydrocharitaceae", Family2 == "Hydrocharitaceae")  

PosiPD = distances %>% group_by(Family1, Family2) %>% 
  summarise(meanPD = mean(disfil), SEPD = sd(disfil)) %>%
  filter(Family1 == "Posidoniaceae", Family2 == "Posidoniaceae")   

CymoPD = distances %>% group_by(Family1, Family2) %>% 
  summarise(meanPD = mean(disfil), SEPD = sd(disfil)) %>%
  filter(Family1 == "Cymodoceaceae", Family2 == "Cymodoceaceae")   

PD_fam = rbind(CymoPD, PosiPD, HydroPD, zosteraPD)

ggplot(PD_fam, aes(x=Family1, y=meanPD, fill=Family1)) + 
  geom_bar(stat="identity") + geom_errorbar(aes(ymin = meanPD - SEPD, ymax = meanPD + SEPD)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) + 
  xlab(" ") + ylab("Phylogenetic distances") +
  theme_bw(base_size = 16) 

ggsave("phylodistancesfam.tiff")

#GAM assessing conection between multi-trait and phylogenetic distances
model_gam=gam(distrais ~ s(disfil, k = 3, bs = "cr"),
                    family=poisson(link = log), data=distances)
summary(model_gam)
plot(model_gam)
gam.check(model_gam)

#get correlations between PCoA scores and traits
correlation=read.table("./data/correlations_PCOA.txt", sep="\t", header=TRUE) %>% glimpse()
library(corrplot)

correlation1=correlation %>% select(-Dim2, -Species)
correlation2=correlation %>% select(-Dim1, -Species) 
cor1=cor(correlation1) %>% round(2) %>% view()
cor2=cor(correlation2) %>% round(2) %>% view()
write.csv(cor1, "cor1.csv")
write.csv(cor2, "cor2.csv")