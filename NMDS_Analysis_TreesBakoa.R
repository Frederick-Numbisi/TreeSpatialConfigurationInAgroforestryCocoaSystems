
library(vegan)
library(ggplot2)
library(grid)



#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
# NDMS for tree density with sampled CAFS

BaTDen = read.csv(file = "D:/PhD Manuscript/PhD Article Manuscripts/ArticleShapefiles/BakoaShapefiles/BakoaHyperframes/BTDens_nmds.csv", 
                  header = TRUE, check.names = FALSE) #, row.names = 1
names(BaTDen)
#View(BaTDen)




BaTDen$Plots = as.factor(BaTDen$Plots)
BaTDen$Group = as.factor(BaTDen$Group)
BaTDen$Origin = as.factor(BaTDen$Origin)
BaTDen$AgeOri = as.factor(BaTDen$AgeOri)

# Reorder the levels of the factor variable
BaTDen$Group = ordered(BaTDen$Group, 
                       levels = c("<= 10yrs","11 - 20yrs","21 - 40yrs","41 - 60yrs","  > 60yrs  "))

#BaTDen = data.frame(BaTDen)

# B) NMDS of Composition of tree structure

BaTDen_c = as.data.frame(BaTDen[6:11])

# Computing the z scores for the dataframe (remove outliers with Z-scores >= 3)
library(matrixStats)
BaTDen_c_zscores = (BaTDen_c  - rowMeans(BaTDen_c ))/(rowSds(as.matrix(BaTDen_c )))[row(BaTDen_c )] 
View(BaTDen_c_zscores)







# Rarify data down to the lowest tree density ("lowest-common-denomenator") so that we can compare evenly
# between samples regardles of sampling depth artefacts.

min_depth_D = min(colSums(BaTDen_c)) # Get the density for the lowest abundance sample 
min_depth_D
BaTDen_c_rarefied = as.data.frame(round(BaTDen_c, min_depth_D))# Rarefy the data to the lowest abundance sample


# Determining the best method for calculating the distance matrix from our data

#library(vegan)

sqrt_BaTDen_c_rarefied = sqrt(BaTDen_c_rarefied)
rank.BaTDen = rankindex(as.matrix(sqrt_BaTDen_c_rarefied), BaTDen_c_rarefied, indices = c("bray", "euclid", "manhattan", "horn"), method = "spearman")
print(paste("The highest rand was give by the", names(sort(rank.BaTDen, decreasing = TRUE)[1]), "method."))

# Compute the community distance matrix using the "Bray-Curtis" method

BaTDen_dis = as.matrix((vegdist(BaTDen_c_rarefied, "bray")))

BaTDen_dis.mds = metaMDS(BaTDen_dis, k=2) # Perform NMDS
BaTDen_dis.mds$stress




# # Etract scores
# sites = scores(BaTDen_dis.mds, display = "sites")
# species = scores(BaTDen_dis.mds, display = "species")
# 
# # Compute the axis ranges
# 
# xlim = range(sites[,1], species[,1])
# ylim = range(sites[,2], species[,2])
# 
# # set up a blank plotting window
# plot(BaTDen_dis.mds, type = "n", xlim = xlim, ylim = ylim)
# 
# # now adding the points (plots) and species (variables) on the blank plot
# 
# points(sites, col="blue", pch=16, cex=1.2)
# text(species, labels = rownames(species), col = 'red', cex = 0.6)

#------------------------------------------------------
# Build a dataframe with NMDS coordinates and metadata

MDS1 = BaTDen_dis.mds$points[,1]
MDS2 = BaTDen_dis.mds$points[,2]

NMDS_Den = data.frame(MDS1 = MDS1, MDS2 = MDS2, Group = BaTDen$Group, Origin = BaTDen$Origin, AgeOri = BaTDen$AgeOri)

# data.scores = as.data.frame(scores(BaTDen_dis.mds)) # extract the site scores and convert to a dataframe
# data.scores$site = rownames(data.scores) # create a column of site names
# data.scores$group = BaTDen$Group # add the Age Group Variable
# head(data.scores) #look at the data
# # 
# species.scores = as.data.frame(scores(BaTDen_dis.mds, "species")) # extract the species scores and convert to dataframe
# species.scores$species = rownames(species.scores) # create a column of species
# head(species.scores) #look at the data

# data.scores = as.data.frame(scores(BaTDen_dis.mds, display = c("sites", "species")))
# data.scores
# 
# #vec.sp = envfit(BaTDen_dis.mds$points, BaTDen_dis, permutations = 1000)
# vec.sp = envfit(BaTDen_dis.mds$points, BaTDen_dis)
# #vec.sp.df = as.data.frame(vec.sp$vector$arrows*sqrt(vec.sp$vectors$r))
# vec.sp.df = as.data.frame(vec.sp$vector$arrows*1)
# vec.sp.df
# #vec.sp.df$species = rownames(vec.spdf)
# vec.sp$species = rownames(vec.sp)
# vec.sp



# plot the NMDS
ggplot() + 
  geom_point(NMDS_Den, aes(x=MDS1, y=MDS2, color = Group)) + 
  #geom_point(aes(data = NMDS_Den, color = NMDS_Den$Group)) +
  #geom_segment(data = vec.sp.df, aes(x=0, xend=MDS1, y=0, yend=MDS2),
   #            arrow = arrow(length = unit(0.5, "cm")), colour = "grey", inherit.aes = FALSE) +
  geom_text(data = vec.sp, aes(x=MDS1, y=MDS2, label = species)) + 
  #geom_text(data = data.scores, aes(x=NMDS1, y=NMDS2, label = species), cex=0.8, col= "black") +
  theme_bw() + #plot points without ellipse
  labs(x="NMDS1", y="NMDS2")+
  theme(panel.background = element_blank(), # remove panel background colour
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


ggplot(NMDS_Den, aes(x=MDS1, y=MDS2, col=Group)) + stat_ellipse() + theme_bw() + #ellipses of age groups
  labs(x="NMDS1", y="NMDS2")+
  theme(panel.background = element_blank(), # remove panel background colour
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


ggplot(NMDS_Den, aes(x=MDS1, y=MDS2, col=Origin)) + geom_point() + theme_bw() +
  labs(x="NMDS1", y="NMDS2")+
  theme(panel.background = element_blank(), # remove panel background colour
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot(NMDS_Den, aes(x=MDS1, y=MDS2, col=Origin)) + stat_ellipse() + theme_bw() +
  labs(x="NMDS1", y="NMDS2")+
  theme(panel.background = element_blank(), # remove panel background colour
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

anosim_Den = anosim(BaTDen_dis, BaTDen$Group)
#anosim_Den2 = anosim(BaTDen_dis, BaTDen$Origin)
anosim_Den # view the results
summary(anosim_Den)
#summary(anosim_Den2)
plot(anosim_Den)
#plot(anosim_Den2)


# BaTDen_f = subset(BaTDen, Origin=="fCAFS") # subset plots created from forest
# BaTDen_s = subset(BaTDen, Origin=="sCAFS") # subset plots created from savannah
# 

# --------------------------------------------------------------------------

library(vegan)
library(ggplot2)

BacoaAss = read.csv(file = "D:...../BTDens_nmds.csv", 
                    header = TRUE, check.names = FALSE, row.names = 1)

# Reorder the levels of the factor variable
BacoaAss$Group = ordered(BacoaAss$Group, 
                       levels = c("<= 10yrs","11 - 20yrs","21 - 40yrs","41 - 60yrs","  > 60yrs  "))

#summary(BacoaAss)
#View(BacoaAss)
names(BacoaAss)

#names(BaTDis)

BacoaAss$Group = as.factor(BacoaAss$Group)
BacoaAss$Origin = as.factor(BacoaAss$Origin)
BacoaAss$AgeOri = as.factor(BacoaAss$AgeOri)


BacoaAss_c = BacoaAss[-11][5:11] # subset just parameters for structural composition of farms
#BaTDis_c = BaTDis[6:15]
names(BacoaAss_c)

# NMDS of Plot parameters
#BacoaAss_c.mds = metaMDS(comm = BacoaAss_c, distance = "bray", autotransform = FALSE)
BacoaAss_c.mds = metaMDS(comm = BacoaAss_c, distance = "bray", k = 2, autotransform = FALSE)
BacoaAss_c.mds$stress
plot(BacoaAss_c.mds$points)
plot(BacoaAss_c.mds, type = "t")


colgroup = c("lightcoral", "gold3", "springgreen3", "skyblue1", "violet")
colorigin = c("lightcoral", "skyblue1")

with(BacoaAss, levels(Group)) # load data and use the Group variable
with(BacoaAss, levels(Origin)) # load data and use the Origin
scl = 3 # scalling = 3


data.scores = as.data.frame(scores(BacoaAss_c.mds)) # extract the site scores and convert to a dataframe
data.scores$site = rownames(data.scores) # create a column of site names
data.scores$group = BacoaAss$Group # add the Age Group Variable
data.scores$origin = BacoaAss$Origin # add the Origin Variable
head(data.scores) #look at the data

species.scores = as.data.frame(scores(BacoaAss_c.mds, "species")) # extract the species scores and convert to dataframe
species.scores$species = rownames(species.scores) # create a column of species
head(species.scores) #look at the data


# plot the NMDS by Farm Age
ggplot() + 
  geom_point(data = data.scores, aes(x=NMDS1, y=NMDS2, colour = group), size=2) +
  geom_text(data = species.scores, aes(x=NMDS1, y=NMDS2, label=species), cex=3, col= "black") + 
  #geom_text(data = data.scores, aes(x=NMDS1, y=NMDS2, label = site)) + 
  #geom_text(data = data.scores, aes(x=NMDS1, y=NMDS2, label = species), cex=0.8, col= "black") +
  #coord_equal() +
  theme_bw() + #plot points without ellipse
  labs(x="NMDS1", y="NMDS2")+
  theme(panel.background = element_blank(), # remove panel background colour
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


# compute centroids (means) for each age group and 95% CI ellipse 
library(ellipse)
centroid_g = aggregate(cbind(NMDS1,NMDS2)~ group,data.scores, mean)
conf_g.rgn = do.call(rbind,lapply(unique(data.scores$group), function(t)
  data.frame(group=as.character(t),
             ellipse(cov(data.scores[data.scores$group==t,1:2]),
                     centre=as.matrix(centroid_g[t,2:3]),
                     level=0.95),
             stringsAsFactors = FALSE)))

ggplot(data = data.scores, aes(x=NMDS1, y=NMDS2, colour = group)) + 
  #geom_path(data = conf_g.rgn, cex=0.8, linetype=2)+ 
  geom_path(data = conf_g.rgn, cex=0.8)+
  theme_bw() + #plot points without ellipse
  labs(x="NMDS1", y="NMDS2")+
  theme(panel.background = element_blank(), # remove panel background colour
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



# plot the NMDS by Farm oring
ggplot() + 
  geom_point(data = data.scores, aes(x=NMDS1, y=NMDS2, colour = origin), size=2) +
  geom_text(data = species.scores, aes(x=NMDS1, y=NMDS2, label=species), cex=3, col= "black") + 
  #geom_text(data = data.scores, aes(x=NMDS1, y=NMDS2, label = site)) + 
  #geom_text(data = data.scores, aes(x=NMDS1, y=NMDS2, label = species), cex=0.8, col= "black") +
  #coord_equal() +
  theme_bw() + #plot points without ellipse
  labs(x="NMDS1", y="NMDS2")+
  theme(panel.background = element_blank(), # remove panel background colour
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


# compute centroids (means) for farm type and 95% CI ellipse 

centroid_o = aggregate(cbind(NMDS1,NMDS2)~ origin,data.scores, mean)
conf_o.rgn = do.call(rbind,lapply(unique(data.scores$origin), function(t)
  data.frame(origin=as.character(t),
             ellipse(cov(data.scores[data.scores$origin==t,1:2]),
                     centre=as.matrix(centroid_o[t,2:3]),
                     level=0.95),
             stringsAsFactors = FALSE)))

ggplot(data = data.scores, aes(x=NMDS1, y=NMDS2, colour = origin)) + 
  #geom_path(data = conf_o.rgn, cex=0.8, linetype=2)+ 
  geom_path(data = conf_o.rgn, cex=0.8 )+ 
  theme_bw() + #plot points without ellipse
  labs(x="NMDS1", y="NMDS2")+
  theme(panel.background = element_blank(), # remove panel background colour
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


# 
# 
# 
# plot(BacoaAss_c.mds, type = "n")
# with(BacoaAss, points(BacoaAss_c.mds, display="sites", cex=1.2, col=colorigin[Origin], # Plot Origin
#                       pch=21, bg = colorigin[Origin]))
# text(BacoaAss_c.mds, display = "species", cex=0.8, col= "black")
# 
# 
# plot(BacoaAss_c.mds, type = "n")
# #points(BaTDen_c.mds, display = "sites", col=colvec[BaTDen$Group])
# text(BacoaAss_c.mds, display = "sites", cex=0.8, col=colgroup[BacoaAss$Group]) # For Plot Age
# text(BacoaAss_c.mds, display = "species", cex=0.8, col= "black")
# 
# plot(BacoaAss_c.mds, type = "n")
# #points(BaTDen_c.mds, display = "sites", col=colvec[BaTDen$Group])
# text(BacoaAss_c.mds, display = "sites", cex=0.8, col=colorigin[BacoaAss$Origin]) # For Plot Origin
# text(BacoaAss_c.mds, display = "species", cex=0.8, col= "black")
# 
# 
# BacMDS_xy = data.frame(BacoaAss_c.mds$points)
# BacMDS_xy$Group = BacoaAss$Group
# BacMDS_xy$Origin = BacoaAss$Origin
# BacMDS_xy$AgeOri = BacoaAss$AgeOri
# BacMDS_xy$nAsT = BacoaAss$nAsT
# 
# 
# 
# ggplot(BacMDS_xy, aes(MDS1, MDS2, color = Group)) + geom_point() + stat_ellipse() + theme_bw()
# 
# anosim_group = anosim(BacoaAss_c.mds, BacMDS_xy$Group)
# anosim_group # take a look at the result
# summary(anosim_group)
# plot(anosim_group)
# 











#-------------------------------------------------------
#-------------------------------------------------------
# B) NMDS of Composition of Live Woody Biomass


BacoaAss_b = BacoaAss[-11][12:15] # subset just parameters for structural composition of farms
#BaTDis_b = BaTDis[6:15]
names(BacoaAss_b)

# NMDS of Plot parameters
#BacoaAss_b.mds = metaMDS(comm = BacoaAss_b, distance = "bray", autotransform = FALSE)
BacoaAss_b.mds = metaMDS(comm = BacoaAss_b, distance = "bray", k = 2, autotransform = FALSE)
BacoaAss_b.mds$stress
plot(BacoaAss_b.mds$points)
plot(BacoaAss_b.mds, type = "t")


data.scores = as.data.frame(scores(BacoaAss_b.mds)) # extract the site scores and convert to a dataframe
data.scores$site = rownames(data.scores) # create a column of site names
data.scores$group = BacoaAss$Group # add the Age Group Variable
data.scores$origin = BacoaAss$Origin # add the Origin Variable
head(data.scores) #look at the data

species.scores = as.data.frame(scores(BacoaAss_b.mds, "species")) # extract the species scores and convert to dataframe
species.scores$species = rownames(species.scores) # create a column of species
head(species.scores) #look at the data


# plot the NMDS by Farm Age
ggplot() + 
  geom_point(data = data.scores, aes(x=NMDS1, y=NMDS2, colour = group), size=2) +
  geom_text(data = species.scores, aes(x=NMDS1, y=NMDS2, label=species), cex=3, col= "black") + 
  #geom_text(data = data.scores, aes(x=NMDS1, y=NMDS2, label = site)) + 
  #geom_text(data = data.scores, aes(x=NMDS1, y=NMDS2, label = species), cex=0.8, col= "black") +
  #coord_equal() +
  theme_bw() + #plot points without ellipse
  labs(x="NMDS1", y="NMDS2")+
  theme(panel.background = element_blank(), # remove panel background colour
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


# compute centroids (means) for each age group and 95% CI ellipse 
library(ellipse)
centroid_b_g = aggregate(cbind(NMDS1,NMDS2)~ group,data.scores, mean)
conf_b_g.rgn = do.call(rbind,lapply(unique(data.scores$group), function(t)
  data.frame(group=as.character(t),
             ellipse(cov(data.scores[data.scores$group==t,1:2]),
                     centre=as.matrix(centroid_b_g[t,2:3]),
                     level=0.95),
             stringsAsFactors = FALSE)))

ggplot(data = data.scores, aes(x=NMDS1, y=NMDS2, colour = group)) + 
  #geom_path(data = conf_b_g.rgn, cex=0.8, linetype=2)+ 
  geom_path(data = conf_b_g.rgn, cex=0.8)+
  theme_bw() + #plot points without ellipse
  labs(x="NMDS1", y="NMDS2")+
  theme(panel.background = element_blank(), # remove panel background colour
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



# plot the NMDS by Farm oring
ggplot() + 
  geom_point(data = data.scores, aes(x=NMDS1, y=NMDS2, colour = origin), size=2) +
  geom_text(data = species.scores, aes(x=NMDS1, y=NMDS2, label=species), cex=3, col= "black") + 
  #geom_text(data = data.scores, aes(x=NMDS1, y=NMDS2, label = site)) + 
  #geom_text(data = data.scores, aes(x=NMDS1, y=NMDS2, label = species), cex=0.8, col= "black") +
  #coord_equal() +
  theme_bw() + #plot points without ellipse
  labs(x="NMDS1", y="NMDS2")+
  theme(panel.background = element_blank(), # remove panel background colour
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


# compute centroids (means) for farm type and 95% CI ellipse 

centroid_b_o = aggregate(cbind(NMDS1,NMDS2)~ origin,data.scores, mean)
conf_b_o.rgn = do.call(rbind,lapply(unique(data.scores$origin), function(t)
  data.frame(origin=as.character(t),
             ellipse(cov(data.scores[data.scores$origin==t,1:2]),
                     centre=as.matrix(centroid_b_o[t,2:3]),
                     level=0.95),
             stringsAsFactors = FALSE)))

ggplot(data = data.scores, aes(x=NMDS1, y=NMDS2, colour = origin)) + 
  #geom_path(data = conf_b_o.rgn, cex=0.8, linetype=2)+ 
  geom_path(data = conf_b_o.rgn, cex=0.8 )+ 
  theme_bw() + #plot points without ellipse
  labs(x="NMDS1", y="NMDS2")+
  theme(panel.background = element_blank(), # remove panel background colour
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


