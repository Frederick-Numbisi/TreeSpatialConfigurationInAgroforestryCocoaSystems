
library(sp)
library(spatstat)
library(maptools)
library(ggplot2)

setwd("D:/...../...../Bakoa_TreesClass4") # Select location of dataset as the working directory

# -----------------------------------------
# Analysis procedure for data in each plot
# -----------------------------------------
# Sample code for plot

Plot23_Data = read.csv("Plot23_tclass4.csv") # Read the data sete
Plot23_Data$TTP = as.factor(Plot23_Data$TTP) # Convert the TTP variable to factor 

Plot23_T4  <- readShapePoints("Plot23_tclass4.shp") # Read the spatial data (shapefile) of tree point locations
SP <- as(Plot23_T4, "SpatialPoints") # Convert the spatial data to spatial point patterns
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot23_Data[5]

# Get Plot231 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot23.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot23 windown (polygon) to trees in Plot23
P$window <- W
plot(P, pch=16, main="Plot23 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot23 All Trees")
plot(density(P))
#legend(0.25, 0.5, names(a), pch=a)

#Set measurement unit of window()
#unitname(P) <- "meter"

# The Kdot(r) inhomogenous function for each tree group within stand

# method (1): estimate intensities by nonparametric smoothing
Co <- split(P)$Co
AsT <- split(P)$AsT
All = unmark(P)
lambdaAll = density.ppp(All, sigma=0.15, at ="points")
lambdaCo <- density.ppp(Co, sigma=0.15, at="points")
lambdaAsT <- density.ppp(AsT, sigma=0.15, at="points")


# Fit model (second order quadratic distribution)
fit <- ppm(P ~ marks * polynom(x,y,2))
L_AsT_plot23 <- Ldot.inhom(P, "AsT", lambdaX=fit, update=FALSE)
plot(Lij_Plot23, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot23 Associated Trees", legendpos="bottom")
L_AsT_plot23$iso  <- L_AsT_plot23$iso  - L_AsT_plot23$r
L_AsT_plot23$theo <- L_AsT_plot23$theo - L_AsT_plot23$r
write.csv(L_AsT_plot23, file = "Li_AsT_Plot23.csv", row.names=TRUE)

#Li (nonhomogenous) for Cocoa Trees2 in stand 

L_Co2_plot23 <- Ldot.inhom(P, "Co", lambdaX=fit, update=FALSE) #using fit computed for all tree

plot(Lij_Plot23, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot23 Cocoa Trees3", legendpos="bottom")
L_AsT_plot23$iso  <- L_AsT_plot23$iso  - L_AsT_plot23$r
L_AsT_plot23$theo <- L_AsT_plot23$theo - L_AsT_plot23$r
write.csv(L_AsT_plot23, file = "Li_Co3_Plot23.csv", row.names=TRUE)


#Bivariate(cross) inhomogenous Kij(r) Ripley function
#The Kcross(r) test, which examines distribution a tree group relative to another group within stand

# method (3): fit parametric intensity model
fit <- ppm(P, ~marks * polynom(x,y,2))
# evaluate fitted intensities at data points
# (these are the intensities of the sub-processes of each type)
inten <- fitted(fit, dataonly=TRUE)
# split according to types of points
lambda2 <- split(inten, P$marks)
Kij_plot23 <- Kcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))
plot(Kij_Plot23, xlab="Distance (m)", main="Kij (2nd order moment) Co3 to Ast")

write.csv(Kij_Plot23, file = "Kij_bivariate_inhom_Plot23_Co3-AsT.csv", row.names=TRUE)



#Bivariate Lij(r) stabilized distance variance of inhomogenous Ripley function

Co <- split(P)$Co
AsT <- split(P)$AsT

# method (1): estimate intensities by nonparametric smoothing
#lambdaCo <- density.ppp(Co, sigma=0.15, at="points")
#lambdaAsT <- density.ppp(AsT, sigma=0.15, at="points")
#L <- Lcross.inhom(P, "Co", "AsT", lambdaCo, lambdaAsT)

# method (2): fit parametric intensity model
fit <- ppm(P, ~marks * polynom(x,y,2))
# evaluate fitted intensities at data points
# (these are the intensities of the sub-processes of each type)
inten <- fitted(fit, dataonly=TRUE)
# split intensity according to types of points
lambda2 <- split(inten, P$marks)
Lij_Plot23 <- Lcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))

plot(Lij_Plot23, . - r ~ r, xlab="Distance (m)", ylab="Lij(r)", main="Lij inhom. Function, Plot23 Cocoa3-AsT", legendpos="bottom")


# We will need to modify the plot graphics parameters to accommodate the
# custom y-label; hence the call to `par`.
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Plot23, . -r ~ r, ylab=expression(hat("Lij")), xlab = "d (m)", main="Lij(r), Plot23 Cocoa Trees", legendmath=TRUE, legendpos="bottom")
#par(OP)

#We can make permanent changes to the Lij object's values as follows:
Lij_Plot23$iso  <- Lij_Plot23$iso  - Lij_Plot23$r
Lij_Plot23$theo <- Lij_Plot23$theo - Lij_Plot23$r

write.csv(Lij_Plot23, file = "Lij_Bivariate_inhom_Plot23_Co3_to_AsT.csv", row.names=TRUE)

#MONTE CARLO (MC) SIMULATION

#To test for clustering/inhibition in the presence of spatial inhomogeneity, 
#the null hypothesis should be an inhomogeneous Poisson process. 
#So we estimate the inhomogeneous intensity functions of the two types of points,
# and then generate simulated point patterns according to these intensities.

#First fit a Poisson point process model to the observed data
#fit2 <- ppm(P, ~marks * polynom(x,y,3))

#n <- 199 #number of iterations for Monte carlo simulation
#p  <- 0.01 # Desired p significance level to display

#Lij_Env_Plot23 <- envelope(fit2, Lcross.inhom, lambdaX=fit2, nsim=n, rank=(p * (n + 1)))
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Env_Plot23, . - r ~ r, ylab=expression(hat("L")),xlab="Distance (m)", main=paste("Lij(r) Plot23 Co3-AsT,  p= ", p),legendmath=TRUE)
#par(OP)
#write.csv(Lij_Env_Plot23, file = "Lij_199MCEnv_Plot23_Co3_to_AsT.csv", row.names=TRUE)



# COMPUTING THE g (Pair Correlation Fucntion) OF COCOA AND ASSOCIATED TREES
#The function value (fv) table method:we will use computed inhomogenous cross K function

# fv object
#Kij_plot23 is the fv object computed above (inhomogenous cross K from Co to AsT)
#method "b" apply smoothing to Y(r) = K(r)/(2 * pi * r) constraining Y(0) = 0, 
#estimate the derivative of Y, and solve;
#method "b" used due to likely inhibition at small distances of r spar=0.8,

gij_plot23 <- pcf(Kij_plot23, spar=0.8,  method="b") 
plot(gij_plot23, main="Plot23 Co3 - AsT g(r) Pair correlation function", xlab="Distance(m)", legendpos="bottom")
write.csv(gij_plot23, file = "pcf_Plot23_Co3_to_AsT.csv", row.names=TRUE)



#-----------------------------------------------
#
# SAVING AND APPENDING THE Ripley-K RESULTS AND FILES IN CSV
#--------------------------------------------------

#A) Appending Multivariate Lij(r) Results

Lij_Results <- c(Lij_Plot2, Lij_Plot3, Lij_Plot4, Lij_Plot5, Lij_Plot6, Lij_Plot8,
                 Lij_Plot10, Lij_Plot11, Lij_Plot12, Lij_Plot13, Lij_Plot21, 
                 Lij_Plot22, Lij_Plot23, Lij_Plot26, Lij_Plot28, Lij_Plot29,
                 Lij_Plot34, Lij_Plot38, Lij_Plot42, Lij_Plot43, Lij_Plot44, Lij_Plot50)


write.csv(Lij_Results, "Lij_Results_Co3-AsT.csv", row.names = TRUE)


#B) Appending pcf-g(r) results

PCFij_Results <- c(gij_Plot2, gij_Plot3, gij_Plot4, gij_Plot5, gij_Plot6, gij_Plot8,
                   gij_Plot10, gij_Plot11, gij_Plot12, gij_Plot13, gij_Plot21,  
                   gij_Plot22, gij_Plot23, gij_Plot26, gij_Plot28, gij_Plot29, 
                   gij_Plot34, pc_Plot38, pc_Plot42, pc_Plot43, pc_Plot44, pc_Plot50)


write.csv(PCFij_Results, "PCFij_Results_Co3-AsT.csv", row.names = TRUE)



#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
#
# CREATIONG OF HYPERFRAME FOR COAF1 AND ASSOCIATED TREES
#
#-------------------------------------------------------------------------

library(sp)
library(spatstat)
library(maptools)
library(ggplot2)

setwd("D:/PhD Manuscript/...../Bakoa_TreesClass4")

# The analysis procedure is replicated for all plots using a function call
# Case analysis for Plot23

Plot23_Data = read.csv("Plot23_tclass4.csv")
Plot23_Data$TTP = as.factor(Plot23_Data$TTP)

Plot23_T4  <- readShapePoints("Plot23_tclass4.shp")
SP <- as(Plot23_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot23_Data[5]

# Get Plot23 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot23.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot23 window (polygon) to trees in Plot23
P$window <- W
plot(P, pch=16, main="Plot23 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot23 All Trees")
plot(density(P))

P23T4 = P
summary(P23T4)

P23T4Co = unmark(split(P23T4)$Co) # To subset just the cocoa trees 
P23TAsT = unmark(split(P23T4)$AsT) # To subset jut Associated trees
P23T4All = unmark(P23T4) # Remove marks (tree type) from point pattern data



#---------------------------------------------------------
#---------------------------------------------------------
# Compute Hyperframe dataset for COAF1 - Assoiated Trees
#
#---------------------------------------------------------

## Create list of spatial objects from spatial point patterns cocoa and associated trees
T4Al = solist(P23T4All,P23bT4All,P23cT4All,P26T4All,P26AaT4All,P26bT4All,P26cT4All,P28T4All,P28AaT4All,P28bT4All,P28cT4All,P38T4All,
              P4T4All,P22T4All,P22AaT4All,P22bT4All,P22cT4All,P6T4All,P6bT4All,P43T4All,P5T4All,P5bT4All,P21T4All,P21AaT4All,
              P21bT4All,P21cT4All,P8T4All,P8bT4All,P34T4All,P34bT4All,P34cT4All,P29T4All,P29bT4All,P42T4All,P42bT4All,P3T4All,
              P3bT4All,P10T4All,P10bT4All,P11T4All,P11bT4All,P12T4All,P2T4All,P2bT4All,P13T4All,P13bT4All,P39T4All,P39bT4All,
              P44T4All,P44bT4All,P44cT4All,P50T4All,P50AaT4All,P50bT4All,P50cT4All,P45aT4All,P48bT4All,P49bT4All,P60aT4All,PNS1T4All,POS3T4All)

# Create list of spatial objects from spatial point patterns for cocoa tree 
T4Co = solist(P23T4Co,P23bT4Co,P23cT4Co,P26T4Co,P26AaT4Co,P26bT4Co,P26cT4Co,P28T4Co,P28AaT4Co,P28bT4Co,P28cT4Co,P38T4Co,
              P4T4Co,P22T4Co,P22AaT4Co,P22bT4Co,P22cT4Co,P6T4Co,P6bT4Co,P43T4Co,P5T4Co,P5bT4Co,P21T4Co,P21AaT4Co,
              P21bT4Co,P21cT4Co,P8T4Co,P8bT4Co,P34T4Co,P34bT4Co,P34cT4Co,P29T4Co,P29bT4Co,P42T4Co,P42bT4Co,P3T4Co,
              P3bT4Co,P10T4Co,P10bT4Co,P11T4Co,P11bT4Co,P12T4Co,P2T4Co,P2bT4Co,P13T4Co,P13bT4Co,P39T4Co,P39bT4Co,
              P44T4Co,P44bT4Co,P44cT4Co,P50T4Co,P50AaT4Co,P50bT4Co,P50cT4Co,P45aT4Co,P48bT4Co,P49bT4Co,P60aT4Co,PNS1T4Co,POS3T4Co)

# Create list of spatial objects from spatial point patterns for Associated trees 
TAsT = solist(P23TAsT,P23bTAsT,P23cTAsT,P26TAsT,P26AaTAsT,P26bTAsT,P26cTAsT,P28TAsT,P28AaTAsT,P28bTAsT,P28cTAsT,P38TAsT,
              P4TAsT,P22TAsT,P22AaTAsT,P22bTAsT,P22cTAsT,P6TAsT,P6bTAsT,P43TAsT,P5TAsT,P5bTAsT,P21TAsT,P21AaTAsT,
              P21bTAsT,P21cTAsT,P8TAsT,P8bTAsT,P34TAsT,P34bTAsT,P34cTAsT,P29TAsT,P29bTAsT,P42TAsT,P42bTAsT,P3TAsT,
              P3bTAsT,P10TAsT,P10bTAsT,P11TAsT,P11bTAsT,P12TAsT,P2TAsT,P2bTAsT,P13TAsT,P13bTAsT,P39TAsT,P39bTAsT,
              P44TAsT,P44bTAsT,P44cTAsT,P50TAsT,P50AaTAsT,P50bTAsT,P50cTAsT,P45aTAsT,P48bTAsT,P49bTAsT,P60aTAsT,PNS1TAsT,POS3TAsT)


sapply(T4Co, npoints) #to get number of Cocoa trees (points) for each spatial pattern (sampled plot)
sapply(TAsT, npoints) #to get number of Associated trees (points) for each spatial pattern (sampled plot)

# CoafT0_D = c(45,51,58,62,67,55,64,74,49,54,57,30,61,41,56,57,38,
#             66,77,53,51,68,46,10,35,44,56,53,62,77,65,50,61,74,
#             125,70,60,30,43,46,59,34,41,61,27,55,46,66,65,83,55,
#             54,50,37,39,75,70,72,52,59,63)
#
# Create list of ID for plots
Plots = as.factor (c('23A','23B','23C','26A','26Aa','26B','26C','28A','28Aa','28B','28C','38A',
                     '4A','22A','22Aa','22B','22C','6A','6B','43A','5A','5B','21A','21Aa',
                     '21B','21C','8A','8B','34A','34B','34C','29A','29B','42A','42B','3A',
                     '3B','10A','10B','11A','11B','12A','2A','2B','13A','13B','39A','39B',
                     '44A','44B','44C','50A','50Aa','50B','50C','45A','48B','49b','60A','NS1','OS3'))
# Create list of code for plots
PID = c(1,2,3,4,5,6,7,8,9,10,11,12,
        13,14,15,16,17,18,19,20,21,22,23,24,
        25,26,27,28,29,30,31,32,33,34,35,36,
        37,38,39,40,41,42,43,44,45,46,47,48,
        49,50,51,52,53,54,55,56,57,58,59,60,61)

AgeG = as.factor(c(1,1,1,2,2,2,2,2,2,2,2,2,
                   3,3,3,3,3,3,3,4,4,4,5,5,
                   5,5,5,5,1,1,1,2,2,2,2,3,
                   3,3,3,4,4,4,4,4,5,5,5,5,
                   5,5,5,5,5,5,5,2,3,4,1,1,5))

OrigP = as.factor(c("Fo","Fo","Fo","Fo","Fo","Fo","Fo","Fo","Fo","Fo","Fo","Fo",
                    "Fo","Fo","Fo","Fo","Fo","Fo","Fo","Fo","Fo","Fo","Fo","Fo",
                    "Fo","Fo","Fo","Fo","Sa","Sa","Sa","Sa","Sa","Sa","Sa","Sa",
                    "Sa","Sa","Sa","Sa","Sa","Sa","Sa","Sa","Sa","Sa","Sa","Sa",
                    "Sa","Sa","Sa","Sa","Sa","Sa","Sa","Sa","Fo","Fo","Fo","Fo","Fo"))

#PID = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)
#AgeG = as.factor(c(1,2,2,2,3,3,3,4,4,5,5,1,2,2,3,3,4,4,4,5,5,5,5))
#OrigP = as.factor(c("Fo","Fo","Fo","Fo","Fo","Fo","Fo","Fo","Fo","Fo","Fo","Sa","Sa","Sa","Sa","Sa","Sa","Sa","Sa","Sa","Sa","Sa","Sa"))

#Create Hyperframe Dataset

BCoa4 = hyperframe(Coa4 = T4Co, AsT = TAsT, Plots = Plots, Code=(PID), Group= AgeG, Origin = OrigP)
BCoa4
plot(BCoa4)
BCoa4$dCoa4 = with(BCoa4, distfun(Coa4)) # Estimate distance function (map) column for cocoa trees
BCoa4$DistMc4 = with(BCoa4, distmap(Coa4)) # Create distance map column for Cocoa trees

BCoa4$dAsT = with(BCoa4, distfun(AsT)) # Estimate distance function (map) column for AsT Trees
BCoa4$DistMAsT = with(BCoa4, distmap(AsT)) # Create distance map column for Associated Trees

BCoa4$nAsT = with(BCoa4, npoints(AsT)) # Create column for No. of AsT trees
BCoa4$nCoa4 = with(BCoa4, npoints(Coa4)) # Create column for No. of Cocoa Trees

# Find the minimum interpoint distance in each point pattern (plot)

BCoa4$MinDC4 = with(BCoa4, min(nndist(Coa4))) # Create column of min. distance between Cocoa trees in each plot

BCoa4$MaxDC4 = with(BCoa4, max(nndist(Coa4))) # Create column of max. distance between Cocoa trees in each plot

BCoa4$MinDAsT = with(BCoa4, min(nndist(AsT)))# Create column of min. distance between AsT trees in each plot
BCoa4$MaxDAsT = with(BCoa4, max(nndist(AsT))) # Create column of max. distance between AsT trees in each plot


#EXPLORATORY ANALYSIS

BCoa4$K4 = with(BCoa4, Kest(Coa4)) # Create column of K-Function for each point pattern

BCoa4$L4 = with(BCoa4, Lest(Coa4)) # Create column of L-Function for each point pattern

BCoa4$P4 = with(BCoa4, pcf(Coa4))  # Create column of Pair Correlation Function (Generates warning for plots with few points/trees)

BCoa4$G4 = with(BCoa4, Gest(Coa4)) # Create the cumulative nearest neighbour distance (G-Function)


BCoa4$area = with(BCoa4, area(Coa4)) # create colummn of surface area for each point pattern



save(BCoa4, file = "D:/....../BCoa4.Rda")



# ----------------------------------------------------
# -----------------------------------------------------

# PLOTTING THE POOLED L-FUNCTION FOR SPATIAL POINT PATTERNS


setwd("D:/...../BakoaHyperframes")
load('BCoa4.Rda')


names(BCoa4)
py4 = BCoa4
py4$area = with(py4, area(Coa4)) # create colummn of area for each point pattern
py4 = as.data.frame(py4, warn = FALSE) # convert hyperframe to dataframe
is.na(py4) <- sapply(py4, is.infinite)
py4


library("dplyr") # the pplyr package very handy for data selection, filtering, summary etc.


write.csv(py4, file = "D:/...../BCoa4.csv")


py4 <- within(py4, {
  Group <- factor(Group, levels=1:5, labels=c("<10 yrs", "10-20 yrs", "21-40 yrs",
                                              "41-60 yrs", ">60 yrs"))
  Origin <- factor(Origin)
})

par(mar = c(5, 4, 4, 4) + 0.3)  # To provide space for the z axis
#plot(sqrt(nCoa4) ~ interaction(Origin,Group), data=py4, xlab="Age Group", main="Cocoa Trees in sampled Plots", legend=TRUE) # Box plot of sqrt (no. points) by age group
plot(sqrt(nCoa4) ~ Group, data=py4, xlab="Age Group", main="Cocoa Trees in sampled Plots") # Box plot of sqrt (no. points) by age group
v = pretty(py4$nCoa4);axis(4, at=sqrt(v),labels=v, ylab="No. Coa4")
mtext("No. of Trees", side=4, line=3)




sapply(split(py4$nCoa4/py4$area, py4$Group), mean)

fitn = glm(nCoa4 ~ offset(log(area)) + Group, family=poisson, data=py4)
newd = data.frame(area=1, Group=levels(py4$Group))
predict(fitn, newdata=newd, type="response")


#boxplot(Sal ~ interaction(group,class), data=test.data)  


# POOLED SUMMARY PLOT OF DISTANCE FUNCTIONS
#-----------------------------------------------------------

py4

Lsplit = split(py4$L4, py4$Group) # Split files based on transformed K-Funtion (L Column)
Lpool = anylapply(Lsplit, pool)
summary(Lpool)
plot(Lpool, cbind(pooliso,pooltheo,hiiso,loiso) - radius  ~ radius , legend=FALSE,
     shade=c("hiiso", "loiso"), xlim=c(0, 6), xlab="radius (m)", equal.scales=TRUE)
plot(Lpool, cbind(pooliso,pooltheo,hiiso,loiso) - radius  ~ radius , legend=T,
     shade=c("hiiso", "loiso"), xlim=c(0, 6), xlab="radius (m)", equal.scales=TRUE)
plot(Lpool[1], cbind(pooliso,pooltheo,hiiso,loiso) - radius  ~ radius , main="Plots < 11yrs",
     shade=c("hiiso", "loiso"), xlim=c(0, 6), xlab="radius (m)", equal.scales=TRUE)
n <- 39 #number of iterations for Monte carlo simulation
# Use MC function called envelope which can be used to compute multiple iterations of the Lest function.
p  <- 0.05 # Desired p significance level to display
py4$Lenv <- with(py4, envelope(Coa4, Lest, nsim=n, rank=(p * (n + 1)), savefuns=TRUE))
EvenL = py4$Lenv


Lenvsplit = pool(py$Lenv, py$Group) # Split files based on transformed K-Funtion (L Column)
Lenvpool = anylapply(Lenvsplit, pool)

summary(Lenvsplit)
Kall = pool(as.anylist(Lenvsplit))
Lenvpool1 = pool(Lenvsplit[[1]])

plot(Lenvpool1)

plot(Lenv_Plot5, . - r ~ r, ylab=expression(hat("L(r)")),xlab="Distance (m)", main=paste("Plot5 nsim=999,  p= ", p),legendmath=TRUE)


plot(Lpool, cbind(pooliso,pooltheo,hiiso,loiso) - r ~ r, legend=FALSE,
     shade=c("hiiso", "loiso"), xlim=c(0, 6), xlab="r(m)", equal.scales=TRUE)

# Plotting individual estimates of the L-function for each plot,
# superimposed on each other within each group of plots (age group)
LeachG = anylapply(Lsplit, collapse.fv, same="theo", different="iso")
plot(LeachG, legend=FALSE, xlim=c(0, 6), ylim=c(0, 6), xlab="radius (m)")


plot(Lpool[1], cbind(pooliso,pooltheo,hiiso,loiso) - r  ~ r , main="Plots < 10yrs",
     shade=c("hiiso", "loiso"), xlim=c(0, 6), xlab="radius (m)", equal.scales=TRUE)

plot(Lpool[2], cbind(pooliso,pooltheo,hiiso,loiso) - r  ~ r , main="Plots 10-20yrs",
     shade=c("hiiso", "loiso"), xlim=c(0, 6), xlab="radius (m)", equal.scales=TRUE, 
     legendpos="float")
plot(Lpool[3], cbind(pooliso,pooltheo,hiiso,loiso) - r  ~ r , main="Plots 21-40yrs",
     shade=c("hiiso", "loiso"), xlim=c(0, 6), xlab="radius (m)", equal.scales=TRUE, 
     legendpos="float")
plot(Lpool[4], cbind(pooliso,pooltheo,hiiso,loiso) - r  ~ r , main="Plots 41-60yrs",
     shade=c("hiiso", "loiso"), xlim=c(0, 6), xlab="radius (m)", equal.scales=TRUE, 
     legendpos="float")
plot(Lpool[5], cbind(pooliso,pooltheo,hiiso,loiso) - r  ~ r , main="Plots > 60yrs",
     shade=c("hiiso", "loiso"), xlim=c(0, 6), xlab="radius (m)", equal.scales=TRUE, 
     legendpos="float")


py4 = BCoa4

names(py4)


Lsplit = split(py4$L4, py4$Group) # Split files based on transformed K-Funtion (L4 Column)
Lpool = anylapply(Lsplit, pool)
summary(Lpool)
Psplit = split(py4$P4, py4$Group) # Split files based on Pair Correlation Function (P4 Column)
Ppool = anylapply(Psplit, pool)

plot(Ppool ,
     shade=c("hiiso", "loiso"), xlim=c(0, 6), xlab="radius (m)", equal.scales=TRUE)


plot(Ppool[1],  main="Plots < 11yrs",
     xlim=c(0, 6), xlab="radius (m)", equal.scales=TRUE)

plot(Ppool[2],  main="Plots 11-20yrs",
     xlim=c(0, 6), xlab="radius (m)", equal.scales=TRUE)

plot(Ppool[3],  main="Plots 21-40yrs",
     xlim=c(0, 6), xlab="radius (m)", equal.scales=TRUE)

plot(Ppool[4],  main="Plots 41-60yrs",
     xlim=c(0, 6), xlab="radius (m)", equal.scales=TRUE)

plot(Ppool[5],  main="Plots > 60yrs",
     xlim=c(0, 6), xlab="radius (m)", equal.scales=TRUE)




par(OP)

plot(Lpool, cbind(pooliso,pooltheo,hiiso,loiso) - r  ~ r , legend=FALSE,
     shade=c("hiiso", "loiso"), xlim=c(0, 6), xlab="radius (m)", equal.scales=TRUE)

plot(Lpool[1], cbind(pooliso,pooltheo,hiiso,loiso) - r  ~ r , main="Plots < 11 yrs",
     shade=c("hiiso", "loiso"), xlim=c(0, 6), xlab="radius (m)", equal.scales=TRUE)

plot(Lpool[2], cbind(pooliso,pooltheo,hiiso,loiso) - r  ~ r , main="Plots 11-20 yrs",
     shade=c("hiiso", "loiso"), xlim=c(0, 6), xlab="radius (m)", equal.scales=TRUE, 
     legendpos="float")
plot(Lpool[3], cbind(pooliso,pooltheo,hiiso,loiso) - r  ~ r , main="Plots 21-40 yrs",
     shade=c("hiiso", "loiso"), xlim=c(0, 6), xlab="radius (m)", equal.scales=TRUE, 
     legendpos="float")
plot(Lpool[4], cbind(pooliso,pooltheo,hiiso,loiso) - r  ~ r , main="Plots 41-60 yrs",
     shade=c("hiiso", "loiso"), xlim=c(0, 6), xlab="radius (m)", equal.scales=TRUE, 
     legendpos="float")
plot(Lpool[5], cbind(pooliso,pooltheo,hiiso,loiso) - r  ~ r , main="Plots > 60 yrs",
     shade=c("hiiso", "loiso"), xlim=c(0, 6), xlab="radius (m)", equal.scales=TRUE, 
     legendpos="float")




#----------------------------------------------------------------------------------------------
# CREATE HYPERFRAME SUBSETS

BCoa4F = subset(BCoa4, (Origin) == "Fo") #A subset hyperframe of only plots of Forest Origin
BCoa4F
BCoa4S = subset(BCoa4, (Origin) == "Sa") #A subset hyperframe of only plots of Savannah Origin
BCoa4S

# Splitting both hyperframe subsets by Plot Age group
LFoSplit = split(BCoa4F$L4, BCoa4F$Group)
LSaSplit = split(BCoa4S$L4, BCoa4S$Group)

# Plot Pooled L-function for plots to Forest Origin
LFpool = anylapply(LFoSplit, pool) # pool split L-Function by Age group
summary(LFpool)

plot(LFpool[1], cbind(pooliso,pooltheo,hiiso,loiso) - r  ~ r , main="Fo: < 11 yrs",
     shade=c("hiiso", "loiso"), xlim=c(0, 6), xlab="radius (m)", equal.scales=TRUE)

plot(LFpool[2], cbind(pooliso,pooltheo,hiiso,loiso) - r  ~ r , main="Fo: 11-20 yrs",
     shade=c("hiiso", "loiso"), xlim=c(0, 6), xlab="radius (m)", equal.scales=TRUE, 
     legendpos="float")
plot(LFpool[3], cbind(pooliso,pooltheo,hiiso,loiso) - r  ~ r , main="Fo: 21-40 yrs",
     shade=c("hiiso", "loiso"), xlim=c(0, 6), xlab="radius (m)", equal.scales=TRUE, 
     legendpos="float")
plot(LFpool[4], cbind(pooliso,pooltheo,hiiso,loiso) - r  ~ r , main="Fo: 41-60 yrs",
     shade=c("hiiso", "loiso"), xlim=c(0, 6), xlab="radius (m)", equal.scales=TRUE, 
     legendpos="float")
plot(LFpool[5], cbind(pooliso,pooltheo,hiiso,loiso) - r  ~ r , main="Fo: > 60yrs",
     shade=c("hiiso", "loiso"), xlim=c(0, 6), xlab="radius (m)", equal.scales=TRUE, 
     legendpos="float")


# Plot Pooled L-function for plots to Savannah Origin
LSpool = anylapply(LSaSplit, pool)# pool split L-Function by Age group
summary(LSpool)

plot(LSpool[1], cbind(pooliso,pooltheo,hiiso,loiso) - r  ~ r , main="Sa: < 11 yrs",
     shade=c("hiiso", "loiso"), xlim=c(0, 6), xlab="radius (m)", equal.scales=TRUE)

plot(LSpool[2], cbind(pooliso,pooltheo,hiiso,loiso) - r  ~ r , main="Sa: 11-20 yrs",
     shade=c("hiiso", "loiso"), xlim=c(0, 6), xlab="radius (m)", equal.scales=TRUE, 
     legendpos="float")
plot(LSpool[3], cbind(pooliso,pooltheo,hiiso,loiso) - r  ~ r , main="Sa: 21-40 yrs",
     shade=c("hiiso", "loiso"), xlim=c(0, 6), xlab="radius (m)", equal.scales=TRUE, 
     legendpos="float")
plot(LSpool[4], cbind(pooliso,pooltheo,hiiso,loiso) - r  ~ r , main="Sa: 41-60 yrs",
     shade=c("hiiso", "loiso"), xlim=c(0, 6), xlab="radius (m)", equal.scales=TRUE, 
     legendpos="float")
plot(LSpool[5], cbind(pooliso,pooltheo,hiiso,loiso) - r  ~ r , main="Sa: > 60 yrs",
     shade=c("hiiso", "loiso"), xlim=c(0, 6), xlab="radius (m)", equal.scales=TRUE, 
     legendpos="float")



# SUBSETTING THE BCoa4 HYPERFRAME
load("D:/....../BCoa4.Rda")


#Subset and eclude rows where Plots '26Aa', '34B', '34C' and 'NS1'
# "!" is (the negation operator) used to exclude the Plots in the list

#BCoa4b = subset(BCoa4, !(Plots) %in% c('26Aa', '34B', '34C' ,'NS1'))
BCoa4b = subset(BCoa4, !(Plots) ==  '34B') #Exlcude plot 34B, lacking Coa3 Trees
BCoa4b

save(BCoa4b, file = "D:/PhD Manuscript/...../BCoa4b.Rda")
load("D:/PhD Manuscript/...../BCoa4b.Rda")

# Save CSV file version

write.csv(BCoa4b, file = "D:/....../BCoa4b.csv")

BCoa4b

plot(BCoa4b$DistMc4)
plot(BCoa4b$dCoa4, main="Distance Function Map, Cocoa Trees" )

plot(BCoa4b$DistMAsT) # Distance Map of Associated Trees
plot(BCoa4b$dAsT, main="Distance Function Map, Associated Trees") # Distance function of Associated trees







sapply(T4Al, npoints) #to get total number of trees for each spatial point

sapply(T4Co, npoints) #to get number of cocoa trees for each spatial point



# Plotting the K-funtions of each of the patterns in the Coa4 dataset( or hyperframe)
#plot(Coa4, quote(plot(Kest(Coa))))
#plot(Coa4, quote(plot(Lest(Coa))))
# Add a new column containing pixel images of the distance mapps for each point pattern
BCoa4$DistM = with(BCoa4, distmap(Coa4))

# Create a list(m), representing a vector of nearest-neighbour distances
BCoa4$MM = with(BCoa4, nndist(Coa4))
BCoa4$NT = with(BCoa4, npoints(Coa4)) # create a column for the No. of point in each patten in Coa Column of Coa4 dataset

# Find the minimum interpoint distance in each point pattern (plot)
with(BCoa4, min(nndist(Coa4)))


# Mixed model fitting for all plots

mppm(Coa ~ 1, data = Coa4, random = ~1|id)

Coa4$K = with(Coa4, Kest(Coa))

Coa4$L = with(Coa4, Lest(Coa))

mppm(Coa ~ Group, data = Coa4, random = ~1|id)

#fitCoa4 = mppm(Coa ~ Group, data = Coa4, random = ~1|id)
Coa4
fitCoa4 = mppm(Coa ~ 1, data = Coa4, random = ~ 1 | id/Group)
anova(fitCoa4)
ranef(fitCoa4)
coef(fitCoa4)

# PERMUTATION TEST FOR DIFFERENCE BETWEEN AGE GROUPS


Ksplit = split(Coa4$K, Coa4$Group) # Split files based on K-Funtion 
Kpool = anylapply(Ksplit, pool)
plot(Kpool, cbind(pooliso,pooltheo,hiiso,loiso) ~ r,
     shade=c("hiiso", "loiso"), xlim=c(0, 6), xlab="r(m)", equal.scales=TRUE)

Lsplit = split(Coa4$L, Coa4$Group) # Split files based on transformed K-Funtion (L Column)
Lpool = anylapply(Lsplit, pool)
plot(Lpool, cbind(pooliso,pooltheo,hiiso,loiso) - r ~ r,
     shade=c("hiiso", "loiso"), xlim=c(0, 6), xlab="r(m)", equal.scales=TRUE)


# Plot estimated K-function for coaf stand point patterns

L_Coa4 = as.anylist(lapply(Coa4, Lest))
plot(L_Coa4[1:6], main="", main.panel=letters[1:6], xlab="r(m)", legend=FALSE)


#radii = with(Coa4, mean(nndist(Coa))) # Error cannot coerce type 'closure' to vector of type 'double'

coaffit= mppm(Coa ~ 1, data=Coa4, random = ~1|id)
anova(coaffit)

Kcoaf = lapply(Coa4, Kest, ratio=TRUE)
Kcoaf[[1]]

# To pool the computed K-Functions

K = pool(as.anylist(Kcoaf))

plot(K, cbind(pooliso, pooltheo, loiso, hiiso) ~ r,
     shade=c("loiso", "hiiso"))




