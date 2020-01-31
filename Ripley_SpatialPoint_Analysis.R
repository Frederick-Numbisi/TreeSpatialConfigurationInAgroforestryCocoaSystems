
library(sp)
library(spatstat)
library(maptools)
library(ggplot2)

setwd("D:/PhD Manuscript/...../Bakoa_TreesClass4")


Plot23_Data = read.csv("Plot23_tclass4.csv")
Plot23_Data$TTP = as.factor(Plot23_Data$TTP)

Plot23_T4  <- readShapePoints("Plot23_tclass4.shp")
SP <- as(Plot23_T4, "SpatialPoints")
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


#--------------------------------------------------

#PLOT 26


Plot26_Data = read.csv("Plot26_tclass4.csv")
Plot26_Data$TTP = as.factor(Plot26_Data$TTP)

Plot26_T4  <- readShapePoints("Plot26_tclass4.shp")
SP <- as(Plot26_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot26_Data[5]

# Get Plot261 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot26.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot26 windown (polygon) to trees in Plot26
P$window <- W
plot(P, pch=16, main="Plot26 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot26 All Trees")
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
L_AsT_Plot26 <- Ldot.inhom(P, "AsT", lambdaX=fit, update=FALSE)
plot(Lij_Plot26, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot26 Associated Trees", legendpos="bottom")
L_AsT_Plot26$iso  <- L_AsT_Plot26$iso  - L_AsT_Plot26$r
L_AsT_Plot26$theo <- L_AsT_Plot26$theo - L_AsT_Plot26$r
write.csv(L_AsT_Plot26, file = "Li_AsT_Plot26.csv", row.names=TRUE)

#Li (nonhomogenous) for Cocoa Trees2 in stand 

L_Co2_Plot26 <- Ldot.inhom(P, "Co", lambdaX=fit, update=FALSE) #using fit computed for all tree

plot(Lij_Plot26, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot26 Cocoa Trees3", legendpos="bottom")
L_AsT_Plot26$iso  <- L_AsT_Plot26$iso  - L_AsT_Plot26$r
L_AsT_Plot26$theo <- L_AsT_Plot26$theo - L_AsT_Plot26$r
write.csv(L_AsT_Plot26, file = "Li_Co3_Plot26.csv", row.names=TRUE)


#Bivariate(cross) inhomogenous Kij(r) Ripley function
#The Kcross(r) test, which examines distribution a tree group relative to another group within stand

# method (3): fit parametric intensity model
fit <- ppm(P, ~marks * polynom(x,y,2))
# evaluate fitted intensities at data points
# (these are the intensities of the sub-processes of each type)
inten <- fitted(fit, dataonly=TRUE)
# split according to types of points
lambda2 <- split(inten, P$marks)
Kij_Plot26 <- Kcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))
plot(Kij_Plot26, xlab="Distance (m)", main="Kij (2nd order moment) Co3 to Ast")

write.csv(Kij_Plot26, file = "Kij_bivariate_inhom_Plot26_Co3-AsT.csv", row.names=TRUE)



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
Lij_Plot26 <- Lcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))

plot(Lij_Plot26, . - r ~ r, xlab="Distance (m)", ylab="Lij(r)", main="Lij inhom. Function, Plot26 Cocoa3-AsT", legendpos="bottom")


# We will need to modify the plot graphics parameters to accommodate the
# custom y-label; hence the call to `par`.
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Plot26, . -r ~ r, ylab=expression(hat("Lij")), xlab = "d (m)", main="Lij(r), Plot26 Cocoa Trees", legendmath=TRUE, legendpos="bottom")
#par(OP)

#We can make permanent changes to the Lij object's values as follows:
Lij_Plot26$iso  <- Lij_Plot26$iso  - Lij_Plot26$r
Lij_Plot26$theo <- Lij_Plot26$theo - Lij_Plot26$r

write.csv(Lij_Plot26, file = "Lij_Bivariate_inhom_Plot26_Co3_to_AsT.csv", row.names=TRUE)

#MONTE CARLO (MC) SIMULATION

#To test for clustering/inhibition in the presence of spatial inhomogeneity, 
#the null hypothesis should be an inhomogeneous Poisson process. 
#So we estimate the inhomogeneous intensity functions of the two types of points,
# and then generate simulated point patterns according to these intensities.

#First fit a Poisson point process model to the observed data
#fit2 <- ppm(P, ~marks * polynom(x,y,3))

#n <- 199 #number of iterations for Monte carlo simulation
#p  <- 0.01 # Desired p significance level to display

#Lij_Env_Plot26 <- envelope(fit2, Lcross.inhom, lambdaX=fit2, nsim=n, rank=(p * (n + 1)))
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Env_Plot26, . - r ~ r, ylab=expression(hat("L")),xlab="Distance (m)", main=paste("Lij(r) Plot26 Co3-AsT,  p= ", p),legendmath=TRUE)
#par(OP)
#write.csv(Lij_Env_Plot26, file = "Lij_199MCEnv_Plot26_Co3_to_AsT.csv", row.names=TRUE)



# COMPUTING THE g (Pair Correlation Fucntion) OF COCOA AND ASSOCIATED TREES
#The function value (fv) table method:we will use computed inhomogenous cross K function

# fv object
#Kij_Plot26 is the fv object computed above (inhomogenous cross K from Co to AsT)
#method "b" apply smoothing to Y(r) = K(r)/(2 * pi * r) constraining Y(0) = 0, 
#estimate the derivative of Y, and solve;
#method "b" used due to likely inhibition at small distances of r spar=0.8,

gij_Plot26 <- pcf(Kij_Plot26, spar=0.8,  method="b") 
plot(gij_Plot26, main="Plot26 Co3 - AsT g(r) Pair correlation function", xlab="Distance(m)", legendpos="bottom")
write.csv(gij_Plot26, file = "pcf_Plot26_Co3_to_AsT.csv", row.names=TRUE)


#--------------------------------------------------

#PLOT 28


Plot28_Data = read.csv("Plot28_tclass4.csv")
Plot28_Data$TTP = as.factor(Plot28_Data$TTP)

Plot28_T4  <- readShapePoints("Plot28_tclass4.shp")
SP <- as(Plot28_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot28_Data[5]

# Get Plot281 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot28.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot28 windown (polygon) to trees in Plot28
P$window <- W
plot(P, pch=16, main="Plot28 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot28 All Trees")
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
L_AsT_Plot28 <- Ldot.inhom(P, "AsT", lambdaX=fit, update=FALSE)
plot(Lij_Plot28, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot28 Associated Trees", legendpos="bottom")
L_AsT_Plot28$iso  <- L_AsT_Plot28$iso  - L_AsT_Plot28$r
L_AsT_Plot28$theo <- L_AsT_Plot28$theo - L_AsT_Plot28$r
write.csv(L_AsT_Plot28, file = "Li_AsT_Plot28.csv", row.names=TRUE)

#Li (nonhomogenous) for Cocoa Trees2 in stand 

L_Co2_Plot28 <- Ldot.inhom(P, "Co", lambdaX=fit, update=FALSE) #using fit computed for all tree

plot(Lij_Plot28, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot28 Cocoa Trees3", legendpos="bottom")
L_AsT_Plot28$iso  <- L_AsT_Plot28$iso  - L_AsT_Plot28$r
L_AsT_Plot28$theo <- L_AsT_Plot28$theo - L_AsT_Plot28$r
write.csv(L_AsT_Plot28, file = "Li_Co3_Plot28.csv", row.names=TRUE)


#Bivariate(cross) inhomogenous Kij(r) Ripley function
#The Kcross(r) test, which examines distribution a tree group relative to another group within stand

# method (3): fit parametric intensity model
fit <- ppm(P, ~marks * polynom(x,y,2))
# evaluate fitted intensities at data points
# (these are the intensities of the sub-processes of each type)
inten <- fitted(fit, dataonly=TRUE)
# split according to types of points
lambda2 <- split(inten, P$marks)
Kij_Plot28 <- Kcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))
plot(Kij_Plot28, xlab="Distance (m)", main="Kij (2nd order moment) Co3 to Ast")

write.csv(Kij_Plot28, file = "Kij_bivariate_inhom_Plot28_Co3-AsT.csv", row.names=TRUE)



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
Lij_Plot28 <- Lcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))

plot(Lij_Plot28, . - r ~ r, xlab="Distance (m)", ylab="Lij(r)", main="Lij inhom. Function, Plot28 Cocoa3-AsT", legendpos="bottom")


# We will need to modify the plot graphics parameters to accommodate the
# custom y-label; hence the call to `par`.
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Plot28, . -r ~ r, ylab=expression(hat("Lij")), xlab = "d (m)", main="Lij(r), Plot28 Cocoa Trees", legendmath=TRUE, legendpos="bottom")
#par(OP)

#We can make permanent changes to the Lij object's values as follows:
Lij_Plot28$iso  <- Lij_Plot28$iso  - Lij_Plot28$r
Lij_Plot28$theo <- Lij_Plot28$theo - Lij_Plot28$r

write.csv(Lij_Plot28, file = "Lij_Bivariate_inhom_Plot28_Co3_to_AsT.csv", row.names=TRUE)

#MONTE CARLO (MC) SIMULATION

#To test for clustering/inhibition in the presence of spatial inhomogeneity, 
#the null hypothesis should be an inhomogeneous Poisson process. 
#So we estimate the inhomogeneous intensity functions of the two types of points,
# and then generate simulated point patterns according to these intensities.

#First fit a Poisson point process model to the observed data
#fit2 <- ppm(P, ~marks * polynom(x,y,3))

#n <- 199 #number of iterations for Monte carlo simulation
#p  <- 0.01 # Desired p significance level to display

#Lij_Env_Plot28 <- envelope(fit2, Lcross.inhom, lambdaX=fit2, nsim=n, rank=(p * (n + 1)))
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Env_Plot28, . - r ~ r, ylab=expression(hat("L")),xlab="Distance (m)", main=paste("Lij(r) Plot28 Co3-AsT,  p= ", p),legendmath=TRUE)
#par(OP)
#write.csv(Lij_Env_Plot28, file = "Lij_199MCEnv_Plot28_Co3_to_AsT.csv", row.names=TRUE)



# COMPUTING THE g (Pair Correlation Fucntion) OF COCOA AND ASSOCIATED TREES
#The function value (fv) table method:we will use computed inhomogenous cross K function

# fv object
#Kij_Plot28 is the fv object computed above (inhomogenous cross K from Co to AsT)
#method "b" apply smoothing to Y(r) = K(r)/(2 * pi * r) constraining Y(0) = 0, 
#estimate the derivative of Y, and solve;
#method "b" used due to likely inhibition at small distances of r spar=0.8,

gij_Plot28 <- pcf(Kij_Plot28, spar=0.8,  method="b") 
plot(gij_Plot28, main="Plot28 Co3 - AsT g(r) Pair correlation function", xlab="Distance(m)", legendpos="bottom")
write.csv(gij_Plot28, file = "pcf_Plot28_Co3_to_AsT.csv", row.names=TRUE)


#--------------------------------------------------

#PLOT 38


Plot38_Data = read.csv("Plot38_tclass4.csv")
Plot38_Data$TTP = as.factor(Plot38_Data$TTP)

Plot38_T4  <- readShapePoints("Plot38_tclass4.shp")
SP <- as(Plot38_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot38_Data[5]

# Get Plot381 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot38.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot38 windown (polygon) to trees in Plot38
P$window <- W
plot(P, pch=16, main="Plot38 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot38 All Trees")
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
L_AsT_Plot38 <- Ldot.inhom(P, "AsT", lambdaX=fit, update=FALSE)
plot(Lij_Plot38, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot38 Associated Trees", legendpos="bottom")
L_AsT_Plot38$iso  <- L_AsT_Plot38$iso  - L_AsT_Plot38$r
L_AsT_Plot38$theo <- L_AsT_Plot38$theo - L_AsT_Plot38$r
write.csv(L_AsT_Plot38, file = "Li_AsT_Plot38.csv", row.names=TRUE)

#Li (nonhomogenous) for Cocoa Trees2 in stand 

L_Co2_Plot38 <- Ldot.inhom(P, "Co", lambdaX=fit, update=FALSE) #using fit computed for all tree

plot(Lij_Plot38, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot38 Cocoa Trees3", legendpos="bottom")
L_AsT_Plot38$iso  <- L_AsT_Plot38$iso  - L_AsT_Plot38$r
L_AsT_Plot38$theo <- L_AsT_Plot38$theo - L_AsT_Plot38$r
write.csv(L_AsT_Plot38, file = "Li_Co3_Plot38.csv", row.names=TRUE)


#Bivariate(cross) inhomogenous Kij(r) Ripley function
#The Kcross(r) test, which examines distribution a tree group relative to another group within stand

# method (3): fit parametric intensity model
fit <- ppm(P, ~marks * polynom(x,y,2))
# evaluate fitted intensities at data points
# (these are the intensities of the sub-processes of each type)
inten <- fitted(fit, dataonly=TRUE)
# split according to types of points
lambda2 <- split(inten, P$marks)
Kij_Plot38 <- Kcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))
plot(Kij_Plot38, xlab="Distance (m)", main="Kij (2nd order moment) Co3 to Ast")

write.csv(Kij_Plot38, file = "Kij_bivariate_inhom_Plot38_Co3-AsT.csv", row.names=TRUE)



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
Lij_Plot38 <- Lcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))

plot(Lij_Plot38, . - r ~ r, xlab="Distance (m)", ylab="Lij(r)", main="Lij inhom. Function, Plot38 Cocoa3-AsT", legendpos="bottom")


# We will need to modify the plot graphics parameters to accommodate the
# custom y-label; hence the call to `par`.
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Plot38, . -r ~ r, ylab=expression(hat("Lij")), xlab = "d (m)", main="Lij(r), Plot38 Cocoa Trees", legendmath=TRUE, legendpos="bottom")
#par(OP)

#We can make permanent changes to the Lij object's values as follows:
Lij_Plot38$iso  <- Lij_Plot38$iso  - Lij_Plot38$r
Lij_Plot38$theo <- Lij_Plot38$theo - Lij_Plot38$r

write.csv(Lij_Plot38, file = "Lij_Bivariate_inhom_Plot38_Co3_to_AsT.csv", row.names=TRUE)

#MONTE CARLO (MC) SIMULATION

#To test for clustering/inhibition in the presence of spatial inhomogeneity, 
#the null hypothesis should be an inhomogeneous Poisson process. 
#So we estimate the inhomogeneous intensity functions of the two types of points,
# and then generate simulated point patterns according to these intensities.

#First fit a Poisson point process model to the observed data
#fit2 <- ppm(P, ~marks * polynom(x,y,3))

#n <- 199 #number of iterations for Monte carlo simulation
#p  <- 0.01 # Desired p significance level to display

#Lij_Env_Plot38 <- envelope(fit2, Lcross.inhom, lambdaX=fit2, nsim=n, rank=(p * (n + 1)))
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Env_Plot38, . - r ~ r, ylab=expression(hat("L")),xlab="Distance (m)", main=paste("Lij(r) Plot38 Co3-AsT,  p= ", p),legendmath=TRUE)
#par(OP)
#write.csv(Lij_Env_Plot38, file = "Lij_199MCEnv_Plot38_Co3_to_AsT.csv", row.names=TRUE)



# COMPUTING THE g (Pair Correlation Fucntion) OF COCOA AND ASSOCIATED TREES
#The function value (fv) table method:we will use computed inhomogenous cross K function

# fv object
#Kij_Plot38 is the fv object computed above (inhomogenous cross K from Co to AsT)
#method "b" apply smoothing to Y(r) = K(r)/(2 * pi * r) constraining Y(0) = 0, 
#estimate the derivative of Y, and solve;
#method "b" used due to likely inhibition at small distances of r spar=0.8,

gij_Plot38 <- pcf(Kij_Plot38, spar=0.8,  method="b") 
plot(gij_Plot38, main="Plot38 Co3 - AsT g(r) Pair correlation function", xlab="Distance(m)", legendpos="bottom")
write.csv(gij_Plot38, file = "pcf_Plot38_Co3_to_AsT.csv", row.names=TRUE)


#--------------------------------------------------

#PLOT 4


Plot4_Data = read.csv("Plot4_tclass4.csv")
Plot4_Data$TTP = as.factor(Plot4_Data$TTP)

Plot4_T4  <- readShapePoints("Plot4_tclass4.shp")
SP <- as(Plot4_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot4_Data[5]

# Get Plot41 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot4.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot4 windown (polygon) to trees in Plot4
P$window <- W
plot(P, pch=16, main="Plot4 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot4 All Trees")
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
L_AsT_Plot4 <- Ldot.inhom(P, "AsT", lambdaX=fit, update=FALSE)
plot(Lij_Plot4, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot4 Associated Trees", legendpos="bottom")
L_AsT_Plot4$iso  <- L_AsT_Plot4$iso  - L_AsT_Plot4$r
L_AsT_Plot4$theo <- L_AsT_Plot4$theo - L_AsT_Plot4$r
write.csv(L_AsT_Plot4, file = "Li_AsT_Plot4.csv", row.names=TRUE)

#Li (nonhomogenous) for Cocoa Trees2 in stand 

L_Co2_Plot4 <- Ldot.inhom(P, "Co", lambdaX=fit, update=FALSE) #using fit computed for all tree

plot(Lij_Plot4, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot4 Cocoa Trees3", legendpos="bottom")
L_AsT_Plot4$iso  <- L_AsT_Plot4$iso  - L_AsT_Plot4$r
L_AsT_Plot4$theo <- L_AsT_Plot4$theo - L_AsT_Plot4$r
write.csv(L_AsT_Plot4, file = "Li_Co3_Plot4.csv", row.names=TRUE)


#Bivariate(cross) inhomogenous Kij(r) Ripley function
#The Kcross(r) test, which examines distribution a tree group relative to another group within stand

# method (3): fit parametric intensity model
fit <- ppm(P, ~marks * polynom(x,y,2))
# evaluate fitted intensities at data points
# (these are the intensities of the sub-processes of each type)
inten <- fitted(fit, dataonly=TRUE)
# split according to types of points
lambda2 <- split(inten, P$marks)
Kij_Plot4 <- Kcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))
plot(Kij_Plot4, xlab="Distance (m)", main="Kij (2nd order moment) Co3 to Ast")

write.csv(Kij_Plot4, file = "Kij_bivariate_inhom_Plot4_Co3-AsT.csv", row.names=TRUE)



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
Lij_Plot4 <- Lcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))

plot(Lij_Plot4, . - r ~ r, xlab="Distance (m)", ylab="Lij(r)", main="Lij inhom. Function, Plot4 Cocoa3-AsT", legendpos="bottom")


# We will need to modify the plot graphics parameters to accommodate the
# custom y-label; hence the call to `par`.
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Plot4, . -r ~ r, ylab=expression(hat("Lij")), xlab = "d (m)", main="Lij(r), Plot4 Cocoa Trees", legendmath=TRUE, legendpos="bottom")
#par(OP)

#We can make permanent changes to the Lij object's values as follows:
Lij_Plot4$iso  <- Lij_Plot4$iso  - Lij_Plot4$r
Lij_Plot4$theo <- Lij_Plot4$theo - Lij_Plot4$r

write.csv(Lij_Plot4, file = "Lij_Bivariate_inhom_Plot4_Co3_to_AsT.csv", row.names=TRUE)

#MONTE CARLO (MC) SIMULATION

#To test for clustering/inhibition in the presence of spatial inhomogeneity, 
#the null hypothesis should be an inhomogeneous Poisson process. 
#So we estimate the inhomogeneous intensity functions of the two types of points,
# and then generate simulated point patterns according to these intensities.

#First fit a Poisson point process model to the observed data
#fit2 <- ppm(P, ~marks * polynom(x,y,3))

#n <- 199 #number of iterations for Monte carlo simulation
#p  <- 0.01 # Desired p significance level to display

#Lij_Env_Plot4 <- envelope(fit2, Lcross.inhom, lambdaX=fit2, nsim=n, rank=(p * (n + 1)))
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Env_Plot4, . - r ~ r, ylab=expression(hat("L")),xlab="Distance (m)", main=paste("Lij(r) Plot4 Co3-AsT,  p= ", p),legendmath=TRUE)
#par(OP)
#write.csv(Lij_Env_Plot4, file = "Lij_199MCEnv_Plot4_Co3_to_AsT.csv", row.names=TRUE)



# COMPUTING THE g (Pair Correlation Fucntion) OF COCOA AND ASSOCIATED TREES
#The function value (fv) table method:we will use computed inhomogenous cross K function

# fv object
#Kij_Plot4 is the fv object computed above (inhomogenous cross K from Co to AsT)
#method "b" apply smoothing to Y(r) = K(r)/(2 * pi * r) constraining Y(0) = 0, 
#estimate the derivative of Y, and solve;
#method "b" used due to likely inhibition at small distances of r spar=0.8,

gij_Plot4 <- pcf(Kij_Plot4, spar=0.8,  method="b") 
plot(gij_Plot4, main="Plot4 Co3 - AsT g(r) Pair correlation function", xlab="Distance(m)", legendpos="bottom")
write.csv(gij_Plot4, file = "pcf_Plot4_Co3_to_AsT.csv", row.names=TRUE)


#------------------------------------------------

#PLOT 22


Plot22_Data = read.csv("Plot22_tclass4.csv")
Plot22_Data$TTP = as.factor(Plot22_Data$TTP)

Plot22_T4  <- readShapePoints("Plot22_tclass4.shp")
SP <- as(Plot22_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot22_Data[5]

# Get Plot221 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot22.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot22 windown (polygon) to trees in Plot22
P$window <- W
plot(P, pch=16, main="Plot22 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot22 All Trees")
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
L_AsT_Plot22 <- Ldot.inhom(P, "AsT", lambdaX=fit, update=FALSE)
plot(Lij_Plot22, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot22 Associated Trees", legendpos="bottom")
L_AsT_Plot22$iso  <- L_AsT_Plot22$iso  - L_AsT_Plot22$r
L_AsT_Plot22$theo <- L_AsT_Plot22$theo - L_AsT_Plot22$r
write.csv(L_AsT_Plot22, file = "Li_AsT_Plot22.csv", row.names=TRUE)

#Li (nonhomogenous) for Cocoa Trees2 in stand 

L_Co2_Plot22 <- Ldot.inhom(P, "Co", lambdaX=fit, update=FALSE) #using fit computed for all tree

plot(Lij_Plot22, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot22 Cocoa Trees3", legendpos="bottom")
L_AsT_Plot22$iso  <- L_AsT_Plot22$iso  - L_AsT_Plot22$r
L_AsT_Plot22$theo <- L_AsT_Plot22$theo - L_AsT_Plot22$r
write.csv(L_AsT_Plot22, file = "Li_Co3_Plot22.csv", row.names=TRUE)


#Bivariate(cross) inhomogenous Kij(r) Ripley function
#The Kcross(r) test, which examines distribution a tree group relative to another group within stand

# method (3): fit parametric intensity model
fit <- ppm(P, ~marks * polynom(x,y,2))
# evaluate fitted intensities at data points
# (these are the intensities of the sub-processes of each type)
inten <- fitted(fit, dataonly=TRUE)
# split according to types of points
lambda2 <- split(inten, P$marks)
Kij_Plot22 <- Kcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))
plot(Kij_Plot22, xlab="Distance (m)", main="Kij (2nd order moment) Co3 to Ast")

write.csv(Kij_Plot22, file = "Kij_bivariate_inhom_Plot22_Co3-AsT.csv", row.names=TRUE)



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
Lij_Plot22 <- Lcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))

plot(Lij_Plot22, . - r ~ r, xlab="Distance (m)", ylab="Lij(r)", main="Lij inhom. Function, Plot22 Cocoa3-AsT", legendpos="bottom")


# We will need to modify the plot graphics parameters to accommodate the
# custom y-label; hence the call to `par`.
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Plot22, . -r ~ r, ylab=expression(hat("Lij")), xlab = "d (m)", main="Lij(r), Plot22 Cocoa Trees", legendmath=TRUE, legendpos="bottom")
#par(OP)

#We can make permanent changes to the Lij object's values as follows:
Lij_Plot22$iso  <- Lij_Plot22$iso  - Lij_Plot22$r
Lij_Plot22$theo <- Lij_Plot22$theo - Lij_Plot22$r

write.csv(Lij_Plot22, file = "Lij_Bivariate_inhom_Plot22_Co3_to_AsT.csv", row.names=TRUE)

#MONTE CARLO (MC) SIMULATION

#To test for clustering/inhibition in the presence of spatial inhomogeneity, 
#the null hypothesis should be an inhomogeneous Poisson process. 
#So we estimate the inhomogeneous intensity functions of the two types of points,
# and then generate simulated point patterns according to these intensities.

#First fit a Poisson point process model to the observed data
#fit2 <- ppm(P, ~marks * polynom(x,y,3))

#n <- 199 #number of iterations for Monte carlo simulation
#p  <- 0.01 # Desired p significance level to display

#Lij_Env_Plot22 <- envelope(fit2, Lcross.inhom, lambdaX=fit2, nsim=n, rank=(p * (n + 1)))
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Env_Plot22, . - r ~ r, ylab=expression(hat("L")),xlab="Distance (m)", main=paste("Lij(r) Plot22 Co3-AsT,  p= ", p),legendmath=TRUE)
#par(OP)
#write.csv(Lij_Env_Plot22, file = "Lij_199MCEnv_Plot22_Co3_to_AsT.csv", row.names=TRUE)



# COMPUTING THE g (Pair Correlation Fucntion) OF COCOA AND ASSOCIATED TREES
#The function value (fv) table method:we will use computed inhomogenous cross K function

# fv object
#Kij_Plot22 is the fv object computed above (inhomogenous cross K from Co to AsT)
#method "b" apply smoothing to Y(r) = K(r)/(2 * pi * r) constraining Y(0) = 0, 
#estimate the derivative of Y, and solve;
#method "b" used due to likely inhibition at small distances of r spar=0.8,

gij_Plot22 <- pcf(Kij_Plot22, spar=0.8,  method="b") 
plot(gij_Plot22, main="Plot22 Co3 - AsT g(r) Pair correlation function", xlab="Distance(m)", legendpos="bottom")
write.csv(gij_Plot22, file = "pcf_Plot22_Co3_to_AsT.csv", row.names=TRUE)


#--------------------------------------------------

#PLOT 6


Plot6_Data = read.csv("Plot6_tclass4.csv")
Plot6_Data$TTP = as.factor(Plot6_Data$TTP)

Plot6_T4  <- readShapePoints("Plot6_tclass4.shp")
SP <- as(Plot6_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot6_Data[5]

# Get Plot61 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot6.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot6 windown (polygon) to trees in Plot6
P$window <- W
plot(P, pch=16, main="Plot6 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot6 All Trees")
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
L_AsT_Plot6 <- Ldot.inhom(P, "AsT", lambdaX=fit, update=FALSE)
plot(Lij_Plot6, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot6 Associated Trees", legendpos="bottom")
L_AsT_Plot6$iso  <- L_AsT_Plot6$iso  - L_AsT_Plot6$r
L_AsT_Plot6$theo <- L_AsT_Plot6$theo - L_AsT_Plot6$r
write.csv(L_AsT_Plot6, file = "Li_AsT_Plot6.csv", row.names=TRUE)

#Li (nonhomogenous) for Cocoa Trees2 in stand 

L_Co2_Plot6 <- Ldot.inhom(P, "Co", lambdaX=fit, update=FALSE) #using fit computed for all tree

plot(Lij_Plot6, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot6 Cocoa Trees3", legendpos="bottom")
L_AsT_Plot6$iso  <- L_AsT_Plot6$iso  - L_AsT_Plot6$r
L_AsT_Plot6$theo <- L_AsT_Plot6$theo - L_AsT_Plot6$r
write.csv(L_AsT_Plot6, file = "Li_Co3_Plot6.csv", row.names=TRUE)


#Bivariate(cross) inhomogenous Kij(r) Ripley function
#The Kcross(r) test, which examines distribution a tree group relative to another group within stand

# method (3): fit parametric intensity model
fit <- ppm(P, ~marks * polynom(x,y,2))
# evaluate fitted intensities at data points
# (these are the intensities of the sub-processes of each type)
inten <- fitted(fit, dataonly=TRUE)
# split according to types of points
lambda2 <- split(inten, P$marks)
Kij_Plot6 <- Kcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))
plot(Kij_Plot6, xlab="Distance (m)", main="Kij (2nd order moment) Co3 to Ast")

write.csv(Kij_Plot6, file = "Kij_bivariate_inhom_Plot6_Co3-AsT.csv", row.names=TRUE)



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
Lij_Plot6 <- Lcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))

plot(Lij_Plot6, . - r ~ r, xlab="Distance (m)", ylab="Lij(r)", main="Lij inhom. Function, Plot6 Cocoa3-AsT", legendpos="bottom")


# We will need to modify the plot graphics parameters to accommodate the
# custom y-label; hence the call to `par`.
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Plot6, . -r ~ r, ylab=expression(hat("Lij")), xlab = "d (m)", main="Lij(r), Plot6 Cocoa Trees", legendmath=TRUE, legendpos="bottom")
#par(OP)

#We can make permanent changes to the Lij object's values as follows:
Lij_Plot6$iso  <- Lij_Plot6$iso  - Lij_Plot6$r
Lij_Plot6$theo <- Lij_Plot6$theo - Lij_Plot6$r

write.csv(Lij_Plot6, file = "Lij_Bivariate_inhom_Plot6_Co3_to_AsT.csv", row.names=TRUE)

#MONTE CARLO (MC) SIMULATION

#To test for clustering/inhibition in the presence of spatial inhomogeneity, 
#the null hypothesis should be an inhomogeneous Poisson process. 
#So we estimate the inhomogeneous intensity functions of the two types of points,
# and then generate simulated point patterns according to these intensities.

#First fit a Poisson point process model to the observed data
#fit2 <- ppm(P, ~marks * polynom(x,y,3))

#n <- 199 #number of iterations for Monte carlo simulation
#p  <- 0.01 # Desired p significance level to display

#Lij_Env_Plot6 <- envelope(fit2, Lcross.inhom, lambdaX=fit2, nsim=n, rank=(p * (n + 1)))
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Env_Plot6, . - r ~ r, ylab=expression(hat("L")),xlab="Distance (m)", main=paste("Lij(r) Plot6 Co3-AsT,  p= ", p),legendmath=TRUE)
#par(OP)
#write.csv(Lij_Env_Plot6, file = "Lij_199MCEnv_Plot6_Co3_to_AsT.csv", row.names=TRUE)



# COMPUTING THE g (Pair Correlation Fucntion) OF COCOA AND ASSOCIATED TREES
#The function value (fv) table method:we will use computed inhomogenous cross K function

# fv object
#Kij_Plot6 is the fv object computed above (inhomogenous cross K from Co to AsT)
#method "b" apply smoothing to Y(r) = K(r)/(2 * pi * r) constraining Y(0) = 0, 
#estimate the derivative of Y, and solve;
#method "b" used due to likely inhibition at small distances of r spar=0.8,

gij_Plot6 <- pcf(Kij_Plot6, spar=0.8,  method="b") 
plot(gij_Plot6, main="Plot6 Co3 - AsT g(r) Pair correlation function", xlab="Distance(m)", legendpos="bottom")
write.csv(gij_Plot6, file = "pcf_Plot6_Co3_to_AsT.csv", row.names=TRUE)


#------------------------------------------------

#PLOT 43


Plot43_Data = read.csv("Plot43_tclass4.csv")
Plot43_Data$TTP = as.factor(Plot43_Data$TTP)

Plot43_T4  <- readShapePoints("Plot43_tclass4.shp")
SP <- as(Plot43_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot43_Data[5]

# Get Plot431 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot43.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot43 windown (polygon) to trees in Plot43
P$window <- W
plot(P, pch=16, main="Plot43 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot43 All Trees")
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
L_AsT_Plot43 <- Ldot.inhom(P, "AsT", lambdaX=fit, update=FALSE)
plot(Lij_Plot43, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot43 Associated Trees", legendpos="bottom")
L_AsT_Plot43$iso  <- L_AsT_Plot43$iso  - L_AsT_Plot43$r
L_AsT_Plot43$theo <- L_AsT_Plot43$theo - L_AsT_Plot43$r
write.csv(L_AsT_Plot43, file = "Li_AsT_Plot43.csv", row.names=TRUE)

#Li (nonhomogenous) for Cocoa Trees2 in stand 

L_Co2_Plot43 <- Ldot.inhom(P, "Co", lambdaX=fit, update=FALSE) #using fit computed for all tree

plot(Lij_Plot43, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot43 Cocoa Trees3", legendpos="bottom")
L_AsT_Plot43$iso  <- L_AsT_Plot43$iso  - L_AsT_Plot43$r
L_AsT_Plot43$theo <- L_AsT_Plot43$theo - L_AsT_Plot43$r
write.csv(L_AsT_Plot43, file = "Li_Co3_Plot43.csv", row.names=TRUE)


#Bivariate(cross) inhomogenous Kij(r) Ripley function
#The Kcross(r) test, which examines distribution a tree group relative to another group within stand

# method (3): fit parametric intensity model
fit <- ppm(P, ~marks * polynom(x,y,2))
# evaluate fitted intensities at data points
# (these are the intensities of the sub-processes of each type)
inten <- fitted(fit, dataonly=TRUE)
# split according to types of points
lambda2 <- split(inten, P$marks)
Kij_Plot43 <- Kcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))
plot(Kij_Plot43, xlab="Distance (m)", main="Kij (2nd order moment) Co3 to Ast")

write.csv(Kij_Plot43, file = "Kij_bivariate_inhom_Plot43_Co3-AsT.csv", row.names=TRUE)



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
Lij_Plot43 <- Lcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))

plot(Lij_Plot43, . - r ~ r, xlab="Distance (m)", ylab="Lij(r)", main="Lij inhom. Function, Plot43 Cocoa3-AsT", legendpos="bottom")


# We will need to modify the plot graphics parameters to accommodate the
# custom y-label; hence the call to `par`.
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Plot43, . -r ~ r, ylab=expression(hat("Lij")), xlab = "d (m)", main="Lij(r), Plot43 Cocoa Trees", legendmath=TRUE, legendpos="bottom")
#par(OP)

#We can make permanent changes to the Lij object's values as follows:
Lij_Plot43$iso  <- Lij_Plot43$iso  - Lij_Plot43$r
Lij_Plot43$theo <- Lij_Plot43$theo - Lij_Plot43$r

write.csv(Lij_Plot43, file = "Lij_Bivariate_inhom_Plot43_Co3_to_AsT.csv", row.names=TRUE)

#MONTE CARLO (MC) SIMULATION

#To test for clustering/inhibition in the presence of spatial inhomogeneity, 
#the null hypothesis should be an inhomogeneous Poisson process. 
#So we estimate the inhomogeneous intensity functions of the two types of points,
# and then generate simulated point patterns according to these intensities.

#First fit a Poisson point process model to the observed data
#fit2 <- ppm(P, ~marks * polynom(x,y,3))

#n <- 199 #number of iterations for Monte carlo simulation
#p  <- 0.01 # Desired p significance level to display

#Lij_Env_Plot43 <- envelope(fit2, Lcross.inhom, lambdaX=fit2, nsim=n, rank=(p * (n + 1)))
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Env_Plot43, . - r ~ r, ylab=expression(hat("L")),xlab="Distance (m)", main=paste("Lij(r) Plot43 Co3-AsT,  p= ", p),legendmath=TRUE)
#par(OP)
#write.csv(Lij_Env_Plot43, file = "Lij_199MCEnv_Plot43_Co3_to_AsT.csv", row.names=TRUE)



# COMPUTING THE g (Pair Correlation Fucntion) OF COCOA AND ASSOCIATED TREES
#The function value (fv) table method:we will use computed inhomogenous cross K function

# fv object
#Kij_Plot43 is the fv object computed above (inhomogenous cross K from Co to AsT)
#method "b" apply smoothing to Y(r) = K(r)/(2 * pi * r) constraining Y(0) = 0, 
#estimate the derivative of Y, and solve;
#method "b" used due to likely inhibition at small distances of r spar=0.8,

gij_Plot43 <- pcf(Kij_Plot43, spar=0.8,  method="b") 
plot(gij_Plot43, main="Plot43 Co3 - AsT g(r) Pair correlation function", xlab="Distance(m)", legendpos="bottom")
write.csv(gij_Plot43, file = "pcf_Plot43_Co3_to_AsT.csv", row.names=TRUE)


#--------------------------------------------------

#PLOT 5


Plot5_Data = read.csv("Plot5_tclass4.csv")
Plot5_Data$TTP = as.factor(Plot5_Data$TTP)

Plot5_T4  <- readShapePoints("Plot5_tclass4.shp")
SP <- as(Plot5_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot5_Data[5]

# Get Plot51 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot5.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot5 windown (polygon) to trees in Plot5
P$window <- W
plot(P, pch=16, main="Plot5 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot5 All Trees")
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
L_AsT_Plot5 <- Ldot.inhom(P, "AsT", lambdaX=fit, update=FALSE)
plot(Lij_Plot5, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot5 Associated Trees", legendpos="bottom")
L_AsT_Plot5$iso  <- L_AsT_Plot5$iso  - L_AsT_Plot5$r
L_AsT_Plot5$theo <- L_AsT_Plot5$theo - L_AsT_Plot5$r
write.csv(L_AsT_Plot5, file = "Li_AsT_Plot5.csv", row.names=TRUE)

#Li (nonhomogenous) for Cocoa Trees2 in stand 

L_Co2_Plot5 <- Ldot.inhom(P, "Co", lambdaX=fit, update=FALSE) #using fit computed for all tree

plot(Lij_Plot5, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot5 Cocoa Trees3", legendpos="bottom")
L_AsT_Plot5$iso  <- L_AsT_Plot5$iso  - L_AsT_Plot5$r
L_AsT_Plot5$theo <- L_AsT_Plot5$theo - L_AsT_Plot5$r
write.csv(L_AsT_Plot5, file = "Li_Co3_Plot5.csv", row.names=TRUE)


#Bivariate(cross) inhomogenous Kij(r) Ripley function
#The Kcross(r) test, which examines distribution a tree group relative to another group within stand

# method (3): fit parametric intensity model
fit <- ppm(P, ~marks * polynom(x,y,2))
# evaluate fitted intensities at data points
# (these are the intensities of the sub-processes of each type)
inten <- fitted(fit, dataonly=TRUE)
# split according to types of points
lambda2 <- split(inten, P$marks)
Kij_Plot5 <- Kcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))
plot(Kij_Plot5, xlab="Distance (m)", main="Kij (2nd order moment) Co3 to Ast")

write.csv(Kij_Plot5, file = "Kij_bivariate_inhom_Plot5_Co3-AsT.csv", row.names=TRUE)



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
Lij_Plot5 <- Lcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))

plot(Lij_Plot5, . - r ~ r, xlab="Distance (m)", ylab="Lij(r)", main="Lij inhom. Function, Plot5 Cocoa3-AsT", legendpos="bottom")


# We will need to modify the plot graphics parameters to accommodate the
# custom y-label; hence the call to `par`.
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Plot5, . -r ~ r, ylab=expression(hat("Lij")), xlab = "d (m)", main="Lij(r), Plot5 Cocoa Trees", legendmath=TRUE, legendpos="bottom")
#par(OP)

#We can make permanent changes to the Lij object's values as follows:
Lij_Plot5$iso  <- Lij_Plot5$iso  - Lij_Plot5$r
Lij_Plot5$theo <- Lij_Plot5$theo - Lij_Plot5$r

write.csv(Lij_Plot5, file = "Lij_Bivariate_inhom_Plot5_Co3_to_AsT.csv", row.names=TRUE)

#MONTE CARLO (MC) SIMULATION

#To test for clustering/inhibition in the presence of spatial inhomogeneity, 
#the null hypothesis should be an inhomogeneous Poisson process. 
#So we estimate the inhomogeneous intensity functions of the two types of points,
# and then generate simulated point patterns according to these intensities.

#First fit a Poisson point process model to the observed data
#fit2 <- ppm(P, ~marks * polynom(x,y,3))

#n <- 199 #number of iterations for Monte carlo simulation
#p  <- 0.01 # Desired p significance level to display

#Lij_Env_Plot5 <- envelope(fit2, Lcross.inhom, lambdaX=fit2, nsim=n, rank=(p * (n + 1)))
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Env_Plot5, . - r ~ r, ylab=expression(hat("L")),xlab="Distance (m)", main=paste("Lij(r) Plot5 Co3-AsT,  p= ", p),legendmath=TRUE)
#par(OP)
#write.csv(Lij_Env_Plot5, file = "Lij_199MCEnv_Plot5_Co3_to_AsT.csv", row.names=TRUE)



# COMPUTING THE g (Pair Correlation Fucntion) OF COCOA AND ASSOCIATED TREES
#The function value (fv) table method:we will use computed inhomogenous cross K function

# fv object
#Kij_Plot5 is the fv object computed above (inhomogenous cross K from Co to AsT)
#method "b" apply smoothing to Y(r) = K(r)/(2 * pi * r) constraining Y(0) = 0, 
#estimate the derivative of Y, and solve;
#method "b" used due to likely inhibition at small distances of r spar=0.8,

gij_Plot5 <- pcf(Kij_Plot5, spar=0.8,  method="b") 
plot(gij_Plot5, main="Plot5 Co3 - AsT g(r) Pair correlation function", xlab="Distance(m)", legendpos="bottom")
write.csv(gij_Plot5, file = "pcf_Plot5_Co3_to_AsT.csv", row.names=TRUE)


#------------------------------------------------

#PLOT 21


Plot21_Data = read.csv("Plot21_tclass4.csv")
Plot21_Data$TTP = as.factor(Plot21_Data$TTP)

Plot21_T4  <- readShapePoints("Plot21_tclass4.shp")
SP <- as(Plot21_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot21_Data[5]

# Get Plot211 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot21.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot21 windown (polygon) to trees in Plot21
P$window <- W
plot(P, pch=16, main="Plot21 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot21 All Trees")
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
L_AsT_Plot21 <- Ldot.inhom(P, "AsT", lambdaX=fit, update=FALSE)
plot(Lij_Plot21, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot21 Associated Trees", legendpos="bottom")
L_AsT_Plot21$iso  <- L_AsT_Plot21$iso  - L_AsT_Plot21$r
L_AsT_Plot21$theo <- L_AsT_Plot21$theo - L_AsT_Plot21$r
write.csv(L_AsT_Plot21, file = "Li_AsT_Plot21.csv", row.names=TRUE)

#Li (nonhomogenous) for Cocoa Trees2 in stand 

L_Co2_Plot21 <- Ldot.inhom(P, "Co", lambdaX=fit, update=FALSE) #using fit computed for all tree

plot(Lij_Plot21, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot21 Cocoa Trees3", legendpos="bottom")
L_AsT_Plot21$iso  <- L_AsT_Plot21$iso  - L_AsT_Plot21$r
L_AsT_Plot21$theo <- L_AsT_Plot21$theo - L_AsT_Plot21$r
write.csv(L_AsT_Plot21, file = "Li_Co3_Plot21.csv", row.names=TRUE)


#Bivariate(cross) inhomogenous Kij(r) Ripley function
#The Kcross(r) test, which examines distribution a tree group relative to another group within stand

# method (3): fit parametric intensity model
fit <- ppm(P, ~marks * polynom(x,y,2))
# evaluate fitted intensities at data points
# (these are the intensities of the sub-processes of each type)
inten <- fitted(fit, dataonly=TRUE)
# split according to types of points
lambda2 <- split(inten, P$marks)
Kij_Plot21 <- Kcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))
plot(Kij_Plot21, xlab="Distance (m)", main="Kij (2nd order moment) Co3 to Ast")

write.csv(Kij_Plot21, file = "Kij_bivariate_inhom_Plot21_Co3-AsT.csv", row.names=TRUE)



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
Lij_Plot21 <- Lcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))

plot(Lij_Plot21, . - r ~ r, xlab="Distance (m)", ylab="Lij(r)", main="Lij inhom. Function, Plot21 Cocoa3-AsT", legendpos="bottom")


# We will need to modify the plot graphics parameters to accommodate the
# custom y-label; hence the call to `par`.
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Plot21, . -r ~ r, ylab=expression(hat("Lij")), xlab = "d (m)", main="Lij(r), Plot21 Cocoa Trees", legendmath=TRUE, legendpos="bottom")
#par(OP)

#We can make permanent changes to the Lij object's values as follows:
Lij_Plot21$iso  <- Lij_Plot21$iso  - Lij_Plot21$r
Lij_Plot21$theo <- Lij_Plot21$theo - Lij_Plot21$r

write.csv(Lij_Plot21, file = "Lij_Bivariate_inhom_Plot21_Co3_to_AsT.csv", row.names=TRUE)

#MONTE CARLO (MC) SIMULATION

#To test for clustering/inhibition in the presence of spatial inhomogeneity, 
#the null hypothesis should be an inhomogeneous Poisson process. 
#So we estimate the inhomogeneous intensity functions of the two types of points,
# and then generate simulated point patterns according to these intensities.

#First fit a Poisson point process model to the observed data
#fit2 <- ppm(P, ~marks * polynom(x,y,3))

#n <- 199 #number of iterations for Monte carlo simulation
#p  <- 0.01 # Desired p significance level to display

#Lij_Env_Plot21 <- envelope(fit2, Lcross.inhom, lambdaX=fit2, nsim=n, rank=(p * (n + 1)))
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Env_Plot21, . - r ~ r, ylab=expression(hat("L")),xlab="Distance (m)", main=paste("Lij(r) Plot21 Co3-AsT,  p= ", p),legendmath=TRUE)
#par(OP)
#write.csv(Lij_Env_Plot21, file = "Lij_199MCEnv_Plot21_Co3_to_AsT.csv", row.names=TRUE)



# COMPUTING THE g (Pair Correlation Fucntion) OF COCOA AND ASSOCIATED TREES
#The function value (fv) table method:we will use computed inhomogenous cross K function

# fv object
#Kij_Plot21 is the fv object computed above (inhomogenous cross K from Co to AsT)
#method "b" apply smoothing to Y(r) = K(r)/(2 * pi * r) constraining Y(0) = 0, 
#estimate the derivative of Y, and solve;
#method "b" used due to likely inhibition at small distances of r spar=0.8,

gij_Plot21 <- pcf(Kij_Plot21, spar=0.8,  method="b") 
plot(gij_Plot21, main="Plot21 Co3 - AsT g(r) Pair correlation function", xlab="Distance(m)", legendpos="bottom")
write.csv(gij_Plot21, file = "pcf_Plot21_Co3_to_AsT.csv", row.names=TRUE)


#--------------------------------------------------

#PLOT 8


Plot8_Data = read.csv("Plot8_tclass4.csv")
Plot8_Data$TTP = as.factor(Plot8_Data$TTP)

Plot8_T4  <- readShapePoints("Plot8_tclass4.shp")
SP <- as(Plot8_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot8_Data[5]

# Get Plot81 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot8.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot8 windown (polygon) to trees in Plot8
P$window <- W
plot(P, pch=16, main="Plot8 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot8 All Trees")
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
L_AsT_Plot8 <- Ldot.inhom(P, "AsT", lambdaX=fit, update=FALSE)
plot(Lij_Plot8, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot8 Associated Trees", legendpos="bottom")
L_AsT_Plot8$iso  <- L_AsT_Plot8$iso  - L_AsT_Plot8$r
L_AsT_Plot8$theo <- L_AsT_Plot8$theo - L_AsT_Plot8$r
write.csv(L_AsT_Plot8, file = "Li_AsT_Plot8.csv", row.names=TRUE)

#Li (nonhomogenous) for Cocoa Trees2 in stand 

L_Co2_Plot8 <- Ldot.inhom(P, "Co", lambdaX=fit, update=FALSE) #using fit computed for all tree

plot(Lij_Plot8, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot8 Cocoa Trees3", legendpos="bottom")
L_AsT_Plot8$iso  <- L_AsT_Plot8$iso  - L_AsT_Plot8$r
L_AsT_Plot8$theo <- L_AsT_Plot8$theo - L_AsT_Plot8$r
write.csv(L_AsT_Plot8, file = "Li_Co3_Plot8.csv", row.names=TRUE)


#Bivariate(cross) inhomogenous Kij(r) Ripley function
#The Kcross(r) test, which examines distribution a tree group relative to another group within stand

# method (3): fit parametric intensity model
fit <- ppm(P, ~marks * polynom(x,y,2))
# evaluate fitted intensities at data points
# (these are the intensities of the sub-processes of each type)
inten <- fitted(fit, dataonly=TRUE)
# split according to types of points
lambda2 <- split(inten, P$marks)
Kij_Plot8 <- Kcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))
plot(Kij_Plot8, xlab="Distance (m)", main="Kij (2nd order moment) Co3 to Ast")

write.csv(Kij_Plot8, file = "Kij_bivariate_inhom_Plot8_Co3-AsT.csv", row.names=TRUE)



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
Lij_Plot8 <- Lcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))

plot(Lij_Plot8, . - r ~ r, xlab="Distance (m)", ylab="Lij(r)", main="Lij inhom. Function, Plot8 Cocoa3-AsT", legendpos="bottom")


# We will need to modify the plot graphics parameters to accommodate the
# custom y-label; hence the call to `par`.
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Plot8, . -r ~ r, ylab=expression(hat("Lij")), xlab = "d (m)", main="Lij(r), Plot8 Cocoa Trees", legendmath=TRUE, legendpos="bottom")
#par(OP)

#We can make permanent changes to the Lij object's values as follows:
Lij_Plot8$iso  <- Lij_Plot8$iso  - Lij_Plot8$r
Lij_Plot8$theo <- Lij_Plot8$theo - Lij_Plot8$r

write.csv(Lij_Plot8, file = "Lij_Bivariate_inhom_Plot8_Co3_to_AsT.csv", row.names=TRUE)

#MONTE CARLO (MC) SIMULATION

#To test for clustering/inhibition in the presence of spatial inhomogeneity, 
#the null hypothesis should be an inhomogeneous Poisson process. 
#So we estimate the inhomogeneous intensity functions of the two types of points,
# and then generate simulated point patterns according to these intensities.

#First fit a Poisson point process model to the observed data
#fit2 <- ppm(P, ~marks * polynom(x,y,3))

#n <- 199 #number of iterations for Monte carlo simulation
#p  <- 0.01 # Desired p significance level to display

#Lij_Env_Plot8 <- envelope(fit2, Lcross.inhom, lambdaX=fit2, nsim=n, rank=(p * (n + 1)))
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Env_Plot8, . - r ~ r, ylab=expression(hat("L")),xlab="Distance (m)", main=paste("Lij(r) Plot8 Co3-AsT,  p= ", p),legendmath=TRUE)
#par(OP)
#write.csv(Lij_Env_Plot8, file = "Lij_199MCEnv_Plot8_Co3_to_AsT.csv", row.names=TRUE)



# COMPUTING THE g (Pair Correlation Fucntion) OF COCOA AND ASSOCIATED TREES
#The function value (fv) table method:we will use computed inhomogenous cross K function

# fv object
#Kij_Plot8 is the fv object computed above (inhomogenous cross K from Co to AsT)
#method "b" apply smoothing to Y(r) = K(r)/(2 * pi * r) constraining Y(0) = 0, 
#estimate the derivative of Y, and solve;
#method "b" used due to likely inhibition at small distances of r spar=0.8,

gij_Plot8 <- pcf(Kij_Plot8, spar=0.8,  method="b") 
plot(gij_Plot8, main="Plot8 Co3 - AsT g(r) Pair correlation function", xlab="Distance(m)", legendpos="bottom")
write.csv(gij_Plot8, file = "pcf_Plot8_Co3_to_AsT.csv", row.names=TRUE)


#------------------------------------------------

#PLOT 34


Plot34_Data = read.csv("Plot34_tclass4.csv")
Plot34_Data$TTP = as.factor(Plot34_Data$TTP)

Plot34_T4  <- readShapePoints("Plot34_tclass4.shp")
SP <- as(Plot34_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot34_Data[5]

# Get Plot341 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot34.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot34 windown (polygon) to trees in Plot34
P$window <- W
plot(P, pch=16, main="Plot34 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot34 All Trees")
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
L_AsT_Plot34 <- Ldot.inhom(P, "AsT", lambdaX=fit, update=FALSE)
plot(Lij_Plot34, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot34 Associated Trees", legendpos="bottom")
L_AsT_Plot34$iso  <- L_AsT_Plot34$iso  - L_AsT_Plot34$r
L_AsT_Plot34$theo <- L_AsT_Plot34$theo - L_AsT_Plot34$r
write.csv(L_AsT_Plot34, file = "Li_AsT_Plot34.csv", row.names=TRUE)

#Li (nonhomogenous) for Cocoa Trees2 in stand 

L_Co2_Plot34 <- Ldot.inhom(P, "Co", lambdaX=fit, update=FALSE) #using fit computed for all tree

plot(Lij_Plot34, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot34 Cocoa Trees3", legendpos="bottom")
L_AsT_Plot34$iso  <- L_AsT_Plot34$iso  - L_AsT_Plot34$r
L_AsT_Plot34$theo <- L_AsT_Plot34$theo - L_AsT_Plot34$r
write.csv(L_AsT_Plot34, file = "Li_Co3_Plot34.csv", row.names=TRUE)


#Bivariate(cross) inhomogenous Kij(r) Ripley function
#The Kcross(r) test, which examines distribution a tree group relative to another group within stand

# method (3): fit parametric intensity model
fit <- ppm(P, ~marks * polynom(x,y,2))
# evaluate fitted intensities at data points
# (these are the intensities of the sub-processes of each type)
inten <- fitted(fit, dataonly=TRUE)
# split according to types of points
lambda2 <- split(inten, P$marks)
Kij_Plot34 <- Kcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))
plot(Kij_Plot34, xlab="Distance (m)", main="Kij (2nd order moment) Co3 to Ast")

write.csv(Kij_Plot34, file = "Kij_bivariate_inhom_Plot34_Co3-AsT.csv", row.names=TRUE)



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
Lij_Plot34 <- Lcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))

plot(Lij_Plot34, . - r ~ r, xlab="Distance (m)", ylab="Lij(r)", main="Lij inhom. Function, Plot34 Cocoa3-AsT", legendpos="bottom")


# We will need to modify the plot graphics parameters to accommodate the
# custom y-label; hence the call to `par`.
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Plot34, . -r ~ r, ylab=expression(hat("Lij")), xlab = "d (m)", main="Lij(r), Plot34 Cocoa Trees", legendmath=TRUE, legendpos="bottom")
#par(OP)

#We can make permanent changes to the Lij object's values as follows:
Lij_Plot34$iso  <- Lij_Plot34$iso  - Lij_Plot34$r
Lij_Plot34$theo <- Lij_Plot34$theo - Lij_Plot34$r

write.csv(Lij_Plot34, file = "Lij_Bivariate_inhom_Plot34_Co3_to_AsT.csv", row.names=TRUE)

#MONTE CARLO (MC) SIMULATION

#To test for clustering/inhibition in the presence of spatial inhomogeneity, 
#the null hypothesis should be an inhomogeneous Poisson process. 
#So we estimate the inhomogeneous intensity functions of the two types of points,
# and then generate simulated point patterns according to these intensities.

#First fit a Poisson point process model to the observed data
#fit2 <- ppm(P, ~marks * polynom(x,y,3))

#n <- 199 #number of iterations for Monte carlo simulation
#p  <- 0.01 # Desired p significance level to display

#Lij_Env_Plot34 <- envelope(fit2, Lcross.inhom, lambdaX=fit2, nsim=n, rank=(p * (n + 1)))
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Env_Plot34, . - r ~ r, ylab=expression(hat("L")),xlab="Distance (m)", main=paste("Lij(r) Plot34 Co3-AsT,  p= ", p),legendmath=TRUE)
#par(OP)
#write.csv(Lij_Env_Plot34, file = "Lij_199MCEnv_Plot34_Co3_to_AsT.csv", row.names=TRUE)



# COMPUTING THE g (Pair Correlation Fucntion) OF COCOA AND ASSOCIATED TREES
#The function value (fv) table method:we will use computed inhomogenous cross K function

# fv object
#Kij_Plot34 is the fv object computed above (inhomogenous cross K from Co to AsT)
#method "b" apply smoothing to Y(r) = K(r)/(2 * pi * r) constraining Y(0) = 0, 
#estimate the derivative of Y, and solve;
#method "b" used due to likely inhibition at small distances of r spar=0.8,

gij_Plot34 <- pcf(Kij_Plot34, spar=0.8,  method="b") 
plot(gij_Plot34, main="Plot34 Co3 - AsT g(r) Pair correlation function", xlab="Distance(m)", legendpos="bottom")
write.csv(gij_Plot34, file = "pcf_Plot34_Co3_to_AsT.csv", row.names=TRUE)


#--------------------------------------------------

#PLOT 29


Plot29_Data = read.csv("Plot29_tclass4.csv")
Plot29_Data$TTP = as.factor(Plot29_Data$TTP)

Plot29_T4  <- readShapePoints("Plot29_tclass4.shp")
SP <- as(Plot29_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot29_Data[5]

# Get Plot291 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot29.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot29 windown (polygon) to trees in Plot29
P$window <- W
plot(P, pch=16, main="Plot29 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot29 All Trees")
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
L_AsT_Plot29 <- Ldot.inhom(P, "AsT", lambdaX=fit, update=FALSE)
plot(Lij_Plot29, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot29 Associated Trees", legendpos="bottom")
L_AsT_Plot29$iso  <- L_AsT_Plot29$iso  - L_AsT_Plot29$r
L_AsT_Plot29$theo <- L_AsT_Plot29$theo - L_AsT_Plot29$r
write.csv(L_AsT_Plot29, file = "Li_AsT_Plot29.csv", row.names=TRUE)

#Li (nonhomogenous) for Cocoa Trees2 in stand 

L_Co2_Plot29 <- Ldot.inhom(P, "Co", lambdaX=fit, update=FALSE) #using fit computed for all tree

plot(Lij_Plot29, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot29 Cocoa Trees3", legendpos="bottom")
L_AsT_Plot29$iso  <- L_AsT_Plot29$iso  - L_AsT_Plot29$r
L_AsT_Plot29$theo <- L_AsT_Plot29$theo - L_AsT_Plot29$r
write.csv(L_AsT_Plot29, file = "Li_Co3_Plot29.csv", row.names=TRUE)


#Bivariate(cross) inhomogenous Kij(r) Ripley function
#The Kcross(r) test, which examines distribution a tree group relative to another group within stand

# method (3): fit parametric intensity model
fit <- ppm(P, ~marks * polynom(x,y,2))
# evaluate fitted intensities at data points
# (these are the intensities of the sub-processes of each type)
inten <- fitted(fit, dataonly=TRUE)
# split according to types of points
lambda2 <- split(inten, P$marks)
Kij_Plot29 <- Kcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))
plot(Kij_Plot29, xlab="Distance (m)", main="Kij (2nd order moment) Co3 to Ast")

write.csv(Kij_Plot29, file = "Kij_bivariate_inhom_Plot29_Co3-AsT.csv", row.names=TRUE)



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
Lij_Plot29 <- Lcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))

plot(Lij_Plot29, . - r ~ r, xlab="Distance (m)", ylab="Lij(r)", main="Lij inhom. Function, Plot29 Cocoa3-AsT", legendpos="bottom")


# We will need to modify the plot graphics parameters to accommodate the
# custom y-label; hence the call to `par`.
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Plot29, . -r ~ r, ylab=expression(hat("Lij")), xlab = "d (m)", main="Lij(r), Plot29 Cocoa Trees", legendmath=TRUE, legendpos="bottom")
#par(OP)

#We can make permanent changes to the Lij object's values as follows:
Lij_Plot29$iso  <- Lij_Plot29$iso  - Lij_Plot29$r
Lij_Plot29$theo <- Lij_Plot29$theo - Lij_Plot29$r

write.csv(Lij_Plot29, file = "Lij_Bivariate_inhom_Plot29_Co3_to_AsT.csv", row.names=TRUE)

#MONTE CARLO (MC) SIMULATION

#To test for clustering/inhibition in the presence of spatial inhomogeneity, 
#the null hypothesis should be an inhomogeneous Poisson process. 
#So we estimate the inhomogeneous intensity functions of the two types of points,
# and then generate simulated point patterns according to these intensities.

#First fit a Poisson point process model to the observed data
#fit2 <- ppm(P, ~marks * polynom(x,y,3))

#n <- 199 #number of iterations for Monte carlo simulation
#p  <- 0.01 # Desired p significance level to display

#Lij_Env_Plot29 <- envelope(fit2, Lcross.inhom, lambdaX=fit2, nsim=n, rank=(p * (n + 1)))
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Env_Plot29, . - r ~ r, ylab=expression(hat("L")),xlab="Distance (m)", main=paste("Lij(r) Plot29 Co3-AsT,  p= ", p),legendmath=TRUE)
#par(OP)
#write.csv(Lij_Env_Plot29, file = "Lij_199MCEnv_Plot29_Co3_to_AsT.csv", row.names=TRUE)



# COMPUTING THE g (Pair Correlation Fucntion) OF COCOA AND ASSOCIATED TREES
#The function value (fv) table method:we will use computed inhomogenous cross K function

# fv object
#Kij_Plot29 is the fv object computed above (inhomogenous cross K from Co to AsT)
#method "b" apply smoothing to Y(r) = K(r)/(2 * pi * r) constraining Y(0) = 0, 
#estimate the derivative of Y, and solve;
#method "b" used due to likely inhibition at small distances of r spar=0.8,

gij_Plot29 <- pcf(Kij_Plot29, spar=0.8,  method="b") 
plot(gij_Plot29, main="Plot29 Co3 - AsT g(r) Pair correlation function", xlab="Distance(m)", legendpos="bottom")
write.csv(gij_Plot29, file = "pcf_Plot29_Co3_to_AsT.csv", row.names=TRUE)


#--------------------------------------------------

#PLOT 42


Plot42_Data = read.csv("Plot42_tclass4.csv")
Plot42_Data$TTP = as.factor(Plot42_Data$TTP)

Plot42_T4  <- readShapePoints("Plot42_tclass4.shp")
SP <- as(Plot42_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot42_Data[5]

# Get Plot421 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot42.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot42 windown (polygon) to trees in Plot42
P$window <- W
plot(P, pch=16, main="Plot42 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot42 All Trees")
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
L_AsT_Plot42 <- Ldot.inhom(P, "AsT", lambdaX=fit, update=FALSE)
plot(Lij_Plot42, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot42 Associated Trees", legendpos="bottom")
L_AsT_Plot42$iso  <- L_AsT_Plot42$iso  - L_AsT_Plot42$r
L_AsT_Plot42$theo <- L_AsT_Plot42$theo - L_AsT_Plot42$r
write.csv(L_AsT_Plot42, file = "Li_AsT_Plot42.csv", row.names=TRUE)

#Li (nonhomogenous) for Cocoa Trees2 in stand 

L_Co2_Plot42 <- Ldot.inhom(P, "Co", lambdaX=fit, update=FALSE) #using fit computed for all tree

plot(Lij_Plot42, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot42 Cocoa Trees3", legendpos="bottom")
L_AsT_Plot42$iso  <- L_AsT_Plot42$iso  - L_AsT_Plot42$r
L_AsT_Plot42$theo <- L_AsT_Plot42$theo - L_AsT_Plot42$r
write.csv(L_AsT_Plot42, file = "Li_Co3_Plot42.csv", row.names=TRUE)


#Bivariate(cross) inhomogenous Kij(r) Ripley function
#The Kcross(r) test, which examines distribution a tree group relative to another group within stand

# method (3): fit parametric intensity model
fit <- ppm(P, ~marks * polynom(x,y,2))
# evaluate fitted intensities at data points
# (these are the intensities of the sub-processes of each type)
inten <- fitted(fit, dataonly=TRUE)
# split according to types of points
lambda2 <- split(inten, P$marks)
Kij_Plot42 <- Kcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))
plot(Kij_Plot42, xlab="Distance (m)", main="Kij (2nd order moment) Co3 to Ast")

write.csv(Kij_Plot42, file = "Kij_bivariate_inhom_Plot42_Co3-AsT.csv", row.names=TRUE)



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
Lij_Plot42 <- Lcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))

plot(Lij_Plot42, . - r ~ r, xlab="Distance (m)", ylab="Lij(r)", main="Lij inhom. Function, Plot42 Cocoa3-AsT", legendpos="bottom")


# We will need to modify the plot graphics parameters to accommodate the
# custom y-label; hence the call to `par`.
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Plot42, . -r ~ r, ylab=expression(hat("Lij")), xlab = "d (m)", main="Lij(r), Plot42 Cocoa Trees", legendmath=TRUE, legendpos="bottom")
#par(OP)

#We can make permanent changes to the Lij object's values as follows:
Lij_Plot42$iso  <- Lij_Plot42$iso  - Lij_Plot42$r
Lij_Plot42$theo <- Lij_Plot42$theo - Lij_Plot42$r

write.csv(Lij_Plot42, file = "Lij_Bivariate_inhom_Plot42_Co3_to_AsT.csv", row.names=TRUE)

#MONTE CARLO (MC) SIMULATION

#To test for clustering/inhibition in the presence of spatial inhomogeneity, 
#the null hypothesis should be an inhomogeneous Poisson process. 
#So we estimate the inhomogeneous intensity functions of the two types of points,
# and then generate simulated point patterns according to these intensities.

#First fit a Poisson point process model to the observed data
#fit2 <- ppm(P, ~marks * polynom(x,y,3))

#n <- 199 #number of iterations for Monte carlo simulation
#p  <- 0.01 # Desired p significance level to display

#Lij_Env_Plot42 <- envelope(fit2, Lcross.inhom, lambdaX=fit2, nsim=n, rank=(p * (n + 1)))
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Env_Plot42, . - r ~ r, ylab=expression(hat("L")),xlab="Distance (m)", main=paste("Lij(r) Plot42 Co3-AsT,  p= ", p),legendmath=TRUE)
#par(OP)
#write.csv(Lij_Env_Plot42, file = "Lij_199MCEnv_Plot42_Co3_to_AsT.csv", row.names=TRUE)



# COMPUTING THE g (Pair Correlation Fucntion) OF COCOA AND ASSOCIATED TREES
#The function value (fv) table method:we will use computed inhomogenous cross K function

# fv object
#Kij_Plot42 is the fv object computed above (inhomogenous cross K from Co to AsT)
#method "b" apply smoothing to Y(r) = K(r)/(2 * pi * r) constraining Y(0) = 0, 
#estimate the derivative of Y, and solve;
#method "b" used due to likely inhibition at small distances of r spar=0.8,

gij_Plot42 <- pcf(Kij_Plot42, spar=0.8,  method="b") 
plot(gij_Plot42, main="Plot42 Co3 - AsT g(r) Pair correlation function", xlab="Distance(m)", legendpos="bottom")
write.csv(gij_Plot42, file = "pcf_Plot42_Co3_to_AsT.csv", row.names=TRUE)


#--------------------------------------------------

#PLOT 42


Plot3_Data = read.csv("Plot3_tclass4.csv")
Plot3_Data$TTP = as.factor(Plot3_Data$TTP)

Plot3_T4  <- readShapePoints("Plot3_tclass4.shp")
SP <- as(Plot3_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot3_Data[5]

# Get Plot31 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot3.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot3 windown (polygon) to trees in Plot3
P$window <- W
plot(P, pch=16, main="Plot3 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot3 All Trees")
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
L_AsT_Plot3 <- Ldot.inhom(P, "AsT", lambdaX=fit, update=FALSE)
plot(Lij_Plot3, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot3 Associated Trees", legendpos="bottom")
L_AsT_Plot3$iso  <- L_AsT_Plot3$iso  - L_AsT_Plot3$r
L_AsT_Plot3$theo <- L_AsT_Plot3$theo - L_AsT_Plot3$r
write.csv(L_AsT_Plot3, file = "Li_AsT_Plot3.csv", row.names=TRUE)

#Li (nonhomogenous) for Cocoa Trees2 in stand 

L_Co2_Plot3 <- Ldot.inhom(P, "Co", lambdaX=fit, update=FALSE) #using fit computed for all tree

plot(Lij_Plot3, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot3 Cocoa Trees3", legendpos="bottom")
L_AsT_Plot3$iso  <- L_AsT_Plot3$iso  - L_AsT_Plot3$r
L_AsT_Plot3$theo <- L_AsT_Plot3$theo - L_AsT_Plot3$r
write.csv(L_AsT_Plot3, file = "Li_Co3_Plot3.csv", row.names=TRUE)


#Bivariate(cross) inhomogenous Kij(r) Ripley function
#The Kcross(r) test, which examines distribution a tree group relative to another group within stand

# method (3): fit parametric intensity model
fit <- ppm(P, ~marks * polynom(x,y,2))
# evaluate fitted intensities at data points
# (these are the intensities of the sub-processes of each type)
inten <- fitted(fit, dataonly=TRUE)
# split according to types of points
lambda2 <- split(inten, P$marks)
Kij_Plot3 <- Kcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))
plot(Kij_Plot3, xlab="Distance (m)", main="Kij (2nd order moment) Co3 to Ast")

write.csv(Kij_Plot3, file = "Kij_bivariate_inhom_Plot3_Co3-AsT.csv", row.names=TRUE)



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
Lij_Plot3 <- Lcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))

plot(Lij_Plot3, . - r ~ r, xlab="Distance (m)", ylab="Lij(r)", main="Lij inhom. Function, Plot3 Cocoa3-AsT", legendpos="bottom")


# We will need to modify the plot graphics parameters to accommodate the
# custom y-label; hence the call to `par`.
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Plot3, . -r ~ r, ylab=expression(hat("Lij")), xlab = "d (m)", main="Lij(r), Plot3 Cocoa Trees", legendmath=TRUE, legendpos="bottom")
#par(OP)

#We can make permanent changes to the Lij object's values as follows:
Lij_Plot3$iso  <- Lij_Plot3$iso  - Lij_Plot3$r
Lij_Plot3$theo <- Lij_Plot3$theo - Lij_Plot3$r

write.csv(Lij_Plot3, file = "Lij_Bivariate_inhom_Plot3_Co3_to_AsT.csv", row.names=TRUE)

#MONTE CARLO (MC) SIMULATION

#To test for clustering/inhibition in the presence of spatial inhomogeneity, 
#the null hypothesis should be an inhomogeneous Poisson process. 
#So we estimate the inhomogeneous intensity functions of the two types of points,
# and then generate simulated point patterns according to these intensities.

#First fit a Poisson point process model to the observed data
#fit2 <- ppm(P, ~marks * polynom(x,y,3))

#n <- 199 #number of iterations for Monte carlo simulation
#p  <- 0.01 # Desired p significance level to display

#Lij_Env_Plot3 <- envelope(fit2, Lcross.inhom, lambdaX=fit2, nsim=n, rank=(p * (n + 1)))
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Env_Plot3, . - r ~ r, ylab=expression(hat("L")),xlab="Distance (m)", main=paste("Lij(r) Plot3 Co3-AsT,  p= ", p),legendmath=TRUE)
#par(OP)
#write.csv(Lij_Env_Plot3, file = "Lij_199MCEnv_Plot3_Co3_to_AsT.csv", row.names=TRUE)



# COMPUTING THE g (Pair Correlation Fucntion) OF COCOA AND ASSOCIATED TREES
#The function value (fv) table method:we will use computed inhomogenous cross K function

# fv object
#Kij_Plot3 is the fv object computed above (inhomogenous cross K from Co to AsT)
#method "b" apply smoothing to Y(r) = K(r)/(2 * pi * r) constraining Y(0) = 0, 
#estimate the derivative of Y, and solve;
#method "b" used due to likely inhibition at small distances of r spar=0.8,

gij_Plot3 <- pcf(Kij_Plot3, spar=0.8,  method="b") 
plot(gij_Plot3, main="Plot3 Co3 - AsT g(r) Pair correlation function", xlab="Distance(m)", legendpos="bottom")
write.csv(gij_Plot3, file = "pcf_Plot3_Co3_to_AsT.csv", row.names=TRUE)


#------------------------------------------------

#PLOT 10


Plot10_Data = read.csv("Plot10_tclass4.csv")
Plot10_Data$TTP = as.factor(Plot10_Data$TTP)

Plot10_T4  <- readShapePoints("Plot10_tclass4.shp")
SP <- as(Plot10_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot10_Data[5]

# Get Plot101 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot10.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot10 windown (polygon) to trees in Plot10
P$window <- W
plot(P, pch=16, main="Plot10 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot10 All Trees")
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
L_AsT_Plot10 <- Ldot.inhom(P, "AsT", lambdaX=fit, update=FALSE)
plot(Lij_Plot10, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot10 Associated Trees", legendpos="bottom")
L_AsT_Plot10$iso  <- L_AsT_Plot10$iso  - L_AsT_Plot10$r
L_AsT_Plot10$theo <- L_AsT_Plot10$theo - L_AsT_Plot10$r
write.csv(L_AsT_Plot10, file = "Li_AsT_Plot10.csv", row.names=TRUE)

#Li (nonhomogenous) for Cocoa Trees2 in stand 

L_Co2_Plot10 <- Ldot.inhom(P, "Co", lambdaX=fit, update=FALSE) #using fit computed for all tree

plot(Lij_Plot10, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot10 Cocoa Trees3", legendpos="bottom")
L_AsT_Plot10$iso  <- L_AsT_Plot10$iso  - L_AsT_Plot10$r
L_AsT_Plot10$theo <- L_AsT_Plot10$theo - L_AsT_Plot10$r
write.csv(L_AsT_Plot10, file = "Li_Co3_Plot10.csv", row.names=TRUE)


#Bivariate(cross) inhomogenous Kij(r) Ripley function
#The Kcross(r) test, which examines distribution a tree group relative to another group within stand

# method (3): fit parametric intensity model
fit <- ppm(P, ~marks * polynom(x,y,2))
# evaluate fitted intensities at data points
# (these are the intensities of the sub-processes of each type)
inten <- fitted(fit, dataonly=TRUE)
# split according to types of points
lambda2 <- split(inten, P$marks)
Kij_Plot10 <- Kcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))
plot(Kij_Plot10, xlab="Distance (m)", main="Kij (2nd order moment) Co3 to Ast")

write.csv(Kij_Plot10, file = "Kij_bivariate_inhom_Plot10_Co3-AsT.csv", row.names=TRUE)



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
Lij_Plot10 <- Lcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))

plot(Lij_Plot10, . - r ~ r, xlab="Distance (m)", ylab="Lij(r)", main="Lij inhom. Function, Plot10 Cocoa3-AsT", legendpos="bottom")


# We will need to modify the plot graphics parameters to accommodate the
# custom y-label; hence the call to `par`.
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Plot10, . -r ~ r, ylab=expression(hat("Lij")), xlab = "d (m)", main="Lij(r), Plot10 Cocoa Trees", legendmath=TRUE, legendpos="bottom")
#par(OP)

#We can make permanent changes to the Lij object's values as follows:
Lij_Plot10$iso  <- Lij_Plot10$iso  - Lij_Plot10$r
Lij_Plot10$theo <- Lij_Plot10$theo - Lij_Plot10$r

write.csv(Lij_Plot10, file = "Lij_Bivariate_inhom_Plot10_Co3_to_AsT.csv", row.names=TRUE)

#MONTE CARLO (MC) SIMULATION

#To test for clustering/inhibition in the presence of spatial inhomogeneity, 
#the null hypothesis should be an inhomogeneous Poisson process. 
#So we estimate the inhomogeneous intensity functions of the two types of points,
# and then generate simulated point patterns according to these intensities.

#First fit a Poisson point process model to the observed data
#fit2 <- ppm(P, ~marks * polynom(x,y,3))

#n <- 199 #number of iterations for Monte carlo simulation
#p  <- 0.01 # Desired p significance level to display

#Lij_Env_Plot10 <- envelope(fit2, Lcross.inhom, lambdaX=fit2, nsim=n, rank=(p * (n + 1)))
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Env_Plot10, . - r ~ r, ylab=expression(hat("L")),xlab="Distance (m)", main=paste("Lij(r) Plot10 Co3-AsT,  p= ", p),legendmath=TRUE)
#par(OP)
#write.csv(Lij_Env_Plot10, file = "Lij_199MCEnv_Plot10_Co3_to_AsT.csv", row.names=TRUE)



# COMPUTING THE g (Pair Correlation Fucntion) OF COCOA AND ASSOCIATED TREES
#The function value (fv) table method:we will use computed inhomogenous cross K function

# fv object
#Kij_Plot10 is the fv object computed above (inhomogenous cross K from Co to AsT)
#method "b" apply smoothing to Y(r) = K(r)/(2 * pi * r) constraining Y(0) = 0, 
#estimate the derivative of Y, and solve;
#method "b" used due to likely inhibition at small distances of r spar=0.8,

gij_Plot10 <- pcf(Kij_Plot10, spar=0.8,  method="b") 
plot(gij_Plot10, main="Plot10 Co3 - AsT g(r) Pair correlation function", xlab="Distance(m)", legendpos="bottom")
write.csv(gij_Plot10, file = "pcf_Plot10_Co3_to_AsT.csv", row.names=TRUE)


#--------------------------------------------------

#PLOT 11


Plot11_Data = read.csv("Plot11_tclass4.csv")
Plot11_Data$TTP = as.factor(Plot11_Data$TTP)

Plot11_T4  <- readShapePoints("Plot11_tclass4.shp")
SP <- as(Plot11_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot11_Data[5]

# Get Plot111 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot11.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot11 windown (polygon) to trees in Plot11
P$window <- W
plot(P, pch=16, main="Plot11 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot11 All Trees")
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
L_AsT_Plot11 <- Ldot.inhom(P, "AsT", lambdaX=fit, update=FALSE)
plot(Lij_Plot11, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot11 Associated Trees", legendpos="bottom")
L_AsT_Plot11$iso  <- L_AsT_Plot11$iso  - L_AsT_Plot11$r
L_AsT_Plot11$theo <- L_AsT_Plot11$theo - L_AsT_Plot11$r
write.csv(L_AsT_Plot11, file = "Li_AsT_Plot11.csv", row.names=TRUE)

#Li (nonhomogenous) for Cocoa Trees2 in stand 

L_Co2_Plot11 <- Ldot.inhom(P, "Co", lambdaX=fit, update=FALSE) #using fit computed for all tree

plot(Lij_Plot11, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot11 Cocoa Trees3", legendpos="bottom")
L_AsT_Plot11$iso  <- L_AsT_Plot11$iso  - L_AsT_Plot11$r
L_AsT_Plot11$theo <- L_AsT_Plot11$theo - L_AsT_Plot11$r
write.csv(L_AsT_Plot11, file = "Li_Co3_Plot11.csv", row.names=TRUE)


#Bivariate(cross) inhomogenous Kij(r) Ripley function
#The Kcross(r) test, which examines distribution a tree group relative to another group within stand

# method (3): fit parametric intensity model
fit <- ppm(P, ~marks * polynom(x,y,2))
# evaluate fitted intensities at data points
# (these are the intensities of the sub-processes of each type)
inten <- fitted(fit, dataonly=TRUE)
# split according to types of points
lambda2 <- split(inten, P$marks)
Kij_Plot11 <- Kcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))
plot(Kij_Plot11, xlab="Distance (m)", main="Kij (2nd order moment) Co3 to Ast")

write.csv(Kij_Plot11, file = "Kij_bivariate_inhom_Plot11_Co3-AsT.csv", row.names=TRUE)



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
Lij_Plot11 <- Lcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))

plot(Lij_Plot11, . - r ~ r, xlab="Distance (m)", ylab="Lij(r)", main="Lij inhom. Function, Plot11 Cocoa3-AsT", legendpos="bottom")


# We will need to modify the plot graphics parameters to accommodate the
# custom y-label; hence the call to `par`.
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Plot11, . -r ~ r, ylab=expression(hat("Lij")), xlab = "d (m)", main="Lij(r), Plot11 Cocoa Trees", legendmath=TRUE, legendpos="bottom")
#par(OP)

#We can make permanent changes to the Lij object's values as follows:
Lij_Plot11$iso  <- Lij_Plot11$iso  - Lij_Plot11$r
Lij_Plot11$theo <- Lij_Plot11$theo - Lij_Plot11$r

write.csv(Lij_Plot11, file = "Lij_Bivariate_inhom_Plot11_Co3_to_AsT.csv", row.names=TRUE)

#MONTE CARLO (MC) SIMULATION

#To test for clustering/inhibition in the presence of spatial inhomogeneity, 
#the null hypothesis should be an inhomogeneous Poisson process. 
#So we estimate the inhomogeneous intensity functions of the two types of points,
# and then generate simulated point patterns according to these intensities.

#First fit a Poisson point process model to the observed data
#fit2 <- ppm(P, ~marks * polynom(x,y,3))

#n <- 199 #number of iterations for Monte carlo simulation
#p  <- 0.01 # Desired p significance level to display

#Lij_Env_Plot11 <- envelope(fit2, Lcross.inhom, lambdaX=fit2, nsim=n, rank=(p * (n + 1)))
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Env_Plot11, . - r ~ r, ylab=expression(hat("L")),xlab="Distance (m)", main=paste("Lij(r) Plot11 Co3-AsT,  p= ", p),legendmath=TRUE)
#par(OP)
#write.csv(Lij_Env_Plot11, file = "Lij_199MCEnv_Plot11_Co3_to_AsT.csv", row.names=TRUE)



# COMPUTING THE g (Pair Correlation Fucntion) OF COCOA AND ASSOCIATED TREES
#The function value (fv) table method:we will use computed inhomogenous cross K function

# fv object
#Kij_Plot11 is the fv object computed above (inhomogenous cross K from Co to AsT)
#method "b" apply smoothing to Y(r) = K(r)/(2 * pi * r) constraining Y(0) = 0, 
#estimate the derivative of Y, and solve;
#method "b" used due to likely inhibition at small distances of r spar=0.8,

gij_Plot11 <- pcf(Kij_Plot11, spar=0.8,  method="b") 
plot(gij_Plot11, main="Plot11 Co3 - AsT g(r) Pair correlation function", xlab="Distance(m)", legendpos="bottom")
write.csv(gij_Plot11, file = "pcf_Plot11_Co3_to_AsT.csv", row.names=TRUE)


#--------------------------------------------------

#PLOT 12


Plot12_Data = read.csv("Plot12_tclass4.csv")
Plot12_Data$TTP = as.factor(Plot12_Data$TTP)

Plot12_T4  <- readShapePoints("Plot12_tclass4.shp")
SP <- as(Plot12_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot12_Data[5]

# Get Plot121 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot12.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot12 windown (polygon) to trees in Plot12
P$window <- W
plot(P, pch=16, main="Plot12 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot12 All Trees")
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
L_AsT_Plot12 <- Ldot.inhom(P, "AsT", lambdaX=fit, update=FALSE)
plot(Lij_Plot12, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot12 Associated Trees", legendpos="bottom")
L_AsT_Plot12$iso  <- L_AsT_Plot12$iso  - L_AsT_Plot12$r
L_AsT_Plot12$theo <- L_AsT_Plot12$theo - L_AsT_Plot12$r
write.csv(L_AsT_Plot12, file = "Li_AsT_Plot12.csv", row.names=TRUE)

#Li (nonhomogenous) for Cocoa Trees2 in stand 

L_Co2_Plot12 <- Ldot.inhom(P, "Co", lambdaX=fit, update=FALSE) #using fit computed for all tree

plot(Lij_Plot12, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot12 Cocoa Trees3", legendpos="bottom")
L_AsT_Plot12$iso  <- L_AsT_Plot12$iso  - L_AsT_Plot12$r
L_AsT_Plot12$theo <- L_AsT_Plot12$theo - L_AsT_Plot12$r
write.csv(L_AsT_Plot12, file = "Li_Co3_Plot12.csv", row.names=TRUE)


#Bivariate(cross) inhomogenous Kij(r) Ripley function
#The Kcross(r) test, which examines distribution a tree group relative to another group within stand

# method (3): fit parametric intensity model
fit <- ppm(P, ~marks * polynom(x,y,2))
# evaluate fitted intensities at data points
# (these are the intensities of the sub-processes of each type)
inten <- fitted(fit, dataonly=TRUE)
# split according to types of points
lambda2 <- split(inten, P$marks)
Kij_Plot12 <- Kcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))
plot(Kij_Plot12, xlab="Distance (m)", main="Kij (2nd order moment) Co3 to Ast")

write.csv(Kij_Plot12, file = "Kij_bivariate_inhom_Plot12_Co3-AsT.csv", row.names=TRUE)



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
Lij_Plot12 <- Lcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))

plot(Lij_Plot12, . - r ~ r, xlab="Distance (m)", ylab="Lij(r)", main="Lij inhom. Function, Plot12 Cocoa3-AsT", legendpos="bottom")


# We will need to modify the plot graphics parameters to accommodate the
# custom y-label; hence the call to `par`.
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Plot12, . -r ~ r, ylab=expression(hat("Lij")), xlab = "d (m)", main="Lij(r), Plot12 Cocoa Trees", legendmath=TRUE, legendpos="bottom")
#par(OP)

#We can make permanent changes to the Lij object's values as follows:
Lij_Plot12$iso  <- Lij_Plot12$iso  - Lij_Plot12$r
Lij_Plot12$theo <- Lij_Plot12$theo - Lij_Plot12$r

write.csv(Lij_Plot12, file = "Lij_Bivariate_inhom_Plot12_Co3_to_AsT.csv", row.names=TRUE)

#MONTE CARLO (MC) SIMULATION

#To test for clustering/inhibition in the presence of spatial inhomogeneity, 
#the null hypothesis should be an inhomogeneous Poisson process. 
#So we estimate the inhomogeneous intensity functions of the two types of points,
# and then generate simulated point patterns according to these intensities.

#First fit a Poisson point process model to the observed data
#fit2 <- ppm(P, ~marks * polynom(x,y,3))

#n <- 199 #number of iterations for Monte carlo simulation
#p  <- 0.01 # Desired p significance level to display

#Lij_Env_Plot12 <- envelope(fit2, Lcross.inhom, lambdaX=fit2, nsim=n, rank=(p * (n + 1)))
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Env_Plot12, . - r ~ r, ylab=expression(hat("L")),xlab="Distance (m)", main=paste("Lij(r) Plot12 Co3-AsT,  p= ", p),legendmath=TRUE)
#par(OP)
#write.csv(Lij_Env_Plot12, file = "Lij_199MCEnv_Plot12_Co3_to_AsT.csv", row.names=TRUE)



# COMPUTING THE g (Pair Correlation Fucntion) OF COCOA AND ASSOCIATED TREES
#The function value (fv) table method:we will use computed inhomogenous cross K function

# fv object
#Kij_Plot12 is the fv object computed above (inhomogenous cross K from Co to AsT)
#method "b" apply smoothing to Y(r) = K(r)/(2 * pi * r) constraining Y(0) = 0, 
#estimate the derivative of Y, and solve;
#method "b" used due to likely inhibition at small distances of r spar=0.8,

gij_Plot12 <- pcf(Kij_Plot12, spar=0.8,  method="b") 
plot(gij_Plot12, main="Plot12 Co3 - AsT g(r) Pair correlation function", xlab="Distance(m)", legendpos="bottom")
write.csv(gij_Plot12, file = "pcf_Plot12_Co3_to_AsT.csv", row.names=TRUE)


#--------------------------------------------------

#PLOT 2


Plot2_Data = read.csv("Plot2_tclass4.csv")
Plot2_Data$TTP = as.factor(Plot2_Data$TTP)

Plot2_T4  <- readShapePoints("Plot2_tclass4.shp")
SP <- as(Plot2_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot2_Data[5]

# Get Plot21 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot2.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot2 windown (polygon) to trees in Plot2
P$window <- W
plot(P, pch=16, main="Plot2 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot2 All Trees")
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
L_AsT_Plot2 <- Ldot.inhom(P, "AsT", lambdaX=fit, update=FALSE)
plot(Lij_Plot2, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot2 Associated Trees", legendpos="bottom")
L_AsT_Plot2$iso  <- L_AsT_Plot2$iso  - L_AsT_Plot2$r
L_AsT_Plot2$theo <- L_AsT_Plot2$theo - L_AsT_Plot2$r
write.csv(L_AsT_Plot2, file = "Li_AsT_Plot2.csv", row.names=TRUE)

#Li (nonhomogenous) for Cocoa Trees2 in stand 

L_Co2_Plot2 <- Ldot.inhom(P, "Co", lambdaX=fit, update=FALSE) #using fit computed for all tree

plot(Lij_Plot2, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot2 Cocoa Trees3", legendpos="bottom")
L_AsT_Plot2$iso  <- L_AsT_Plot2$iso  - L_AsT_Plot2$r
L_AsT_Plot2$theo <- L_AsT_Plot2$theo - L_AsT_Plot2$r
write.csv(L_AsT_Plot2, file = "Li_Co3_Plot2.csv", row.names=TRUE)


#Bivariate(cross) inhomogenous Kij(r) Ripley function
#The Kcross(r) test, which examines distribution a tree group relative to another group within stand

# method (3): fit parametric intensity model
fit <- ppm(P, ~marks * polynom(x,y,2))
# evaluate fitted intensities at data points
# (these are the intensities of the sub-processes of each type)
inten <- fitted(fit, dataonly=TRUE)
# split according to types of points
lambda2 <- split(inten, P$marks)
Kij_Plot2 <- Kcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))
plot(Kij_Plot2, xlab="Distance (m)", main="Kij (2nd order moment) Co3 to Ast")

write.csv(Kij_Plot2, file = "Kij_bivariate_inhom_Plot2_Co3-AsT.csv", row.names=TRUE)



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
Lij_Plot2 <- Lcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))

plot(Lij_Plot2, . - r ~ r, xlab="Distance (m)", ylab="Lij(r)", main="Lij inhom. Function, Plot2 Cocoa3-AsT", legendpos="bottom")


# We will need to modify the plot graphics parameters to accommodate the
# custom y-label; hence the call to `par`.
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Plot2, . -r ~ r, ylab=expression(hat("Lij")), xlab = "d (m)", main="Lij(r), Plot2 Cocoa Trees", legendmath=TRUE, legendpos="bottom")
#par(OP)

#We can make permanent changes to the Lij object's values as follows:
Lij_Plot2$iso  <- Lij_Plot2$iso  - Lij_Plot2$r
Lij_Plot2$theo <- Lij_Plot2$theo - Lij_Plot2$r

write.csv(Lij_Plot2, file = "Lij_Bivariate_inhom_Plot2_Co3_to_AsT.csv", row.names=TRUE)

#MONTE CARLO (MC) SIMULATION

#To test for clustering/inhibition in the presence of spatial inhomogeneity, 
#the null hypothesis should be an inhomogeneous Poisson process. 
#So we estimate the inhomogeneous intensity functions of the two types of points,
# and then generate simulated point patterns according to these intensities.

#First fit a Poisson point process model to the observed data
#fit2 <- ppm(P, ~marks * polynom(x,y,3))

#n <- 199 #number of iterations for Monte carlo simulation
#p  <- 0.01 # Desired p significance level to display

#Lij_Env_Plot2 <- envelope(fit2, Lcross.inhom, lambdaX=fit2, nsim=n, rank=(p * (n + 1)))
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Env_Plot2, . - r ~ r, ylab=expression(hat("L")),xlab="Distance (m)", main=paste("Lij(r) Plot2 Co3-AsT,  p= ", p),legendmath=TRUE)
#par(OP)
#write.csv(Lij_Env_Plot2, file = "Lij_199MCEnv_Plot2_Co3_to_AsT.csv", row.names=TRUE)



# COMPUTING THE g (Pair Correlation Fucntion) OF COCOA AND ASSOCIATED TREES
#The function value (fv) table method:we will use computed inhomogenous cross K function

# fv object
#Kij_Plot2 is the fv object computed above (inhomogenous cross K from Co to AsT)
#method "b" apply smoothing to Y(r) = K(r)/(2 * pi * r) constraining Y(0) = 0, 
#estimate the derivative of Y, and solve;
#method "b" used due to likely inhibition at small distances of r spar=0.8,

gij_Plot2 <- pcf(Kij_Plot2, spar=0.8,  method="b") 
plot(gij_Plot2, main="Plot2 Co3 - AsT g(r) Pair correlation function", xlab="Distance(m)", legendpos="bottom")
write.csv(gij_Plot2, file = "pcf_Plot2_Co3_to_AsT.csv", row.names=TRUE)


#------------------------------------------------

#PLOT 13


Plot13_Data = read.csv("Plot13_tclass4.csv")
Plot13_Data$TTP = as.factor(Plot13_Data$TTP)

Plot13_T4  <- readShapePoints("Plot13_tclass4.shp")
SP <- as(Plot13_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot13_Data[5]

# Get Plot131 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot13.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot13 windown (polygon) to trees in Plot13
P$window <- W
plot(P, pch=16, main="Plot13 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot13 All Trees")
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
L_AsT_Plot13 <- Ldot.inhom(P, "AsT", lambdaX=fit, update=FALSE)
plot(Lij_Plot13, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot13 Associated Trees", legendpos="bottom")
L_AsT_Plot13$iso  <- L_AsT_Plot13$iso  - L_AsT_Plot13$r
L_AsT_Plot13$theo <- L_AsT_Plot13$theo - L_AsT_Plot13$r
write.csv(L_AsT_Plot13, file = "Li_AsT_Plot13.csv", row.names=TRUE)

#Li (nonhomogenous) for Cocoa Trees2 in stand 

L_Co2_Plot13 <- Ldot.inhom(P, "Co", lambdaX=fit, update=FALSE) #using fit computed for all tree

plot(Lij_Plot13, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot13 Cocoa Trees3", legendpos="bottom")
L_AsT_Plot13$iso  <- L_AsT_Plot13$iso  - L_AsT_Plot13$r
L_AsT_Plot13$theo <- L_AsT_Plot13$theo - L_AsT_Plot13$r
write.csv(L_AsT_Plot13, file = "Li_Co3_Plot13.csv", row.names=TRUE)


#Bivariate(cross) inhomogenous Kij(r) Ripley function
#The Kcross(r) test, which examines distribution a tree group relative to another group within stand

# method (3): fit parametric intensity model
fit <- ppm(P, ~marks * polynom(x,y,2))
# evaluate fitted intensities at data points
# (these are the intensities of the sub-processes of each type)
inten <- fitted(fit, dataonly=TRUE)
# split according to types of points
lambda2 <- split(inten, P$marks)
Kij_Plot13 <- Kcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))
plot(Kij_Plot13, xlab="Distance (m)", main="Kij (2nd order moment) Co3 to Ast")

write.csv(Kij_Plot13, file = "Kij_bivariate_inhom_Plot13_Co3-AsT.csv", row.names=TRUE)



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
Lij_Plot13 <- Lcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))

plot(Lij_Plot13, . - r ~ r, xlab="Distance (m)", ylab="Lij(r)", main="Lij inhom. Function, Plot13 Cocoa3-AsT", legendpos="bottom")


# We will need to modify the plot graphics parameters to accommodate the
# custom y-label; hence the call to `par`.
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Plot13, . -r ~ r, ylab=expression(hat("Lij")), xlab = "d (m)", main="Lij(r), Plot13 Cocoa Trees", legendmath=TRUE, legendpos="bottom")
#par(OP)

#We can make permanent changes to the Lij object's values as follows:
Lij_Plot13$iso  <- Lij_Plot13$iso  - Lij_Plot13$r
Lij_Plot13$theo <- Lij_Plot13$theo - Lij_Plot13$r

write.csv(Lij_Plot13, file = "Lij_Bivariate_inhom_Plot13_Co3_to_AsT.csv", row.names=TRUE)

#MONTE CARLO (MC) SIMULATION

#To test for clustering/inhibition in the presence of spatial inhomogeneity, 
#the null hypothesis should be an inhomogeneous Poisson process. 
#So we estimate the inhomogeneous intensity functions of the two types of points,
# and then generate simulated point patterns according to these intensities.

#First fit a Poisson point process model to the observed data
#fit2 <- ppm(P, ~marks * polynom(x,y,3))

#n <- 199 #number of iterations for Monte carlo simulation
#p  <- 0.01 # Desired p significance level to display

#Lij_Env_Plot13 <- envelope(fit2, Lcross.inhom, lambdaX=fit2, nsim=n, rank=(p * (n + 1)))
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Env_Plot13, . - r ~ r, ylab=expression(hat("L")),xlab="Distance (m)", main=paste("Lij(r) Plot13 Co3-AsT,  p= ", p),legendmath=TRUE)
#par(OP)
#write.csv(Lij_Env_Plot13, file = "Lij_199MCEnv_Plot13_Co3_to_AsT.csv", row.names=TRUE)



# COMPUTING THE g (Pair Correlation Fucntion) OF COCOA AND ASSOCIATED TREES
#The function value (fv) table method:we will use computed inhomogenous cross K function

# fv object
#Kij_Plot13 is the fv object computed above (inhomogenous cross K from Co to AsT)
#method "b" apply smoothing to Y(r) = K(r)/(2 * pi * r) constraining Y(0) = 0, 
#estimate the derivative of Y, and solve;
#method "b" used due to likely inhibition at small distances of r spar=0.8,

gij_Plot13 <- pcf(Kij_Plot13, spar=0.8,  method="b") 
plot(gij_Plot13, main="Plot13 Co3 - AsT g(r) Pair correlation function", xlab="Distance(m)", legendpos="bottom")
write.csv(gij_Plot13, file = "pcf_Plot13_Co3_to_AsT.csv", row.names=TRUE)


#--------------------------------------------------

#PLOT 39


Plot39_Data = read.csv("Plot39_tclass4.csv")
Plot39_Data$TTP = as.factor(Plot39_Data$TTP)

Plot39_T4  <- readShapePoints("Plot39_tclass4.shp")
SP <- as(Plot39_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot39_Data[5]

# Get Plot391 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot39.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot39 windown (polygon) to trees in Plot39
P$window <- W
plot(P, pch=16, main="Plot39 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot39 All Trees")
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
L_AsT_Plot39 <- Ldot.inhom(P, "AsT", lambdaX=fit, update=FALSE)
plot(Lij_Plot39, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot39 Associated Trees", legendpos="bottom")
L_AsT_Plot39$iso  <- L_AsT_Plot39$iso  - L_AsT_Plot39$r
L_AsT_Plot39$theo <- L_AsT_Plot39$theo - L_AsT_Plot39$r
write.csv(L_AsT_Plot39, file = "Li_AsT_Plot39.csv", row.names=TRUE)

#Li (nonhomogenous) for Cocoa Trees2 in stand 

L_Co2_Plot39 <- Ldot.inhom(P, "Co", lambdaX=fit, update=FALSE) #using fit computed for all tree

plot(Lij_Plot39, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot39 Cocoa Trees3", legendpos="bottom")
L_AsT_Plot39$iso  <- L_AsT_Plot39$iso  - L_AsT_Plot39$r
L_AsT_Plot39$theo <- L_AsT_Plot39$theo - L_AsT_Plot39$r
write.csv(L_AsT_Plot39, file = "Li_Co3_Plot39.csv", row.names=TRUE)


#Bivariate(cross) inhomogenous Kij(r) Ripley function
#The Kcross(r) test, which examines distribution a tree group relative to another group within stand

# method (3): fit parametric intensity model
fit <- ppm(P, ~marks * polynom(x,y,2))
# evaluate fitted intensities at data points
# (these are the intensities of the sub-processes of each type)
inten <- fitted(fit, dataonly=TRUE)
# split according to types of points
lambda2 <- split(inten, P$marks)
Kij_Plot39 <- Kcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))
plot(Kij_Plot39, xlab="Distance (m)", main="Kij (2nd order moment) Co3 to Ast")

write.csv(Kij_Plot39, file = "Kij_bivariate_inhom_Plot39_Co3-AsT.csv", row.names=TRUE)



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
Lij_Plot39 <- Lcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))

plot(Lij_Plot39, . - r ~ r, xlab="Distance (m)", ylab="Lij(r)", main="Lij inhom. Function, Plot39 Cocoa3-AsT", legendpos="bottom")


# We will need to modify the plot graphics parameters to accommodate the
# custom y-label; hence the call to `par`.
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Plot39, . -r ~ r, ylab=expression(hat("Lij")), xlab = "d (m)", main="Lij(r), Plot39 Cocoa Trees", legendmath=TRUE, legendpos="bottom")
#par(OP)

#We can make permanent changes to the Lij object's values as follows:
Lij_Plot39$iso  <- Lij_Plot39$iso  - Lij_Plot39$r
Lij_Plot39$theo <- Lij_Plot39$theo - Lij_Plot39$r

write.csv(Lij_Plot39, file = "Lij_Bivariate_inhom_Plot39_Co3_to_AsT.csv", row.names=TRUE)

#MONTE CARLO (MC) SIMULATION

#To test for clustering/inhibition in the presence of spatial inhomogeneity, 
#the null hypothesis should be an inhomogeneous Poisson process. 
#So we estimate the inhomogeneous intensity functions of the two types of points,
# and then generate simulated point patterns according to these intensities.

#First fit a Poisson point process model to the observed data
#fit2 <- ppm(P, ~marks * polynom(x,y,3))

#n <- 199 #number of iterations for Monte carlo simulation
#p  <- 0.01 # Desired p significance level to display

#Lij_Env_Plot39 <- envelope(fit2, Lcross.inhom, lambdaX=fit2, nsim=n, rank=(p * (n + 1)))
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Env_Plot39, . - r ~ r, ylab=expression(hat("L")),xlab="Distance (m)", main=paste("Lij(r) Plot39 Co3-AsT,  p= ", p),legendmath=TRUE)
#par(OP)
#write.csv(Lij_Env_Plot39, file = "Lij_199MCEnv_Plot39_Co3_to_AsT.csv", row.names=TRUE)



# COMPUTING THE g (Pair Correlation Fucntion) OF COCOA AND ASSOCIATED TREES
#The function value (fv) table method:we will use computed inhomogenous cross K function

# fv object
#Kij_Plot39 is the fv object computed above (inhomogenous cross K from Co to AsT)
#method "b" apply smoothing to Y(r) = K(r)/(2 * pi * r) constraining Y(0) = 0, 
#estimate the derivative of Y, and solve;
#method "b" used due to likely inhibition at small distances of r spar=0.8,

gij_Plot39 <- pcf(Kij_Plot39, spar=0.8,  method="b") 
plot(gij_Plot39, main="Plot39 Co3 - AsT g(r) Pair correlation function", xlab="Distance(m)", legendpos="bottom")
write.csv(gij_Plot39, file = "pcf_Plot39_Co3_to_AsT.csv", row.names=TRUE)


#--------------------------------------------------

#PLOT 44


Plot44_Data = read.csv("Plot44_tclass4.csv")
Plot44_Data$TTP = as.factor(Plot44_Data$TTP)

Plot44_T4  <- readShapePoints("Plot44_tclass4.shp")
SP <- as(Plot44_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot44_Data[5]

# Get Plot441 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot44.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot44 windown (polygon) to trees in Plot44
P$window <- W
plot(P, pch=16, main="Plot44 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot44 All Trees")
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
L_AsT_Plot44 <- Ldot.inhom(P, "AsT", lambdaX=fit, update=FALSE)
plot(Lij_Plot44, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot44 Associated Trees", legendpos="bottom")
L_AsT_Plot44$iso  <- L_AsT_Plot44$iso  - L_AsT_Plot44$r
L_AsT_Plot44$theo <- L_AsT_Plot44$theo - L_AsT_Plot44$r
write.csv(L_AsT_Plot44, file = "Li_AsT_Plot44.csv", row.names=TRUE)

#Li (nonhomogenous) for Cocoa Trees2 in stand 

L_Co2_Plot44 <- Ldot.inhom(P, "Co", lambdaX=fit, update=FALSE) #using fit computed for all tree

plot(Lij_Plot44, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot44 Cocoa Trees3", legendpos="bottom")
L_AsT_Plot44$iso  <- L_AsT_Plot44$iso  - L_AsT_Plot44$r
L_AsT_Plot44$theo <- L_AsT_Plot44$theo - L_AsT_Plot44$r
write.csv(L_AsT_Plot44, file = "Li_Co3_Plot44.csv", row.names=TRUE)


#Bivariate(cross) inhomogenous Kij(r) Ripley function
#The Kcross(r) test, which examines distribution a tree group relative to another group within stand

# method (3): fit parametric intensity model
fit <- ppm(P, ~marks * polynom(x,y,2))
# evaluate fitted intensities at data points
# (these are the intensities of the sub-processes of each type)
inten <- fitted(fit, dataonly=TRUE)
# split according to types of points
lambda2 <- split(inten, P$marks)
Kij_Plot44 <- Kcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))
plot(Kij_Plot44, xlab="Distance (m)", main="Kij (2nd order moment) Co3 to Ast")

write.csv(Kij_Plot44, file = "Kij_bivariate_inhom_Plot44_Co3-AsT.csv", row.names=TRUE)



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
Lij_Plot44 <- Lcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))

plot(Lij_Plot44, . - r ~ r, xlab="Distance (m)", ylab="Lij(r)", main="Lij inhom. Function, Plot44 Cocoa3-AsT", legendpos="bottom")


# We will need to modify the plot graphics parameters to accommodate the
# custom y-label; hence the call to `par`.
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Plot44, . -r ~ r, ylab=expression(hat("Lij")), xlab = "d (m)", main="Lij(r), Plot44 Cocoa Trees", legendmath=TRUE, legendpos="bottom")
#par(OP)

#We can make permanent changes to the Lij object's values as follows:
Lij_Plot44$iso  <- Lij_Plot44$iso  - Lij_Plot44$r
Lij_Plot44$theo <- Lij_Plot44$theo - Lij_Plot44$r

write.csv(Lij_Plot44, file = "Lij_Bivariate_inhom_Plot44_Co3_to_AsT.csv", row.names=TRUE)

#MONTE CARLO (MC) SIMULATION

#To test for clustering/inhibition in the presence of spatial inhomogeneity, 
#the null hypothesis should be an inhomogeneous Poisson process. 
#So we estimate the inhomogeneous intensity functions of the two types of points,
# and then generate simulated point patterns according to these intensities.

#First fit a Poisson point process model to the observed data
#fit2 <- ppm(P, ~marks * polynom(x,y,3))

#n <- 199 #number of iterations for Monte carlo simulation
#p  <- 0.01 # Desired p significance level to display

#Lij_Env_Plot44 <- envelope(fit2, Lcross.inhom, lambdaX=fit2, nsim=n, rank=(p * (n + 1)))
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Env_Plot44, . - r ~ r, ylab=expression(hat("L")),xlab="Distance (m)", main=paste("Lij(r) Plot44 Co3-AsT,  p= ", p),legendmath=TRUE)
#par(OP)
#write.csv(Lij_Env_Plot44, file = "Lij_199MCEnv_Plot44_Co3_to_AsT.csv", row.names=TRUE)



# COMPUTING THE g (Pair Correlation Fucntion) OF COCOA AND ASSOCIATED TREES
#The function value (fv) table method:we will use computed inhomogenous cross K function

# fv object
#Kij_Plot44 is the fv object computed above (inhomogenous cross K from Co to AsT)
#method "b" apply smoothing to Y(r) = K(r)/(2 * pi * r) constraining Y(0) = 0, 
#estimate the derivative of Y, and solve;
#method "b" used due to likely inhibition at small distances of r spar=0.8,

gij_Plot44 <- pcf(Kij_Plot44, spar=0.8,  method="b") 
plot(gij_Plot44, main="Plot44 Co3 - AsT g(r) Pair correlation function", xlab="Distance(m)", legendpos="bottom")
write.csv(gij_Plot44, file = "pcf_Plot44_Co3_to_AsT.csv", row.names=TRUE)


#--------------------------------------------------

#PLOT 50


Plot50_Data = read.csv("Plot50_tclass4.csv")
Plot50_Data$TTP = as.factor(Plot50_Data$TTP)

Plot50_T4  <- readShapePoints("Plot50_tclass4.shp")
SP <- as(Plot50_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot50_Data[5]

# Get Plot501 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot50.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot50 windown (polygon) to trees in Plot50
P$window <- W
plot(P, pch=16, main="Plot50 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot50 All Trees")
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
L_AsT_Plot50 <- Ldot.inhom(P, "AsT", lambdaX=fit, update=FALSE)
plot(Lij_Plot50, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot50 Associated Trees", legendpos="bottom")
L_AsT_Plot50$iso  <- L_AsT_Plot50$iso  - L_AsT_Plot50$r
L_AsT_Plot50$theo <- L_AsT_Plot50$theo - L_AsT_Plot50$r
write.csv(L_AsT_Plot50, file = "Li_AsT_Plot50.csv", row.names=TRUE)

#Li (nonhomogenous) for Cocoa Trees2 in stand 

L_Co2_Plot50 <- Ldot.inhom(P, "Co", lambdaX=fit, update=FALSE) #using fit computed for all tree

plot(Lij_Plot50, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot50 Cocoa Trees3", legendpos="bottom")
L_AsT_Plot50$iso  <- L_AsT_Plot50$iso  - L_AsT_Plot50$r
L_AsT_Plot50$theo <- L_AsT_Plot50$theo - L_AsT_Plot50$r
write.csv(L_AsT_Plot50, file = "Li_Co3_Plot50.csv", row.names=TRUE)


#Bivariate(cross) inhomogenous Kij(r) Ripley function
#The Kcross(r) test, which examines distribution a tree group relative to another group within stand

# method (3): fit parametric intensity model
fit <- ppm(P, ~marks * polynom(x,y,2))
# evaluate fitted intensities at data points
# (these are the intensities of the sub-processes of each type)
inten <- fitted(fit, dataonly=TRUE)
# split according to types of points
lambda2 <- split(inten, P$marks)
Kij_Plot50 <- Kcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))
plot(Kij_Plot50, xlab="Distance (m)", main="Kij (2nd order moment) Co3 to Ast")

write.csv(Kij_Plot50, file = "Kij_bivariate_inhom_Plot50_Co3-AsT.csv", row.names=TRUE)



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
Lij_Plot50 <- Lcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))

plot(Lij_Plot50, . - r ~ r, xlab="Distance (m)", ylab="Lij(r)", main="Lij inhom. Function, Plot50 Cocoa3-AsT", legendpos="bottom")


# We will need to modify the plot graphics parameters to accommodate the
# custom y-label; hence the call to `par`.
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Plot50, . -r ~ r, ylab=expression(hat("Lij")), xlab = "d (m)", main="Lij(r), Plot50 Cocoa Trees", legendmath=TRUE, legendpos="bottom")
#par(OP)

#We can make permanent changes to the Lij object's values as follows:
Lij_Plot50$iso  <- Lij_Plot50$iso  - Lij_Plot50$r
Lij_Plot50$theo <- Lij_Plot50$theo - Lij_Plot50$r

write.csv(Lij_Plot50, file = "Lij_Bivariate_inhom_Plot50_Co3_to_AsT.csv", row.names=TRUE)

#MONTE CARLO (MC) SIMULATION

#To test for clustering/inhibition in the presence of spatial inhomogeneity, 
#the null hypothesis should be an inhomogeneous Poisson process. 
#So we estimate the inhomogeneous intensity functions of the two types of points,
# and then generate simulated point patterns according to these intensities.

#First fit a Poisson point process model to the observed data
#fit2 <- ppm(P, ~marks * polynom(x,y,3))

#n <- 199 #number of iterations for Monte carlo simulation
#p  <- 0.01 # Desired p significance level to display

#Lij_Env_Plot50 <- envelope(fit2, Lcross.inhom, lambdaX=fit2, nsim=n, rank=(p * (n + 1)))
#OP <- par(mar=c(5,5,4,4))
#plot(Lij_Env_Plot50, . - r ~ r, ylab=expression(hat("L")),xlab="Distance (m)", main=paste("Lij(r) Plot50 Co3-AsT,  p= ", p),legendmath=TRUE)
#par(OP)
#write.csv(Lij_Env_Plot50, file = "Lij_199MCEnv_Plot50_Co3_to_AsT.csv", row.names=TRUE)



# COMPUTING THE g (Pair Correlation Fucntion) OF COCOA AND ASSOCIATED TREES
#The function value (fv) table method:we will use computed inhomogenous cross K function

# fv object
#Kij_Plot50 is the fv object computed above (inhomogenous cross K from Co to AsT)
#method "b" apply smoothing to Y(r) = K(r)/(2 * pi * r) constraining Y(0) = 0, 
#estimate the derivative of Y, and solve;
#method "b" used due to likely inhibition at small distances of r spar=0.8,

gij_Plot50 <- pcf(Kij_Plot50, spar=0.8,  method="b") 
plot(gij_Plot50, main="Plot50 Co3 - AsT g(r) Pair correlation function", xlab="Distance(m)", legendpos="bottom")
write.csv(gij_Plot50, file = "pcf_Plot50_Co3_to_AsT.csv", row.names=TRUE)


# SAVING AND APPENDING FILES IN CSV


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

# CREATIONG OF HYPERFRAME FOR COAF1 AND ASSOCIATED TREES

library(sp)
library(spatstat)
library(maptools)
library(ggplot2)

setwd("D:/PhD Manuscript/...../Bakoa_TreesClass4")


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


#----------------------------------------------------------------------------

#PLOT 23B

Plot23b_Data = read.csv("Plot23B_tclass4.csv")
Plot23b_Data$TTP = as.factor(Plot23b_Data$TTP)

Plot23b_T4  <- readShapePoints("Plot23B_tclass4.shp")
SP <- as(Plot23b_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot23b_Data[33]

# Get Plot23b shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot23B.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot23b windown (polygon) to trees in Plot23b
P$window <- W
a=plot(P, pch=16, main="Plot23 Cocoa Trees")
#dev.copy(tiff,'Plot23b Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot23b")

P23bT4 = P
summary(P23bT4)

P23bT4Co = unmark(split(P23bT4)$Co) # To subset the cocoa trees 
P23bTAsT = unmark(split(P23bT4)$AsT) # To subset Associated trees
P23bT4All = unmark(P23bT4) # Remove marks (tree type) from point pattern data

#--------------------------------------------------------------------

#PLOT 23C

Plot23c_Data = read.csv("Plot23C_tclass4.csv")
Plot23c_Data$TTP = as.factor(Plot23c_Data$TTP)


Plot23C_T4  <- readShapePoints("Plot23C_tclass4.shp")
SP <- as(Plot23C_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot23c_Data[33]

# Get Plot23c shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot23C.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot23C windown (polygon) to trees in Plot23c
P$window <- W
a=plot(P, pch=16, main="Plot23c Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot23c")
# Get summary of Cocoa trees in plot 26c
P23cT4 = P
summary(P23cT4)


P23cT4Co = unmark(split(P23cT4)$Co) # To subset just the cocoa trees 
P23cTAsT = unmark(split(P23cT4)$AsT) # To subset jut Associated trees
P23cT4All = unmark(P23cT4) # Remove marks (tree type) from point pattern data



#----------------------------

#PLOT 26

Plot26_Data = read.csv("Plot26_tclass4.csv")
Plot26_Data$TTP = as.factor(Plot26_Data$TTP)

Plot26_T4  <- readShapePoints("Plot26_tclass4.shp")
SP <- as(Plot26_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot26_Data[5]

# Get Plot26 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot28.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot26 window (polygon) to trees in Plot26
P$window <- W
plot(P, pch=16, main="Plot26 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot26 All Trees")
plot(density(P))

P26T4 = P
summary(P26T4)
P26T4Co = unmark(split(P26T4)$Co) # To subset just the cocoa trees 
P26TAsT = unmark(split(P26T4)$AsT) # To subset jut Associated trees
P26T4All = unmark(P26T4) # Remove marks (tree type) from point pattern data

#--------------------------------------------------------------------

#PLOT 26Aa

Plot26Aa_Data = read.csv("Plot26Aa_tclass4.csv")
Plot26Aa_Data$TTP = as.factor(Plot26Aa_Data$TTP)

Plot26Aa_T4  <- readShapePoints("Plot26Aa_tclass4.shp")
SP <- as(Plot26Aa_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot26Aa_Data[33]

# Get Plot26Aa shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot26Aa.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot26Aa windown (polygon) to trees in Plot26Aa
P$window <- W
a=plot(P, pch=16, main="Plot26Aa Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot26Aa")
# Get summary of Cocoa trees in plot 26Aa
P26AaT4 = P
summary(P26AaT4)
P26AaT4Co = unmark(split(P26AaT4)$Co) # To subset just the cocoa trees 
P26AaTAsT = unmark(split(P26AaT4)$AsT) # To subset jut Associated trees
P26AaT4All = unmark(P26AaT4) # Remove marks (tree type) from point pattern data

#--------------------------------------------------------------------

#PLOT 26b

Plot26b_Data = read.csv("Plot26B_tclass4.csv")
Plot26b_Data$TTP = as.factor(Plot26b_Data$TTP)

Plot26b_T4  <- readShapePoints("Plot26B_tclass4.shp")
SP <- as(Plot26b_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot26b_Data[33]

# Get Plot26b shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot26b.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot26b windown (polygon) to trees in Plot26b
P$window <- W
a=plot(P, pch=16, main="Plot26b Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot26b")
# Get summary of Cocoa trees in plot 26b
P26bT4 = P
summary(P26bT4)

P26bT4Co = unmark(split(P26bT4)$Co) # To subset just the cocoa trees 
P26bTAsT = unmark(split(P26bT4)$AsT) # To subset jut Associated trees
P26bT4All = unmark(P26bT4) # Remove marks (tree type) from point pattern data


#--------------------------------------------------------------------

#PLOT 26c

Plot26c_Data = read.csv("Plot26C_tclass4.csv")
Plot26c_Data$TTP = as.factor(Plot26c_Data$TTP)


Plot26c_T4  <- readShapePoints("Plot26C_tclass4.shp")
SP <- as(Plot26c_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot26c_Data[33]

# Get Plot26c shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot26C.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot26c windown (polygon) to trees in Plot26c
P$window <- W
a=plot(P, pch=16, main="Plot26c Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot26c")
# Get summary of Cocoa trees in plot 26c
P26cT4 = P
summary(P26cT4)

P26cT4Co = unmark(split(P26cT4)$Co) # To subset just the cocoa trees 
P26cTAsT = unmark(split(P26cT4)$AsT) # To subset jut Associated trees
P26cT4All = unmark(P26cT4) # Remove marks (tree type) from point pattern data


#---------------------------------
#PLOT 28

Plot28_Data = read.csv("Plot28_tclass4.csv")
Plot28_Data$TTP = as.factor(Plot28_Data$TTP)

Plot28_T4  <- readShapePoints("Plot28_tclass4.shp")
SP <- as(Plot28_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot28_Data[5]

# Get Plot28 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot26.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot28 window (polygon) to trees in Plot28
P$window <- W
plot(P, pch=16, main="Plot28 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot28 All Trees")
plot(density(P))

P28T4 = P
summary(P28T4)
P28T4Co = unmark(split(P28T4)$Co) # To subset the cocoa trees 
P28TAsT = unmark(split(P28T4)$AsT) # To subset Associated trees
P28T4All = unmark(P28T4) # Remove marks (tree type) from point pattern data


#--------------------------------------------------------------------

#PLOT 28Aa

Plot28Aa_Data = read.csv("Plot28Aa_tclass4.csv")
Plot28Aa_Data$TTP = as.factor(Plot28Aa_Data$TTP)

Plot28Aa_T4  <- readShapePoints("Plot28Aa_tclass4.shp")
SP <- as(Plot28Aa_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot28Aa_Data[33] # set mark to point pattern using TTP column

# Get Plot28Aa shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot28Aa.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot28Aa windown (polygon) to trees in Plot28Aa
P$window <- W
summary(P)
plot(split(P))
plot(density(split(P)))
a=plot(P, pch=16, main="Plot28Aa All Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot28Aa")
# Get summary of Cocoa trees in plot 28Aa
P28AaT4 = P
summary(P28AaT4)

P28AaT4Co = unmark(split(P28AaT4)$Co) # To subset the cocoa trees 
P28AaTAsT = unmark(split(P28AaT4)$AsT) # To subset Associated trees
P28AaT4All = unmark(P28AaT4) # Remove marks (tree type) from point pattern data

#--------------------------------------------------------------------

#PLOT 28b

Plot28b_Data = read.csv("Plot28B_tclass4.csv")
Plot28b_Data$TTP = as.factor(Plot28b_Data$TTP)

Plot28b_T4  <- readShapePoints("Plot28B_tclass4.shp")
SP <- as(Plot28b_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot28b_Data[33] # set mark to point pattern using TTP column

# Get Plot28b shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot28b.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot28b windown (polygon) to trees in Plot28b
P$window <- W
a=plot(P, pch=16, main="Plot28b Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot28b")
# Get summary of Cocoa trees in plot 28b
P28bT4 = P
summary(P28bT4)

P28bT4Co = unmark(split(P28bT4)$Co) # To subset the cocoa trees 
P28bTAsT = unmark(split(P28bT4)$AsT) # To subset Associated trees
P28bT4All = unmark(P28bT4) # Remove marks (tree type) from point pattern data

#--------------------------------------------------------------------

#PLOT 28c

Plot28c_Data = read.csv("Plot28C_tclass4.csv")
Plot28c_Data$TTP = as.factor(Plot28c_Data$TTP)

Plot28c_T4  <- readShapePoints("Plot28C_tclass4.shp")
SP <- as(Plot28c_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot28c_Data[33] # set marks to point pattern using TTP Column

# Get Plot28c shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot28C.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot28c windown (polygon) to trees in Plot28c
P$window <- W
a=plot(P, pch=16, main="Plot28c Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot28c")
# Get summary of Cocoa trees in plot 28c
P28cT4 = P
summary(P28cT4)

P28cT4Co = unmark(split(P28cT4)$Co) # To subset the cocoa trees 
P28cTAsT = unmark(split(P28cT4)$AsT) # To subset Associated trees
P28cT4All = unmark(P28cT4) # Remove marks (tree type) from point pattern data


#----------------------------

#PLOT 38


Plot38_Data = read.csv("Plot38_tclass4.csv")
Plot38_Data$TTP = as.factor(Plot38_Data$TTP)

Plot38_T4  <- readShapePoints("Plot38_tclass4.shp")
SP <- as(Plot38_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot38_Data[5]

# Get Plot38 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot38.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot38 window (polygon) to trees in Plot38
P$window <- W
plot(P, pch=16, main="Plot38 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot38 All Trees")
plot(density(P))

P38T4 = P
summary(P38T4)
P38T4Co = unmark(split(P38T4)$Co) # To subset just the cocoa trees 
P38TAsT = unmark(split(P38T4)$AsT) # To subset jut Associated trees
P38T4All = unmark(P38T4) # Remove marks (tree type) from point pattern data


#---------------------
#PLOT 4


Plot4_Data = read.csv("Plot4_tclass4.csv")
Plot4_Data$TTP = as.factor(Plot4_Data$TTP)

Plot4_T4  <- readShapePoints("Plot4_tclass4.shp")
SP <- as(Plot4_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot4_Data[5]

# Get Plot4 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot4.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot4 window (polygon) to trees in Plot4
P$window <- W
plot(P, pch=16, main="Plot4 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot4 All Trees")
plot(density(P))
#legend(0.25, 0.5, names(a), pch=a)

P4T4 = P
summary(P4T4)
P4T4Co = unmark(split(P4T4)$Co) # To subset just the cocoa trees 
P4TAsT = unmark(split(P4T4)$AsT) # To subset jut Associated trees
P4T4All = unmark(P4T4) # Remove marks (tree type) from point pattern data



#-------------------------

#PLOT 22


Plot22_Data = read.csv("Plot22_tclass4.csv")
Plot22_Data$TTP = as.factor(Plot22_Data$TTP)

Plot22_T4  <- readShapePoints("Plot22_tclass4.shp")
SP <- as(Plot22_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot22_Data[5]

# Get Plot22 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot22.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot22 window (polygon) to trees in Plot22
P$window <- W
plot(P, pch=16, main="Plot22 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot22 All Trees")
plot(density(P))

P22T4 = P
summary(P22T4)
P22T4Co = unmark(split(P22T4)$Co) # To subset the cocoa trees 
P22TAsT = unmark(split(P22T4)$AsT) # To subset Associated trees
P22T4All = unmark(P22T4) # Remove marks (tree type) from point pattern data


#--------------------------------------------------------------------


#PLOT 22Aa

Plot22Aa_Data = read.csv("Plot22Aa_tclass4.csv")
Plot22Aa_Data$TTP = as.factor(Plot22Aa_Data$TTP)

Plot22Aa_T4  <- readShapePoints("Plot22Aa_tclass4.shp")
SP <- as(Plot22Aa_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot22Aa_Data[33]

# Get Plot22Aa shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot22Aa.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot22Aa windown (polygon) to trees in Plot22Aa
P$window <- W
a=plot(P, pch=16, main="Plot22Aa Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot22Aa")
# Get summary of Cocoa trees in plot 22Aa
P22AaT4 = P
summary(P22AaT4)

P22AaT4Co = unmark(split(P22AaT4)$Co) # To subset the cocoa trees 
P22AaTAsT = unmark(split(P22AaT4)$AsT) # To subset Associated trees
P22AaT4All = unmark(P22AaT4) # Remove marks (tree type) from point pattern data

#--------------------------------------------------------------------

#PLOT 22b

Plot22b_Data = read.csv("Plot22B_tclass4.csv")
Plot22b_Data$TTP = as.factor(Plot22b_Data$TTP)

Plot22b_T4  <- readShapePoints("Plot22B_tclass4.shp")
SP <- as(Plot22b_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot22b_Data[33]

# Get Plot2b shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot22b.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot22b windown (polygon) to trees in Plot22b
P$window <- W
a=plot(P, pch=16, main="Plot22b Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot22b")
# Get summary of Cocoa trees in plot 22b
P22bT4 = P
summary(P22bT4)

P22bT4Co = unmark(split(P22bT4)$Co) # To subset the cocoa trees 
P22bTAsT = unmark(split(P22bT4)$AsT) # To subset Associated trees
P22bT4All = unmark(P22bT4) # Remove marks (tree type) from point pattern data

#--------------------------------------------------------------------

#PLOT 22c

Plot22c_Data = read.csv("Plot22C_tclass4.csv")
Plot22c_Data$TTP = as.factor(Plot22c_Data$TTP)

Plot22c_T4  <- readShapePoints("Plot22C_tclass4.shp")
SP <- as(Plot22c_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot22c_Data[33]

# Get Plot22c shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot22C.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot22c windown (polygon) to trees in Plot22c
P$window <- W
a=plot(P, pch=16, main="Plot22c Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot22c")
# Get summary of Cocoa trees in plot 22c
P22cT4 = P
summary(P22cT4)

P22cT4Co = unmark(split(P22cT4)$Co) # To subset the cocoa trees 
P22cTAsT = unmark(split(P22cT4)$AsT) # To subset Associated trees
P22cT4All = unmark(P22cT4) # Remove marks (tree type) from point pattern data


#------------------------
#PLOT 6


Plot6_Data = read.csv("Plot6_tclass4.csv")
Plot6_Data$TTP = as.factor(Plot6_Data$TTP)

Plot6_T4  <- readShapePoints("Plot6_tclass4.shp")
SP <- as(Plot6_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot6_Data[5]

# Get Plot6 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot6.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot6 window (polygon) to trees in Plot6
P$window <- W
plot(P, pch=16, main="Plot6 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot6 All Trees")
plot(density(P))

P6T4 = P
summary(P6T4)
P6T4Co = unmark(split(P6T4)$Co) # To subset the cocoa trees 
P6TAsT = unmark(split(P6T4)$AsT) # To subset Associated trees
P6T4All = unmark(P6T4) # Remove marks (tree type) from point pattern data

#--------------------------------------------------------------------

#PLOT 6b

Plot6b_Data = read.csv("Plot6B_tclass4.csv")
Plot6b_Data$TTP = as.factor(Plot6b_Data$TTP)

Plot6b_T4  <- readShapePoints("Plot6B_tclass4.shp")
SP <- as(Plot6b_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot6b_Data[33]

# Get Plot6b shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot6b.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot6b windown (polygon) to trees in Plot6b
P$window <- W
a=plot(P, pch=16, main="Plot6b Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot6b")
# Get summary of Cocoa trees in plot 6b
P6bT4 = P
summary(P6bT4)

P6bT4Co = unmark(split(P6bT4)$Co) # To subset the cocoa trees 
P6bTAsT = unmark(split(P6bT4)$AsT) # To subset Associated trees
P6bT4All = unmark(P6bT4) # Remove marks (tree type) from point pattern data

#-------------------------
#PLOT 43


Plot43_Data = read.csv("Plot43_tclass4.csv")
Plot43_Data$TTP = as.factor(Plot43_Data$TTP)

Plot43_T4  <- readShapePoints("Plot43_tclass4.shp")
SP <- as(Plot43_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot43_Data[5]

# Get Plot43 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot43.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot43 window (polygon) to trees in Plot43
P$window <- W
plot(P, pch=16, main="Plot43 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot43 All Trees")
plot(density(P))

P43T4 = P
summary(P43T4)
P43T4Co = unmark(split(P43T4)$Co) # To subset just the cocoa trees 
P43TAsT = unmark(split(P43T4)$AsT) # To subset jut Associated trees
P43T4All = unmark(P43T4) # Remove marks (tree type) from point pattern data



#------------------------
#PLOT 5


Plot5_Data = read.csv("Plot5_tclass4.csv")
Plot5_Data$TTP = as.factor(Plot5_Data$TTP)

Plot5_T4  <- readShapePoints("Plot5_tclass4.shp")
SP <- as(Plot5_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot5_Data[5]

# Get Plot5 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot5.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot5 window (polygon) to trees in Plot5
P$window <- W
plot(P, pch=16, main="Plot5 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot5 All Trees")
plot(density(P))

P5T4 = P
summary(P5T4)
P5T4Co = unmark(split(P5T4)$Co) # To subset the cocoa trees 
P5TAsT = unmark(split(P5T4)$AsT) # To subset Associated trees
P5T4All = unmark(P5T4) # Remove marks (tree type) from point pattern data


#--------------------------------------------------------------------

#PLOT 5b

Plot5b_Data = read.csv("Plot5B_tclass4.csv")
Plot5b_Data$TTP = as.factor(Plot5b_Data$TTP)


Plot5b_T4  <- readShapePoints("Plot5B_tclass4.shp")
SP <- as(Plot5b_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot5b_Data[33]

# Get Plot5b shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot5b.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot5b windown (polygon) to trees in Plot5b
P$window <- W
summary(P)
plot(P)
plot(split(P))
plot(density(split(P)))
a=plot(P, pch=16, main="Plot5b All Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot5b")
# Get summary of Cocoa trees in plot 5b
P5bT4 = P
summary(P5bT4)

P5bT4Co = unmark(split(P5bT4)$Co) # To subset the cocoa trees 
P5bTAsT = unmark(split(P5bT4)$AsT) # To subset Associated trees
P5bT4All = unmark(P5bT4) # Remove marks (tree type) from point pattern data


#---------------------------
#PLOT 21


Plot21_Data = read.csv("Plot21_tclass4.csv")
Plot21_Data$TTP = as.factor(Plot21_Data$TTP)

Plot21_T4  <- readShapePoints("Plot21_tclass4.shp")
SP <- as(Plot21_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot21_Data[5]

# Get Plot21 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot21.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot21 window (polygon) to trees in Plot21
P$window <- W
plot(P, pch=16, main="Plot21 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot21 All Trees")
plot(density(P))

P21T4 = P
summary(P21T4)
P21T4Co = unmark(split(P21T4)$Co) # To subset the cocoa trees 
P21TAsT = unmark(split(P21T4)$AsT) # To subset Associated trees
P21T4All = unmark(P21T4) # Remove marks (tree type) from point pattern data

#--------------------------------------------------------------------

#PLOT 21Aa

Plot21Aa_Data = read.csv("Plot21Aa_tclass4.csv")
Plot21Aa_Data$TTP = as.factor(Plot21Aa_Data$TTP)

Plot21Aa_T4  <- readShapePoints("Plot21Aa_tclass4.shp")
SP <- as(Plot21Aa_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot21Aa_Data[33]

# Get Plot21Aa shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot21Aa.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot21Aa windown (polygon) to trees in Plot21Aa
P$window <- W
a=plot(P, pch=16, main="Plot21Aa Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot21Aa")
# Get summary of Cocoa trees in plot 21Aa
P21AaT4 = P
summary(P21AaT4)

P21AaT4Co = unmark(split(P21AaT4)$Co) # To subset the cocoa trees 
P21AaTAsT = unmark(split(P21AaT4)$AsT) # To subset Associated trees
P21AaT4All = unmark(P21AaT4) # Remove marks (tree type) from point pattern data


#--------------------------------------------------------------------

#PLOT 21b

Plot21b_Data = read.csv("Plot21B_tclass4.csv")
Plot21b_Data$TTP = as.factor(Plot21b_Data$TTP)

Plot21b_T4  <- readShapePoints("Plot21B_tclass4.shp")
SP <- as(Plot21b_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot21b_Data[33]

# Get Plot21b shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot21B.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot21b windown (polygon) to trees in Plot21b
P$window <- W
a=plot(P, pch=16, main="Plot21b Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot21b")
# Get summary of Cocoa trees in plot 21b
P21bT4 = P
summary(P21bT4)

P21bT4Co = unmark(split(P21bT4)$Co) # To subset the cocoa trees 
P21bTAsT = unmark(split(P21bT4)$AsT) # To subset Associated trees
P21bT4All = unmark(P21bT4) # Remove marks (tree type) from point pattern data

#--------------------------------------------------------------------

#PLOT 21c

Plot21c_Data = read.csv("Plot21C_tclass4.csv")
Plot21c_Data$TTP = as.factor(Plot21c_Data$TTP)

Plot21c_T4  <- readShapePoints("Plot21C_tclass4.shp")
SP <- as(Plot21c_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot21c_Data[33]

# Get Plot21c shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot21C.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot21c windown (polygon) to trees in Plot21c
P$window <- W
summary(P)
plot(split(P))
plot(density(split(P)))
a=plot(P, pch=16, main="Plot21c Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot21c")
# Get summary of Cocoa trees in plot 21c
P21cT4 = P
summary(P21cT4)

P21cT4Co = unmark(split(P21cT4)$Co) # To subset the cocoa trees 
P21cTAsT = unmark(split(P21cT4)$AsT) # To subset Associated trees
P21cT4All = unmark(P21cT4) # Remove marks (tree type) from point pattern data


#----------------------------
#PLOT 8


Plot8_Data = read.csv("Plot8_tclass4.csv")
Plot8_Data$TTP = as.factor(Plot8_Data$TTP)

Plot8_T4  <- readShapePoints("Plot8_tclass4.shp")
SP <- as(Plot8_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot8_Data[5]

# Get Plot8 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot8.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot8 window (polygon) to trees in Plot8
P$window <- W
plot(P, pch=16, main="Plot8 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot8 All Trees")
plot(density(P))

P8T4 = P
summary(P8T4)
P8T4Co = unmark(split(P8T4)$Co) # To subsetthe cocoa trees 
P8TAsT = unmark(split(P8T4)$AsT) # To subset Associated trees
P8T4All = unmark(P8T4) # Remove marks (tree type) from point pattern data

#--------------------------------------------------------------------

#PLOT 8b

Plot8b_Data = read.csv("Plot8B_tclass4.csv")
Plot8b_Data$TTP = as.factor(Plot8b_Data$TTP)

Plot8b_T4  <- readShapePoints("Plot8B_tclass4.shp")
SP <- as(Plot8b_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot8b_Data[33]

# Get Plot8b shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot8b.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot8b windown (polygon) to trees in Plot8b
P$window <- W
a=plot(P, pch=16, main="Plot8b Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot8b")
# Get summary of Cocoa trees in plot 8b
P8bT4 = P
summary(P8bT4)

P8bT4Co = unmark(split(P8bT4)$Co) # To subsetthe cocoa trees 
P8bTAsT = unmark(split(P8bT4)$AsT) # To subset Associated trees
P8bT4All = unmark(P8bT4) # Remove marks (tree type) from point pattern data


#-------------------------------
#PLOT 34


Plot34_Data = read.csv("Plot34_tclass4.csv")
Plot34_Data$TTP = as.factor(Plot34_Data$TTP)

Plot34_T4  <- readShapePoints("Plot34_tclass4.shp")
SP <- as(Plot34_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot34_Data[5]

# Get Plot34 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot34.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot34 window (polygon) to trees in Plot34
P$window <- W
plot(P, pch=16, main="Plot34 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot34 All Trees")
plot(density(P))

P34T4 = P
summary(P)
P34T4Co = unmark(split(P34T4)$Co) # To subset the cocoa trees 
P34TAsT = unmark(split(P34T4)$AsT) # To subset Associated trees
P34T4All = unmark(P34T4) # Remove marks (tree type) from point pattern data


#--------------------------------------------------------------------

#PLOT 34b

Plot34b_Data = read.csv("Plot34B_tclass4.csv")
Plot34b_Data$TTP = as.factor(Plot34b_Data$TTP)

Plot34b_T4  <- readShapePoints("Plot34B_tclass4.shp")
SP <- as(Plot34b_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot34b_Data[33]

# Get Plot34b shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot34b.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot34b windown (polygon) to trees in Plot34b
P$window <- W
a=plot(P, pch=16, main="Plot34b Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot34b")
# Get summary of Cocoa trees in plot 34b
P34bT4 = P
summary(P34bT4)

P34bT4Co = unmark(split(P34bT4)$Co) # To subset the cocoa trees 
P34bTAsT = unmark(split(P34bT4)$AsT) # To subset Associated trees
P34bT4All = unmark(P34bT4) # Remove marks (tree type) from point pattern data


#--------------------------------------------------------------------

#PLOT 34c

Plot34c_Data = read.csv("Plot34C_tclass4.csv")
Plot34c_Data$TTP = as.factor(Plot34c_Data$TTP)

Plot34c_T4  <- readShapePoints("Plot34C_tclass4.shp")
SP <- as(Plot34c_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot34c_Data[33]

# Get Plot34c shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot34C.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot34c windown (polygon) to trees in Plot34c
P$window <- W
a=plot(P, pch=16, main="Plot34c Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot34c")
# Get summary of Cocoa trees in plot 34c
P34cT4 = P
summary(P34cT4)

P34cT4Co = unmark(split(P34cT4)$Co) # To subset the cocoa trees 
P34cTAsT = unmark(split(P34cT4)$AsT) # To subset Associated trees
P34cT4All = unmark(P34cT4) # Remove marks (tree type) from point pattern data


#----------------------------
#PLOT 29


Plot29_Data = read.csv("Plot29_tclass4.csv")
Plot29_Data$TTP = as.factor(Plot29_Data$TTP)

Plot29_T4  <- readShapePoints("Plot29_tclass4.shp")
SP <- as(Plot29_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot29_Data[5]

# Get Plot29 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot29.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot29 window (polygon) to trees in Plot29
P$window <- W
plot(P, pch=16, main="Plot29 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot29 All Trees")
plot(density(P))

P29T4 = P
summary(P29T4)
P29T4Co = unmark(split(P29T4)$Co) # To subset the cocoa trees 
P29TAsT = unmark(split(P29T4)$AsT) # To subset Associated trees
P29T4All = unmark(P29T4) # Remove marks (tree type) from point pattern data


#--------------------------------------------------------------------

#PLOT 29b

Plot29b_Data = read.csv("Plot29B_tclass4.csv")
Plot29b_Data$TTP = as.factor(Plot29b_Data$TTP)


Plot29b_T4  <- readShapePoints("Plot29B_tclass4.shp")
SP <- as(Plot29b_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot29b_Data[33]

# Get Plot29b shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot29b.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot29b windown (polygon) to trees in Plot29b
P$window <- W
a=plot(P, pch=16, main="Plot29b Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot29b")
# Get summary of Cocoa trees in plot 29b
P29bT4 = P
summary(P29bT4)

P29bT4Co = unmark(split(P29bT4)$Co) # To subset the cocoa trees 
P29bTAsT = unmark(split(P29bT4)$AsT) # To subset Associated trees
P29bT4All = unmark(P29bT4) # Remove marks (tree type) from point pattern data


#----------------------------
#PLOT 42


Plot42_Data = read.csv("Plot42_tclass4.csv")
Plot42_Data$TTP = as.factor(Plot42_Data$TTP)

Plot42_T4  <- readShapePoints("Plot42_tclass4.shp")
SP <- as(Plot42_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot42_Data[5]

# Get Plot42 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot42.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot42 window (polygon) to trees in Plot42
P$window <- W
plot(P, pch=16, main="Plot42 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot42 All Trees")
plot(density(P))

P42T4 = P
summary(P42T4)
P42T4Co = unmark(split(P42T4)$Co) # To subset the cocoa trees 
P42TAsT = unmark(split(P42T4)$AsT) # To subset Associated trees
P42T4All = unmark(P42T4) # Remove marks (tree type) from point pattern data

#--------------------------------------------------------------------

#PLOT 42b

Plot42b_Data = read.csv("Plot42B_tclass4.csv")
Plot42b_Data$TTP = as.factor(Plot42b_Data$TTP)

Plot42b_T4  <- readShapePoints("Plot42B_tclass4.shp")
SP <- as(Plot42b_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot42b_Data[33]

# Get Plot42b shapefile
S2  <- readShapePoly("D:/PhD Manuscript/......../Plot42b.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot42b windown (polygon) to trees in Plot42b
P$window <- W
a=plot(P, pch=16, main="Plot42b Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot42b")
# Get summary of Cocoa trees in plot 42b
P42bT4 = P
summary(P42bT4)

P42bT4Co = unmark(split(P42bT4)$Co) # To subset the cocoa trees 
P42bTAsT = unmark(split(P42bT4)$AsT) # To subset Associated trees
P42bT4All = unmark(P42bT4) # Remove marks (tree type) from point pattern data


#------------------------------
#PLOT 3


Plot3_Data = read.csv("Plot3_tclass4.csv")
Plot3_Data$TTP = as.factor(Plot3_Data$TTP)

Plot3_T4  <- readShapePoints("Plot3_tclass4.shp")
SP <- as(Plot3_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot3_Data[5]

# Get Plot3 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/......./Plot3.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot3 window (polygon) to trees in Plot3
P$window <- W
plot(P, pch=16, main="Plot3 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot3 All Trees")
plot(density(P))

P3T4 = P
summary(P3T4)
P3T4Co = unmark(split(P3T4)$Co) # To subset the cocoa trees 
P3TAsT = unmark(split(P3T4)$AsT) # To subset Associated trees
P3T4All = unmark(P3T4) # Remove marks (tree type) from point pattern data


#------------------------------------------------------------

#PLOT 3b

Plot3b_Data = read.csv("Plot3B_tclass4.csv")
Plot3b_Data$TTP = as.factor(Plot3b_Data$TTP)

Plot3b_T4  <- readShapePoints("Plot3B_tclass4.shp")
SP <- as(Plot3b_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot3b_Data[33]

# Get Plot3b shapefile
S2  <- readShapePoly("D:/PhD Manuscript/....../Plot3b.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot3b windown (polygon) to trees in Plot3b
P$window <- W
a=plot(P, pch=16, main="Plot3b Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot3b")
# Get summary of Cocoa trees in plot 3b
P3bT4 = P
summary(P3bT4)

P3bT4Co = unmark(split(P3bT4)$Co) # To subset the cocoa trees 
P3bTAsT = unmark(split(P3bT4)$AsT) # To subset Associated trees
P3bT4All = unmark(P3bT4) # Remove marks (tree type) from point pattern data


#----------------------------
#PLOT 10


Plot10_Data = read.csv("Plot10_tclass4.csv")
Plot10_Data$TTP = as.factor(Plot10_Data$TTP)

Plot10_T4  <- readShapePoints("Plot10_tclass4.shp")
SP <- as(Plot10_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot10_Data[5]

# Get Plot10 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/......../Plot10.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot10 window (polygon) to trees in Plot10
P$window <- W
plot(P, pch=16, main="Plot10 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot10 All Trees")
plot(density(P))

P10T4 = P
summary(P10T4)
P10T4Co = unmark(split(P10T4)$Co) # To subset the cocoa trees 
P10TAsT = unmark(split(P10T4)$AsT) # To subset Associated trees
P10T4All = unmark(P10T4) # Remove marks (tree type) from point pattern data


#--------------------------------------------------------------------

#PLOT 10b

Plot10b_Data = read.csv("Plot10B_tclass4.csv")
Plot10b_Data$TTP = as.factor(Plot10b_Data$TTP)

Plot10b_T4  <- readShapePoints("Plot10B_tclass4.shp")
SP <- as(Plot10b_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot10b_Data[33]

# Get Plot10b shapefile
S2  <- readShapePoly("D:/PhD Manuscript/......./Plot10b.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot10b windown (polygon) to trees in Plot10b
P$window <- W
a=plot(P, pch=16, main="Plot10b Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot10b")
# Get summary of Cocoa trees in plot 10b
P10bT4 = P
summary(P10bT4)

P10bT4Co = unmark(split(P10bT4)$Co) # To subset the cocoa trees 
P10bTAsT = unmark(split(P10bT4)$AsT) # To subset Associated trees
P10bT4All = unmark(P10bT4) # Remove marks (tree type) from point pattern data


#-------------------------------
#PLOT 11


Plot11_Data = read.csv("Plot11_tclass4.csv")
Plot11_Data$TTP = as.factor(Plot11_Data$TTP)

Plot11_T4  <- readShapePoints("Plot11_tclass4.shp")
SP <- as(Plot11_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot11_Data[5]

# Get Plot11 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/....../Plot11.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot11 window (polygon) to trees in Plot11
P$window <- W
plot(P, pch=16, main="Plot11 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot11 All Trees")
plot(density(P))

P11T4 = P
summary(P11T4)
P11T4Co = unmark(split(P11T4)$Co) # To subset the cocoa trees 
P11TAsT = unmark(split(P11T4)$AsT) # To subset Associated trees
P11T4All = unmark(P11T4) # Remove marks (tree type) from point pattern data


#--------------------------------------------------------------------

#PLOT 11b

Plot11b_Data = read.csv("Plot11B_tclass4.csv")
Plot11b_Data$TTP = as.factor(Plot11b_Data$TTP)

Plot11b_T4  <- readShapePoints("Plot11B_tclass4.shp")
SP <- as(Plot11b_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot11b_Data[33]

# Get Plot11b shapefile
S2  <- readShapePoly("D:/PhD Manuscript/....../Plot11b.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot11b windown (polygon) to trees in Plot11b
P$window <- W
a=plot(P, pch=16, main="Plot11b Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot11b")
# Get summary of Cocoa trees in plot 11b
P11bT4 = P
summary(P11bT4)

P11bT4Co = unmark(split(P11bT4)$Co) # To subset the cocoa trees 
P11bTAsT = unmark(split(P11bT4)$AsT) # To subset Associated trees
P11bT4All = unmark(P11bT4) # Remove marks (tree type) from point pattern data


#--------------------------------

#PLOT 12


Plot12_Data = read.csv("Plot12_tclass4.csv")
Plot12_Data$TTP = as.factor(Plot12_Data$TTP)

Plot12_T4  <- readShapePoints("Plot12_tclass4.shp")
SP <- as(Plot12_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot12_Data[5]

# Get Plot12 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/......./Plot12.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot12 window (polygon) to trees in Plot12
P$window <- W
plot(P, pch=16, main="Plot12 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot12 All Trees")
plot(density(P))

P12T4 = P
summary(P12T4)
P12T4Co = unmark(split(P12T4)$Co) # To subset just the cocoa trees 
P12TAsT = unmark(split(P12T4)$AsT) # To subset jut Associated trees
P12T4All = unmark(P12T4) # Remove marks (tree type) from point pattern data


#------------------------------------
#PLOT 2


Plot2_Data = read.csv("Plot2_tclass4.csv")
Plot2_Data$TTP = as.factor(Plot2_Data$TTP)

Plot2_T4  <- readShapePoints("Plot2_tclass4.shp")
SP <- as(Plot2_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot2_Data[5]

# Get Plot2 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/......./Plot2.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot2 window (polygon) to trees in Plot2
P$window <- W
plot(P, pch=16, main="Plot2 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot2 All Trees")
plot(density(P))

P2T4 = P
summary(P2T4)
P2T4Co = unmark(split(P2T4)$Co) # To subset the cocoa trees 
P2TAsT = unmark(split(P2T4)$AsT) # To subset Associated trees
P2T4All = unmark(P2T4) # Remove marks (tree type) from point pattern data


#--------------------------------------------------------------------

#PLOT 2b

Plot2b_Data = read.csv("Plot2B_tclass4.csv")
Plot2b_Data$TTP = as.factor(Plot2b_Data$TTP)

Plot2b_T4  <- readShapePoints("Plot2B_tclass4.shp")
SP <- as(Plot2b_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot2b_Data[33]

# Get Plot2b shapefile
S2  <- readShapePoly("D:/PhD Manuscript/....../Plot2b.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot2b windown (polygon) to trees in Plot2b
P$window <- W
a=plot(P, pch=16, main="Plot2b Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot2b")
# Get summary of Cocoa trees in plot 2b
P2bT4 = P
summary(P2bT4)

P2bT4Co = unmark(split(P2bT4)$Co) # To subset the cocoa trees 
P2bTAsT = unmark(split(P2bT4)$AsT) # To subset Associated trees
P2bT4All = unmark(P2bT4) # Remove marks (tree type) from point pattern data


#------------------------------
#PLOT 13


Plot13_Data = read.csv("Plot13_tclass4.csv")
Plot13_Data$TTP = as.factor(Plot13_Data$TTP)

Plot13_T4  <- readShapePoints("Plot13_tclass4.shp")
SP <- as(Plot13_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot13_Data[5]

# Get Plot13 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/....../Plot13.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot13 window (polygon) to trees in Plot13
P$window <- W
plot(P, pch=16, main="Plot13 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot13 All Trees")
plot(density(P))

P13T4 = P
summary(P13T4)
P13T4Co = unmark(split(P13T4)$Co) # To subset the cocoa trees 
P13TAsT = unmark(split(P13T4)$AsT) # To subset Associated trees
P13T4All = unmark(P13T4) # Remove marks (tree type) from point pattern data

#--------------------------------------------------------------------

#PLOT 13b

Plot13b_Data = read.csv("Plot13B_tclass4.csv")
Plot13b_Data$TTP = as.factor(Plot13b_Data$TTP)

Plot13b_T4  <- readShapePoints("Plot13B_tclass4.shp")
SP <- as(Plot13b_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot13b_Data[33]

# Get Plot13b shapefile
S2  <- readShapePoly("D:/PhD Manuscript/....../Plot13b.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot13b windown (polygon) to trees in Plot13b
P$window <- W
a=plot(P, pch=16, main="Plot13b Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot13b")
# Get summary of Cocoa trees in plot 13b
P13bT4 = P
summary(P13bT4)

P13bT4Co = unmark(split(P13bT4)$Co) # To subset the cocoa trees 
P13bTAsT = unmark(split(P13bT4)$AsT) # To subset Associated trees
P13bT4All = unmark(P13bT4) # Remove marks (tree type) from point pattern data


#------------------------------------
#PLOT 39


Plot39_Data = read.csv("Plot39_tclass4.csv")
Plot39_Data$TTP = as.factor(Plot39_Data$TTP)

Plot39_T4  <- readShapePoints("Plot39_tclass4.shp")
SP <- as(Plot39_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot39_Data[6]

# Get Plot39 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/....../Plot39.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot39 window (polygon) to trees in Plot39
P$window <- W
plot(P, pch=16, main="Plot39 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot39 All Trees")
plot(density(P))

P39T4 = P
summary(P39T4)
P39T4Co = unmark(split(P39T4)$Co) # To subset the cocoa trees 
P39TAsT = unmark(split(P39T4)$AsT) # To subset Associated trees
P39T4All = unmark(P39T4) # Remove marks (tree type) from point pattern data


#--------------------------------------------------------------------

#PLOT 39b

Plot39b_Data = read.csv("Plot39B_tclass4.csv")
Plot39b_Data$TTP = as.factor(Plot39b_Data$TTP)

Plot39b_T4  <- readShapePoints("Plot39B_tclass4.shp")
SP <- as(Plot39b_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot39b_Data[33]

# Get Plot39b shapefile
S2  <- readShapePoly("D:/PhD Manuscript/......./Plot39B.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot39b windown (polygon) to trees in Plot39b
P$window <- W
a=plot(P, pch=16, main="Plot39b Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot39b")
# Get summary of Cocoa trees in plot 39b
P39bT4 = P
summary(P39bT4)

P39bT4Co = unmark(split(P39bT4)$Co) # To subset the cocoa trees 
P39bTAsT = unmark(split(P39bT4)$AsT) # To subset Associated trees
P39bT4All = unmark(P39bT4) # Remove marks (tree type) from point pattern data


#----------------------------------
#PLOT 44


Plot44_Data = read.csv("Plot44_tclass4.csv")
Plot44_Data$TTP = as.factor(Plot44_Data$TTP)

Plot44_T4  <- readShapePoints("Plot44_tclass4.shp")
SP <- as(Plot44_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot44_Data[5]

# Get Plot44 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/......./Plot44.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot44 window (polygon) to trees in Plot44
P$window <- W
plot(P, pch=16, main="Plot44 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot44 All Trees")
plot(density(P))

P44T4 = P
summary(P44T4)
P44T4Co = unmark(split(P44T4)$Co) # To subset the cocoa trees 
P44TAsT = unmark(split(P44T4)$AsT) # To subset Associated trees
P44T4All = unmark(P44T4) # Remove marks (tree type) from point pattern data


#--------------------------------------------------------------------

#PLOT 44b

Plot44b_Data = read.csv("Plot44B_tclass4.csv")
Plot44b_Data$TTP = as.factor(Plot44b_Data$TTP)

Plot44b_T4  <- readShapePoints("Plot44B_tclass4.shp")
SP <- as(Plot44b_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot44b_Data[33]

# Get Plot44b shapefile
S2  <- readShapePoly("D:/PhD Manuscript/....../Plot44b.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot44b windown (polygon) to trees in Plot44b
P$window <- W
a=plot(P, pch=16, main="Plot44b Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot44b")
# Get summary of Cocoa trees in plot 44b
P44bT4 = P
summary(P44bT4)

P44bT4Co = unmark(split(P44bT4)$Co) # To subset the cocoa trees 
P44bTAsT = unmark(split(P44bT4)$AsT) # To subset Associated trees
P44bT4All = unmark(P44bT4) # Remove marks (tree type) from point pattern data


#--------------------------------------------------------------------

#PLOT 44c

Plot44c_Data = read.csv("Plot44C_tclass4.csv")
Plot44c_Data$TTP = as.factor(Plot44c_Data$TTP)

Plot44c_T4  <- readShapePoints("Plot44C_tclass4.shp")
SP <- as(Plot44c_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot44c_Data[33]

# Get Plot44c shapefile
S2  <- readShapePoly("D:/PhD Manuscript/....../Plot44C.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot44c windown (polygon) to trees in Plot44c
P$window <- W
a=plot(P, pch=16, main="Plot44c Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot44c")
# Get summary of Cocoa trees in plot 44c
P44cT4 = P
summary(P44cT4)

P44cT4Co = unmark(split(P44cT4)$Co) # To subset the cocoa trees 
P44cTAsT = unmark(split(P44cT4)$AsT) # To subset Associated trees
P44cT4All = unmark(P44cT4) # Remove marks (tree type) from point pattern data


#-----------------------------------

#PLOT 50


Plot50_Data = read.csv("Plot50_tclass4.csv")
Plot50_Data$TTP = as.factor(Plot50_Data$TTP)

Plot50_T4  <- readShapePoints("Plot50_tclass4.shp")
SP <- as(Plot50_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot50_Data[5]

# Get Plot50 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/....../Plot50.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot50 window (polygon) to trees in Plot50
P$window <- W
plot(P, pch=16, main="Plot50 All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot50 All Trees")
plot(density(P))

P50T4 = P
summary(P50T4)
P50T4Co = unmark(split(P50T4)$Co) # To subset the Cocoa trees 
P50TAsT = unmark(split(P50T4)$AsT) # To subset Associated trees
P50T4All = unmark(P50T4) # Remove marks (tree type) from point pattern data


#--------------------------------------------------------------------

#PLOT 50Aa

Plot50Aa_Data = read.csv("Plot50Aa_tclass4.csv")
Plot50Aa_Data$TTP = as.factor(Plot50Aa_Data$TTP)

Plot50Aa_T4  <- readShapePoints("Plot50Aa_tclass4.shp")
SP <- as(Plot50Aa_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot50Aa_Data[33]

# Get Plot50Aa shapefile
S2  <- readShapePoly("D:/PhD Manuscript/....../Plot50Aa.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot50Aa windown (polygon) to trees in Plot50Aa
P$window <- W
a=plot(P, pch=16, main="Plot50Aa Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot50Aa")
# Get summary of Cocoa trees in plot 50Aa
P50AaT4 = P
summary(P50AaT4)

P50AaT4Co = unmark(split(P50AaT4)$Co) # To subset the Cocoa trees 
P50AaTAsT = unmark(split(P50AaT4)$AsT) # To subset Associated trees
P50AaT4All = unmark(P50AaT4) # Remove marks (tree type) from point pattern data

#--------------------------------------------------------------------

#PLOT 50b

Plot50b_Data = read.csv("Plot50B_tclass4.csv")
Plot50b_Data$TTP = as.factor(Plot50b_Data$TTP)

Plot50b_T4  <- readShapePoints("Plot50B_tclass4.shp")
SP <- as(Plot50b_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot50b_Data[33]

# Get Plot50b shapefile
S2  <- readShapePoly("D:/PhD Manuscript/....../Plot50b.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot50b windown (polygon) to trees in Plot50b
P$window <- W
a=plot(P, pch=16, main="Plot50b Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot50b")
# Get summary of Cocoa trees in plot 50b
P50bT4 = P
summary(P50bT4)

P50bT4Co = unmark(split(P50bT4)$Co) # To subset the Cocoa trees 
P50bTAsT = unmark(split(P50bT4)$AsT) # To subset Associated trees
P50bT4All = unmark(P50bT4) # Remove marks (tree type) from point pattern data

#--------------------------------------------------------------------

#PLOT 50c

Plot50c_Data = read.csv("Plot50C_tclass4.csv")
Plot50c_Data$TTP = as.factor(Plot50c_Data$TTP)

Plot50c_T4  <- readShapePoints("Plot50C_tclass4.shp")
SP <- as(Plot50c_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot50c_Data[33]

# Get Plot50c shapefile
S2  <- readShapePoly("D:/PhD Manuscript/....../Plot50C.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot50c windown (polygon) to trees in Plot50c
P$window <- W
a=plot(P, pch=16, main="Plot50c Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot50c")
# Get summary of Cocoa trees in plot 50c
P50cT4 = P
summary(P50cT4)

P50cT4Co = unmark(split(P50cT4)$Co) # To subset the Cocoa trees 
P50cTAsT = unmark(split(P50cT4)$AsT) # To subset Associated trees
P50cT4All = unmark(P50cT4) # Remove marks (tree type) from point pattern data

#--------------------------------------------------------------------

#PLOT 45a

Plot45a_Data = read.csv("Plot45A_tclass4.csv")
Plot45a_Data$TTP = as.factor(Plot45a_Data$TTP)

Plot45a_T4  <- readShapePoints("Plot45A_tclass4.shp")
SP <- as(Plot45a_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot45a_Data[33]

# Get Plot45a shapefile
S2  <- readShapePoly("D:/PhD Manuscript/....../Plot45A.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot45a windown (polygon) to trees in Plot45a
P$window <- W
a=plot(P, pch=16, main="Plot45a Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot45a")
# Get summary of Cocoa trees in plot 45a
P45aT4 = P
summary(P45aT4)

P45aT4Co = unmark(split(P45aT4)$Co) # To subset the Cocoa trees 
P45aTAsT = unmark(split(P45aT4)$AsT) # To subset Associated trees
P45aT4All = unmark(P45aT4) # Remove marks (tree type) from point pattern data

#--------------------------------------------------------------------

#PLOT 48b

Plot48b_Data = read.csv("Plot48B_tclass4.csv")
Plot48b_Data$TTP = as.factor(Plot48b_Data$TTP)

Plot48b_T4  <- readShapePoints("Plot48B_tclass4.shp")
SP <- as(Plot48b_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot48b_Data[33]

# Get Plot48b shapefile
S2  <- readShapePoly("D:/PhD Manuscript/......./Plot48B.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot48b windown (polygon) to trees in Plot48b
P$window <- W
a=plot(P, pch=16, main="Plot48b Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot48b")
# Get summary of Cocoa trees in plot 48b
P48bT4 = P
summary(P48bT4)

P48bT4Co = unmark(split(P48bT4)$Co) # To subset the Cocoa trees 
P48bTAsT = unmark(split(P48bT4)$AsT) # To subset Associated trees
P48bT4All = unmark(P48bT4) # Remove marks (tree type) from point pattern data


#--------------------------------------------------------------------

#PLOT 49b

Plot49b_Data = read.csv("Plot49B_tclass4.csv")
Plot49b_Data$TTP = as.factor(Plot49b_Data$TTP)

Plot49b_T4  <- readShapePoints("Plot49B_tclass4.shp")
SP <- as(Plot49b_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot49b_Data[33]

# Get Plot49b shapefile
S2  <- readShapePoly("D:/PhD Manuscript/....../Plot49B.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot49b windown (polygon) to trees in Plot49b
P$window <- W
a=plot(P, pch=16, main="Plot49b Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot49b")
# Get summary of Cocoa trees in plot 49b
P49bT4 = P
summary(P49bT4)

P49bT4Co = unmark(split(P49bT4)$Co) # To subset the Cocoa trees 
P49bTAsT = unmark(split(P49bT4)$AsT) # To subset Associated trees
P49bT4All = unmark(P49bT4) # Remove marks (tree type) from point pattern data

#--------------------------------------------------------------------

#PLOT 60a

Plot60a_Data = read.csv("Plot60A_tclass4.csv")
Plot60a_Data$TTP = as.factor(Plot60a_Data$TTP)

Plot60a_T4  <- readShapePoints("Plot60A_tclass4.shp")
SP <- as(Plot60a_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot60a_Data[33]

# Get Plot60a shapefile
S2  <- readShapePoly("D:/PhD Manuscript/....../Plot60A.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot60a windown (polygon) to trees in Plot60a
P$window <- W
a=plot(P, pch=16, main="Plot60a Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, Plot60a")
# Get summary of Cocoa trees in plot 60a
P60aT4 = P
summary(P60aT4)

P60aT4Co = unmark(split(P60aT4)$Co) # To subset the Cocoa trees 
P60aTAsT = unmark(split(P60aT4)$AsT) # To subset Associated trees
P60aT4All = unmark(P60aT4) # Remove marks (tree type) from point pattern data

#--------------------------------------------------------------------

#PLOT NS1

PlotNS1_Data = read.csv("PlotNS1_tclass4.csv")
PlotNS1_Data$TTP = as.factor(PlotNS1_Data$TTP)

PlotNS1_T4  <- readShapePoints("PlotNS1_tclass4.shp")
SP <- as(PlotNS1_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = PlotNS1_Data[33]

# Get PlotNS1 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/....../PlotNS1.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply PlotNS1 windown (polygon) to trees in PlotNS1
P$window <- W
a=plot(P, pch=16, main="PlotNS1 Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, PlotNS1")
# Get summary of Cocoa trees in plot NS1
PNS1T4 = P
summary(PNS1T4)

PNS1T4Co = unmark(split(PNS1T4)$Co) # To subset the Cocoa trees 
PNS1TAsT = unmark(split(PNS1T4)$AsT) # To subset Associated trees
PNS1T4All = unmark(PNS1T4) # Remove marks (tree type) from point pattern data

#--------------------------------------------------------------------

#PLOT OS3

PlotOS3_Data = read.csv("PlotOS3_tclass4.csv")
PlotOS3_Data$TTP = as.factor(PlotOS3_Data$TTP)

PlotOS3_T4  <- readShapePoints("PlotOS3_tclass4.shp")
SP <- as(PlotOS3_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = PlotOS3_Data[33]

# Get PlotOS3 shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../PlotOS3.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply PlotOS3 windown (polygon) to trees in PlotOS3
P$window <- W
a=plot(P, pch=16, main="PlotOS3 Cocoa Trees")
#dev.copy(tiff,'Plot23 Trees.tiff')
#dev.off()
plot(density(P), main="Kernel intensity, PlotOS3")
# Get summary of Cocoa trees in plot OS3
POS3T4 = P
summary(POS3T4)

POS3T4Co = unmark(split(POS3T4)$Co) # To subset the Cocoa trees 
POS3TAsT = unmark(split(POS3T4)$AsT) # To subset Associated trees
POS3T4All = unmark(POS3T4) # Remove marks (tree type) from point pattern data


#------------------------------

# Compute Hyperframe dataset for COAF1 - Assoiated Trees

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

# PLOTTING THE POOLED L-FUNCTION FOR SPATIAL POINT PATTERN


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




