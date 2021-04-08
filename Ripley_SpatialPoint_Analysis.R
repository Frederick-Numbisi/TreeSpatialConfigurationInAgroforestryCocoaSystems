
library(sp)
library(spatstat)
library(maptools)
library(ggplot2)

setwd("D:/...../...../Bakoa_TreesClass4") # Select location of dataset as the working directory

# -----------------------------------------
# Analysis procedure for data in each plot
# -----------------------------------------
# Sample code for plot

Plot_Data = read.csv("Plot_tclass4.csv") # Read the data sete
Plot_Data$TTP = as.factor(Plot23_Data$TTP) # Convert the TTP variable to factor 

Plot_T4  <- readShapePoints("Plot_tclass4.shp") # Read the spatial data (shapefile) of tree point locations
SP <- as(Plot_T4, "SpatialPoints") # Convert the spatial data to spatial point patterns
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot_Data[5]

# Get Plot shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot23.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot window (polygon) to trees in the Plot
P$window <- W
plot(P, pch=16, main="Plot All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot All Trees")
plot(density(P))
#legend(0.25, 0.5, names(a), pch=a)

#Set measurement unit of window()
#unitname(P) <- "meter"

## The Kdot(r) inhomogenous function for each tree group within stand

# method (1): estimate intensities by nonparametric smoothing
Co <- split(P)$Co
AsT <- split(P)$AsT
All = unmark(P)
lambdaAll = density.ppp(All, sigma=0.15, at ="points")
lambdaCo <- density.ppp(Co, sigma=0.15, at="points")
lambdaAsT <- density.ppp(AsT, sigma=0.15, at="points")


# Fit model (second order quadratic distribution)
fit <- ppm(P ~ marks * polynom(x,y,2))
L_AsT_plot <- Ldot.inhom(P, "AsT", lambdaX=fit, update=FALSE)
plot(Lij_Plot, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot Associated Trees", legendpos="bottom")
L_AsT_plot$iso  <- L_AsT_plot$iso  - L_AsT_plot$r
L_AsT_plot$theo <- L_AsT_plot$theo - L_AsT_plot$r
write.csv(L_AsT_plot23, file = "Li_AsT_Plot.csv", row.names=TRUE)

#Li (nonhomogenous) for Cocoa Trees2 in stand 

L_Co2_plot <- Ldot.inhom(P, "Co", lambdaX=fit, update=FALSE) #using fit computed for all tree

plot(Lij_Plot, . - r ~ r, xlab="Distance (m)", ylab="Li(r)", main="Li inhom. Function, Plot Cocoa Trees3", legendpos="bottom")
L_AsT_plot$iso  <- L_AsT_plot$iso  - L_AsT_plot$r
L_AsT_plot$theo <- L_AsT_plot$theo - L_AsT_plot$r
write.csv(L_AsT_plot, file = "Li_Co3_Plot.csv", row.names=TRUE)


#Bivariate(cross) inhomogenous Kij(r) Ripley function
#The Kcross(r) test, which examines distribution a tree group relative to another group within stand

# method (3): fit parametric intensity model
fit <- ppm(P, ~marks * polynom(x,y,2))

# evaluate fitted intensities at data points
# (these are the intensities of the sub-processes of each type)
inten <- fitted(fit, dataonly=TRUE)

# split according to types of points
lambda2 <- split(inten, P$marks)
Kij_plot <- Kcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))
plot(Kij_Plot, xlab="Distance (m)", main="Kij (2nd order moment) Co3 to Ast")

write.csv(Kij_Plot23, file = "Kij_bivariate_inhom_Plot_Co3-AsT.csv", row.names=TRUE)



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
Lij_Plot <- Lcross.inhom(P, "Co", "AsT", lambda2$Co, lambda2$AsT, correction = c("border", "isotropic", "Ripley", "translate"))

plot(Lij_Plot, . - r ~ r, xlab="Distance (m)", ylab="Lij(r)", main="Lij inhom. Function, Plot23 Cocoa3-AsT", legendpos="bottom")


#We can make permanent changes to the Lij object's values as follows:
Lij_Plot$iso  <- Lij_Plot$iso  - Lij_Plot$r
Lij_Plot$theo <- Lij_Plot$theo - Lij_Plot$r

write.csv(Lij_Plot, file = "Lij_Bivariate_inhom_Plot23_Co3_to_AsT.csv", row.names=TRUE)



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
# CREATE HYPERFRAME POINT PATTERN SPATIAL DATA FOR COCOA AND ASSOCIATED TREES
#
#-------------------------------------------------------------------------

library(sp)
library(spatstat)
library(maptools)
library(ggplot2)

setwd("/Path/to/the/...../Working/Directory")

# The analysis procedure is replicated for all plots using a function call
# Case analysis for a Plot

Plot_Data = read.csv("Plot_tclass4.csv")
Plot_Data$TTP = as.factor(Plot_Data$TTP)

Plot_T4  <- readShapePoints("Plot_tclass4.shp")
SP <- as(Plot_T4, "SpatialPoints")
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot_Data[5]

# Get Plot shapefile
S2  <- readShapePoly("/Directory/to/...../Plot.shp")
SP2 <- as(S2, "SpatialPolygons")
W   <- as(SP2, "owin")

# Apply Plot window (polygon) to trees in Plot
P$window <- W
plot(P, pch=16, main="Plot All Trees")
summary(P)
plot(split(P))
plot(density(split(P)))
a = plot(P, main="Plot All Trees")
plot(density(P))

PT4 = P
summary(PT4)

PT4Co = unmark(split(PT4)$Co) # To subset just the cocoa trees 
PT4AsT = unmark(split(PT4)$AsT) # To subset jut Associated trees
PT4All = unmark(PT4) # Remove marks (tree type) from point pattern data




