
library(sp)
library(spatstat)
library(maptools)
library(ggplot2)

setwd("Path/to/the/working/Directory") # Select location of dataset as the working directory

# -----------------------------------------
# Analysis procedure for data in each plot
# -----------------------------------------
# Sample code for analysing trees in a sample plot

Plot_Data = read.csv("Plot_tclass4.csv") # Read the data sete
Plot_Data$TTP = as.factor(Plot_Data$TTP) # Convert the TTP variable to factor 

Plot_T4  <- readShapePoints("Plot_tclass4.shp") # Read the spatial data (shapefile) of tree point locations
SP <- as(Plot_T4, "SpatialPoints") # Convert the spatial data to spatial point patterns
P  <- as(SP, "ppp")
#Add marks to P(data)
marks(P) = Plot_Data[5]

# Get Plot shapefile
S2  <- readShapePoly("D:/PhD Manuscript/...../Plot.shp")
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



## Compute the Kdot(r) inhomogenous function for each tree group within stand

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

write.csv(Lij_Plot, file = "Lij_Bivariate_inhom_Plot_Co3_to_AsT.csv", row.names=TRUE)



#-----------------------------------------------
#
# SAVING AND APPENDING THE Ripley-K RESULTS AND FILES IN CSV
#--------------------------------------------------

#Appending Multivariate Lij(r) Results
file_path_dir <- "/path/to/the/directory/of/Lij(r)Results/"

#Create function to merge results
RipleyKmerge <- function(mydirectory){
  filenames <- list.files(path = mydirectory, full.names = TRUE)
  Lijlist <- lapply(filenames, function(x){
    read.csv(file = x, header = TRUE)}
                   )
  Reduce(function(x,y) {
    merge(x, y, all = TRUE)}, Lijlist)
}

# Run function using directory of results
Lij_Results <- RipleyKmerge(file_path_dir)
write.csv(Lij_Results, "Lij_Results_Cocoa-AsT.csv", row.names = TRUE)





