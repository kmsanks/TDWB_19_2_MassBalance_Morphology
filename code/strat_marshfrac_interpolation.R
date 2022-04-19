# The purpose of this script is to interpolate the marsh fraction across the experimental delta top using Bayesian kriging.
# This script is supplemental to the results of the study "Marsh sedimentation controls delta top morphology, slope, and mass balance." submitted to Geophysical Research Letters in 2022.
# This script will produce results described in the Supporting Information (SI) Text 2.1 and SI Figure 1.
# Please source this script, do not run and install any packages needed before sourcing.

#-------------------------------
# Clean all the stored variables
#-------------------------------
rm(list = ls())
t0 <- Sys.time(); t0

#-----------------------
# Set working directory
#-----------------------
setwd(getSrcDirectory()[1])
#Check working directory using
getwd() #Check the code location is correct

#-------------------
# Load the packages
#-------------------
#install any packages if needed (example install sp using install.packages('sp'))
library(sp);library(maptools);library(spatstat);library(stats);library(MASS);library(splancs) 
library(scatterplot3d);library(RandomFields);library(methods);library(geoR);library(lattice)
library(rgdal);library(gstat);library(RSAGA);library(plotrix)

#-------------------------
# Import data and border
#-------------------------
strat <- read.csv("../data/strat_data.csv", header = T); head(strat, 3); names(strat)
names(strat) <- c('x_coord', 'y_coord', 'marsh fraction', 'marsh thickness')

#-------------------------
# Create prediction grid
#-------------------------
strat.grid <- pred_grid(c(0,800), c(0,800), by=1)

#-------------------------------
# Prepare Covariates and Geodata
#-------------------------------
strat.geodata <- as.geodata(strat, coords.col = 1:2, data.col = 3)
coords <- strat[1:2]
coordinates(coords) <-~x_coord + y_coord

#Variograms and models were created using trends= cte, 1st, and 2nd

#-------------------------------
# Variogram for marsh fraction in stratigraphy
#-------------------------------
strat.variog <- variog(trend = "cte", strat.geodata, uvec = seq(0, 300, 10), lambda = 0)
svg(filename='../figures/marsh_frac_variog.svg')
plot(strat.variog)
title(main = 'Semivariogram for Marsh Fraction (-)')
dev.off()

#-----------------------
# REML Accretion Rate
#-----------------------
set.seed(12345)
strat.var.reml <- likfit(strat.geodata, coords = strat.geodata$coords, data = strat.geodata$data, trend = 'cte' , ini.cov.pars = c(0.6, 30), fix.nugget = F, 
                       nugget = 0, lambda = 1, kappa = 0.5, psiA = 0, psiR = 1, lik.method = "REML", components = F)

strat.env.model <- variog.model.env(strat.geodata, obj.var = strat.variog, model = strat.var.reml)
strat.var.reml; summary(strat.var.reml)

beta0 <- strat.var.reml$beta; names(beta0) <- NULL; beta0
tausq <- strat.var.reml$nugget; tausq
sigmasq <- strat.var.reml$sigmasq; sigmasq
nsratio <- tausq/sigmasq; nsratio
phi <- strat.var.reml$phi; phi
range <- strat.var.reml$practicalRange; range
lambda <- strat.var.reml$lambda; lambda
AIC.ns <- strat.var.reml$nospatial$AIC.ns; AIC.ns
AIC.s <- strat.var.reml$AIC; AIC.s
DAIC <- AIC.ns-AIC.s; DAIC

#---------
# Kriging
#---------
set.seed(268)
strat.MC <- model.control(trend.d = "cte", trend.l = "cte", kappa = 0.5, lambda = lambda)
strat.PC <- prior.control(sigmasq.prior = "sc.inv.chisq", df.sigmasq = 5, sigmasq = seq(0, 10, l = 21), phi.prior = "unif", 
                        phi.discrete = seq(0, 500, l = 21), tausq.rel.prior = "unif", tausq.rel.discrete = seq(0, 1, l = 21))
strat.OC <- output.control(n.post = 1000, simulations.predictive = F, moments = T)
strat.kb <- krige.bayes(strat.geodata, loc = strat.grid, model = strat.MC, prior = strat.PC, output = strat.OC)
strat.kb.summary <- function(x){quantile(x, prob = c(0.5, 0.025, 0.975))}

#------------------
# LOO-X-validation 
#------------------
strat.var.valid <- xvalid(strat.geodata, model = strat.var.reml, reest = F)
pdf('../figures/strat_xvalid.pdf', pointsize=8, bg= "transparent")
par(mfrow=c(5,2))
plot(strat.var.valid)
#title(main = 'Cross Validation for stratigraphic marsh fraction (-)')
dev.off()

#tiff(filename = paste("../figures/strat/cross-validation-",tempvar,".tif"), width = 5, height = 5, 
#units = "in", pointsize = 12, bg = "transparent", res = 300, restoreConsole = T)
par(mfcol = c(5, 2), mar = c(2.3, 2.3, 0.5, 0.5), mgp = c(1.3, 0.6, 0))
plot(strat.var.valid)	
dev.off()

#tiff(filename = paste("../figures/strat/validate-",tempvar,".tif"), width = 5, height = 5, 
#units = "in", pointsize = 12, bg = "transparent", res = 300, restoreConsole = T)
plot(strat.var.valid$data, strat.var.valid$predicted, xlim = c(0, 0.5), ylim = c(0, 0.5), 
     xlab = "Observed", ylab = "Predicted", pch = 21, bg = "grey", cex = 1.2, 
     cex.axis = 1.2, cex.lab = 1.3, tck = 0.03)
strat.lm.var <- lm(strat.var.valid$predicted~strat.var.valid$data); abline(strat.lm.var)	
dev.off()

#-------------------
#Parameter summary
#-------------------
p <- c(0.5, 0.025, 0.05, 0.95, 0.975)
beta.q <-quantile(strat.kb$posterior$sample$beta, p)
phi.q <-quantile(strat.kb$posterior$sample$phi, p)
sigmasq.q <-quantile(strat.kb$posterior$sample$sigmasq, p)
tausq.rel.q <-quantile(strat.kb$posterior$sample$tausq.rel, p)
tausq <- (strat.kb$posterior$sample$tausq.rel)*(strat.kb$posterior$sample$sigmasq) 
tausq.q <- quantile(tausq, p)

beta.bayes <-strat.kb$posterior$sample$beta
sigmasq.bayes <-strat.kb$posterior$sample$sigmasq
tausq.bayes <- tausq
nsratio.bayes <- strat.kb$posterior$sample$tausq.rel
phi.bayes <- strat.kb$posterior$sample$phi

samples.pars  <- cbind(beta.bayes, sigmasq.bayes, tausq.bayes, nsratio.bayes, phi.bayes)
write.table(samples.pars, file = "../data/interpolation/L/samples.pars.txt", append = T, quote = T, sep = " ", eol = "\n", na = "NA", dec = ".", 
            row.names = F, col.names = T, qmethod = c("escape", "double"))

summary.pars <- cbind(beta.q, phi.q, sigmasq.q, tausq.q, tausq.rel.q, strat.R2)
write.table(summary.pars, file = "../data/interpolation/L/summary.txt", append = T, quote = T, sep = " ", eol = "\n", na = "NA", dec = ".", 
            row.names = F, col.names = T, qmethod = c("escape", "double"))

#---------------------------------
# Create maps of mean and variance
#---------------------------------
lai.colors <- colorRampPalette(c("dark orange", "green", "dark green", "purple"))

B.matrix <- data.matrix(strat.grid, rownames.force = NA)
lai.colors <- colorRampPalette(c("dark orange", "green", "dark green", "purple"))
pdf('../figures/strat_interpolation.pdf', bg= "transparent")
image(strat.kb, xlim = c(0,800), ylim = c(0,800),zlim = c(0,1),loc = strat.grid, col=lai.colors(10), xlab="Coord X", ylab="Coord Y")
dev.off()

grd <- strat.grid*1000
output = rasterFromXYZ(cbind(strat.grid, strat.kb$predictive$mean))#,crs = "+proj=utm +north +zone=15 +datum=WGS84")  
lai.colors <- colorRampPalette(c("dark orange", "green", "dark green", "purple"))
pdf('../figures/strat_interpolation.pdf', bg= "transparent")
image(output) 
dev.off()

x <- writeRaster(output, '../data/interpolation/strat_prediction.tif', overwrite=TRUE)

#Save interpolation values and coordinates
df2 <- data.frame(rasterToPoints(output))
df2.x <- write.csv(df2, '../data/interpolation/strat_prediction.csv')

#Variance
output_var = rasterFromXYZ(cbind(strat.grid, strat.kb$predictive$variance))#,crs = "+proj=utm +north +zone=15 +datum=WGS84")  
lai.colors <- colorRampPalette(c("dark orange", "green", "dark green", "purple"))
pdf('../figures/strat_variance.pdf', bg= "transparent")
image(output_var) 
dev.off()

x <- writeRaster(output_var, '../data/interpolation/strat_variance.tif', overwrite=TRUE)

#Save variance values and coordinates
df2.var <- data.frame(rasterToPoints(output_var))
df2.x.var <- write.csv(df2.var, '../data/interpolation/strat_variance.csv')

#Uncertainty
output_unc_dis = rasterFromXYZ(cbind(strat.grid, strat.kb$predictive$distribution))#,crs = "+proj=utm +north +zone=15 +datum=WGS84")  
lai.colors <- colorRampPalette(c("dark orange", "green", "dark green", "purple"))
pdf('../figures/strat_uncertain_distribution.pdf', bg= "transparent")
image(output_unc) 

dev.off()

x <- writeRaster(output_unc_dis, '../data/interpolation/strat_uncertain_distribution.tif', overwrite=TRUE)

#Save uncertainty values and coordinates
df2.unc <- data.frame(rasterToPoints(output_unc_dis))
df2.x.unc <- write.csv(df2.unc, '../data/interpolation/strat_uncertain_distribution.csv')

tiff(filename = "../figures/strat.krige.tif", width = 5, height = 5, units = "in", pointsize = 12, bg = "transparent", res = 300, restoreConsole = T)
image(strat.kb, xlim = c(0, 800), ylim = c(0, 800), zlim = c(0, 1), xlab = " XCoord", ylab = "YCoord", 
      col = lai.colors(100), tck = 10, xaxt = "n", yaxt = "n")
contour(strat.kb, levels = seq(0, 1, by = 0.1), drawlabels = F, lwd = 0.3, col = "dark green", labcex = 0.5, method = "flattest", add = T)
axis(1, at = seq(0, 800, by = 10), tck = 5)
axis(2, at = seq(0, 800, by = 10), tck = 5)
legend.krige(x.leg = c(0, 800), y.leg = c(0, 800), values = c(0, 1), col = lai.colors(100), vertical = F, offset.leg = 2)
dev.off()

tiff(filename = "../figures/strat.krige.var.tif", width = 5, height = 5, units = "in", pointsize = 12, bg = "transparent", res = 300, restoreConsole = T)
image(strat.kb, val = "variance", xlim = c(0, 800), ylim = c(0, 800), zlim = c(0, 1), xlab = " XCoord", ylab = "YCoord", col = lai.colors(100), tck = 10, xaxt = "n", yaxt = "n")
contour(strat.kb, val = "variance", levels = seq(0, 1, by = 0.1), drawlabels = F, lwd = 0.3, col = "brown", labcex = 0.5, method = "flattest", add = T)
axis(1, at = seq(0, 800, by = 10), tck = 5)
axis(2, at = seq(0, 800, by = 10), tck = 5)
legend.krige(x.leg = c(0, 800), y.leg = c(0, 800), values = c(0, 1), col = lai.colors(100), vertical = F, offset.leg = 2)
dev.off()
#---------------------------------
# Plot parameter distribution
#---------------------------------
tiff(filename = "../figures/strat.krige.parms.tif", width = 5, height = 5, units = "in", pointsize = 12, bg = "transparent", res = 300, restoreConsole = T)
par(mfrow = c(2, 2), mar = c(3, 3, 1, 0.5), mgp = c(2, 1, 0))
hist(beta.bayes, col = "dark green", main = "", xlab = expression(beta,0), prob = T); rug (beta.bayes, col = 4); lines(density (beta.bayes))
hist(phi.bayes, col = "green", main = "", xlab = expression(phi), prob = T); rug (phi.bayes, col = 4); lines(density (phi.bayes)) 
hist(sigmasq.bayes, col = "yellow", main = "", xlab = expression(sigma^2), prob = T); rug (sigmasq.bayes, col = 4); lines(density (sigmasq.bayes))
hist(tausq.bayes, col = "red", main = "", xlab = expression(tau^2), prob = T); rug (tausq.bayes, col = 4); lines(density (tausq.bayes))
dev.off()

tiff(filename = "../figures/strat.krige.parmspp.tif", width = 10, height = 5, units = "in", pointsize = 12, bg = "transparent", res = 300, restoreConsole = T)
par(mfrow = c(1,2), mar = c(3, 3, 1, 0.5), mgp = c(2, 1, 0))
plot(strat.kb, col = c("white", "dark green"))
dev.off()

#------------------
# Plot variogram
#------------------
tiff(filename = "../figures/strat.variogram.tif", width = 5, height = 5, units = "in", pointsize = 12, bg = "transparent", res = 300, restoreConsole = T)
plot(strat.variog, pch = 21, bg = "dark green", cex.axis = 1.2, cex.lab = 1.2, xlab = "Distance", ylab = "Semivariance", tck = 0.01)
lines(strat.kb, max.dist = 300, summ = "median", post = "par", lty = 5, lwd = 2, col = "grey")
lines(strat.kb, max.dist = 300, summ = "mode", post = "par", lty = 3, lwd = 2)
lines(strat.kb, summ = strat.kb.summary, ty="l", lty=c(1,1,1), lwd = c(2, 1, 1), col= c("black", "blue","blue"))
legend("bottomright", legend = c( "Post. Mean", "Post. Median", "Post. Mode", "95% CI"), lty = c(1, 5, 3, 1), lwd = c(2, 2, 2, 1), col = c("black", "grey", "black", "blue"), cex = 0.8)
dev.off()

#-----------
#Save output
#-----------
save(list = ls(all = T), file = "../data/interpolation/strat.RData")

print(Sys.time()-t0)

#----
#END#
#----