# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

check_setup()
feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "hbond_null",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

n_pts=100000

n_pts_ball=integer(n_pts * 6/pi) # fraction of volume in ball

null.cube <- data.frame(x=runif(n_pts_ball)-.5, y=runif(n_pts_ball)-.5, z=runif(n_pts_ball)-.5)
null.cube <- transform(null.cube, length = sqrt(x^2 + y^2 + z^2))

null.ball <- null.cube[null.cube$length <= .5,]
null.ball$id <- "ball_null"


dens_raw <- estimate_density_1d(
  data = null.ball,
  id=c("id"),
  variable = "length",
  histogram=TRUE)
dens_raw$normalization <- "Raw"

dens_normalized <- estimate_density_1d(
  data = null.ball,
  id=c("id"),
  variable = "length",
  histogram=TRUE,
  weight_fun = radial_3d_normalization)
dens_normalized$normalization <- "Normalized"

dens <- rbind(dens_raw, dens_normalized)

plot_id = "3D_radial_length_null_model"
ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=normalization)) +
	geom_indicator(aes(indicator=counts, colour=normalization, group=normalization)) +
	ggtitle("3D Radial Length null Model by Normalization.") +
	labs(x=expression(paste('Distance to Origin')),
	     y="FeatureDensity")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


# refine ggplot2::plotmatrix to use smaller sized points
# and use binned rather than smoothed densities
plotmatrix <- function (data, mapping = aes(), colour = "black")
{
    grid <- expand.grid(x = 1:ncol(data), y = 1:ncol(data))
    grid <- subset(grid, x != y)
    all <- do.call("rbind", lapply(1:nrow(grid), function(i) {
        xcol <- grid[i, "x"]
        ycol <- grid[i, "y"]
        data.frame(xvar = names(data)[ycol], yvar = names(data)[xcol],
            x = data[, xcol], y = data[, ycol], data)
    }))
    all$xvar <- factor(all$xvar, levels = names(data))
    all$yvar <- factor(all$yvar, levels = names(data))
    densities <- do.call("rbind", lapply(1:ncol(data), function(i) {
        data.frame(xvar = names(data)[i], yvar = names(data)[i],
            x = data[, i])
    }))
    mapping <- plyr::defaults(mapping, aes_string(x = "x", y = "y"))
    class(mapping) <- "uneval"
    ggplot(all, mapping) + facet_grid(xvar ~ yvar, scales = "free") +
        geom_point(colour = colour, na.rm = TRUE, size=.3) + stat_bin(aes(x = x,
        y = ..ndensity.. * diff(range(x)) + min(x)), data = densities,
        position = "identity", colour = "grey20", geom = "line")
}



plot_id = "3D_ball_null_model"
plotmatrix(data=null.ball[,c("x", "y", "z")]) + theme_bw() +
	ggtitle("3D Ball Null Model.") +
	coord_equal(ratio=1)
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

n_pts <- 100000
gaussian.ball <- data.frame(x=rnorm(n_pts),y=rnorm(n_pts), z=rnorm(n_pts))

null.sphere<-transform(gaussian.ball,
  x=(x/sqrt(x^2+y^2+z^2)),
  y=(y/sqrt(x^2+y^2+z^2)),
  z=(z/sqrt(x^2+y^2+z^2)))

null.sphere$id <- "spherical_null"
null.sphere$angle <- pi - acos(null.sphere$x)
null.sphere <- null.sphere[null.sphere$angle < pi/2,]
dens <- estimate_density_1d(
  data = null.sphere,
  id=c("id"),
  variable = "angle",
	weight_fun=conical_3d_normalization)

plot_id = "hbond_AHD_spherical_null_model"
ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y)) +
	geom_indicator(aes(indicator=counts)) +
	ggtitle("Cosine of spherical Angle Null Model") +
	labs(x=expression(paste('Cosine(Central Angle)')),
	     y="FeatureDensity")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


n_pts=100000
null.cube <- data.frame(x=10*runif(n_pts)-5, y=10*runif(n_pts)-5, z=10*runif(n_pts)-5)

null.ball <- transform(null.cube,
	r=sqrt(x^2 + y^2 + z^2),
	cosTHETA=z/sqrt(x^2 + y^2 + z^2))
null.ball <- null.ball[null.ball$r < 5,]

plot_id = "hbond_cosBAH_vs_r_cubed_spherical_null_model"
ggplot(data=null.ball) + theme_bw() +
	geom_point(aes(x=cosTHETA, y=r^3)) +
	ggtitle("cosTHETA vs r^3 of Spherical Null Model") +
	labs(x=expression(paste('Cos(THETA)')),
	     y="r^3")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)





})) # end FeaturesAnalysis
