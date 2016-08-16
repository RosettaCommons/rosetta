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
id = "angle_kernels",
author = "Matthew O'Meara",
brief_description = "Compare the effects of kernel on estimating the angular projection of a 3d distribution",
long_description = "

",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


# the default bin width does funny things when sd(x) = 0
bw <- .3

test_points <-
	data.frame(
		theta=rep(
			c(10, 20, 30, 40, 50, 60, 70, 80), 100))
test_points$theta_factor <- factor(test_points$theta)
test_points$theta_rad <- test_points$theta * pi/180
test_points$cos_theta <- cos(test_points$theta_rad)

print(summary(test_points))

dens_cos_cos <-
	estimate_density_1d(
		test_points, c("theta_factor"),
		"cos_theta",
		bw=bw,
		sample_domain=c(-1,1),
		refrect_right=FALSE, right_boundary = 1)
dens_cos_cos$method <- "estimate cos(theta), plot cos(theta)"

dens_cos_theta <-
	estimate_density_1d(
		test_points, c("theta_factor"),
		"cos_theta",
		bw=bw,
		sample_domain=c(-1,1),
		refrect_right=FALSE, right_boundary = 1)
dens_cos_theta$method <- "estimate cos(theta), plot theta"
dens_cos_theta$x <- acos(dens_cos_theta$x) * 180/pi



dens_theta_cos <-
	estimate_density_1d(
		test_points, c("theta_factor"),
		"theta",
		bw=bw*50,
		sample_domain=c(0,180))
dens_theta_cos$x <- cos(dens_theta_cos$x *pi/180)
dens_theta_cos$method <- "estimate theta, plot cos(theta)"

dens_theta_theta <-
	estimate_density_1d(
		test_points, c("theta_factor"),
		"theta",
		bw=bw*50,
		sample_domain=c(0,180))
dens_theta_theta$method <- "estimate theta, plot theta"


dens_theta_cos_weight_theta <-
	estimate_density_1d(
		test_points,
		c("theta_factor"),
		"theta_rad",
		weight_fun=conical_3d_normalization,
		bw=bw,
		sample_domain=c(0,pi),
		refrect_left=FALSE, left_boundary = 0)
dens_theta_cos_weight_theta$method <- "estimate theta, weight cos(theta), plot theta"
dens_theta_cos_weight_theta$x <- dens_theta_cos_weight_theta$x * 180/pi

dens_theta_cos_weight_cos <-
	estimate_density_1d(
		test_points,
		c("theta_factor"),
		"theta_rad",
		weight_fun=conical_3d_normalization,
		bw=bw,
		sample_domain=c(0,pi),
		refrect_left=FALSE, left_boundary = 0)
dens_theta_cos_weight_cos$method <- "estimate theta, weight cos(theta), plot cos(theta)"
dens_theta_cos_weight_cos$x <- cos(dens_theta_cos_weight_cos$x)


dens <- rbind(
	dens_cos_cos,
	dens_cos_theta,
	dens_theta_cos,
	dens_theta_theta,
	dens_theta_cos_weight_cos,
	dens_theta_cos_weight_theta)


plot_id <- "angle_kernels"
p <- ggplot(dens) + theme_bw() +
	geom_line(aes(x, y, color=theta_factor, group=theta_factor)) +
	facet_wrap(~method, scales="free", ncol=2) +
	ggtitle("Kernels for Angle Distributions")

save_plots(self, plot_id, sample_sources, output_dir, output_formats)




})) # end FeaturesAnalysis
