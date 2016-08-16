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
id = "OHdonor",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


extract_transform_features <- function(sample_sources){
	sele <-"
SELECT
	geom.AHdist,
FROM
	hbond_geom_coords AS geom,
	hbonds AS hb,
	hbond_sites AS don, hbond_sites AS acc
WHERE
	geom.struct_id = hb.struct_id AND geom.hbond_id = hb.hbond_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	acc.HBChemType != 'hbacc_PBA' AND
	(don.HBChemType = 'hbdon_AHX' OR don.HBChemType = 'hbdon_HXL');"

	query_sample_sources(sample_sources, sele)

plot_id <- "OHdonor_AHdist_all_acceptor_types"
dens <- estimate_density_1d(
  data = all_geom,
  ids = c("sample_source"),
  variable = "AHdist",
  weight_fun = radial_3d_normalization)

p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	ggtitle("Hydroxyl Donor Hydrogen Bonds A-H Distance\n(normalized for equal weight per unit distance)") +
	scale_y_continuous(
		"Feature Density",
		limits=c(0,2.9), breaks=0:2) +
	scale_x_continuous(
		expression(paste('Acceptor -- Hydrogen Distance (', ring(A), ')')),
		limits=c(1.4,2.7), breaks=c(1.6, 1.9, 2.2, 2.6))

save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id <- "OHdonor_cosBAH_all_acceptor_types"
dens <- estimate_density_1d(
  data = all_geom,
  ids = c("sample_source"),
  variable = "cosBAH")
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(aes(x=acos(x)*180/pi, y=y, colour=sample_source)) +
	geom_indicator(indicator=counts, colour=sample_source, group=sample_source)) +
	ggtitle("Hydroxyl Donor Hydrogen Bonds BAH Angle \n(normalized for equal weight per unit distance)") +
	labs(x=expression(paste('Base -- Acceptor -- Hydrogen (degrees)')),
	     y="Feature Density")


save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id <- "OHdonor_cosAHD_all_acceptor_types"
dens <- estimate_density_1d(
  data = all_geom,
  ids = c("sample_source"),
  variable = "cosAHD")
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=acos(x)*180/pi, y=y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	ggtitle("Hydroxyl Donor Hydrogen Bonds BAH Angle \n(normalized for equal weight per unit distance)") +
	labs(x=expression(paste('Acceptor -- Hydrogen -- Donor (degrees)')),
              y="Feature Density")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id <- "OHdonor_chi_all_acceptor_types"
dens <- estimate_density_1d_wrap(
  data = all_geom,
  ids = c("sample_source"),
  variable = "chi")
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=acos(x)*360/pi, y=y, colour=sample_source)) +
	geom_indicator(aes(colour=sample_source, indicator=counts, group=sample_source)) +
	ggtitle("Hydroxyl Donor Hydrogen Bonds BAH Angle \n(normalized for equal weight per unit distance)") +
	labs(x=expression(paste('Acceptor Base -- Acceptor Torsion (degrees)')),
	     y="Feature Density")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)



sidechain_geom <- all_geom[all_geom$acc_chem_type != "hbacc_PBA",]

plot_id <- "OHdonor_AHdist_sidechain_acceptor_types"
dens <- estimate_density_1d(
  data = sidechain_geom,
  ids = c("sample_source"),
  variable = "AHdist",
  weight_fun = radial_3d_normalization)
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	ggtitle("Hydroxyl Donor Hydrogen Bonds A-H Distance to Sidechain Acceptors\n(normalized for equal weight per unit distance)") +
	scale_y_continuous("Feature Density", limits=c(0,2.9), breaks=0:2) +
	scale_x_continuous(
		expression(paste('Acceptor -- Proton Distance (', ring(A), ')')),
		limits=c(1.4,2.7), breaks=c(1.6, 1.9, 2.2, 2.6))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id <- "OHdonor_cosBAH_sidechain_acceptor_types"
dens <- estimate_density_1d(
  data = sidechain_geom,
  ids = c("sample_source"),
  variable = "cosBAH")
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=acos(x)*360/pi, y=y, colour=sample_source)) +
	geom_indicator(aes(colour=sample_source, indicator=counts, group=sample_source)) +
	ggtitle("Hydroxyl Donor Hydrogen Bonds BAH Angle to Sidechain Acceptors\n(normalized for equal weight per unit distance)") +
	labs(x=expression(paste('Base -- Acceptor -- Hydrogen (degrees)')),
              y="Feature Density")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id <- "OHdonor_cosAHD_sidechain_acceptor_types"
dens <- estimate_density_1d(
  data = sidechain_geom,
  ids = c("sample_source"),
  variable = "cosAHD")
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=acos(x)*360/pi, y=y, colour=sample_source)) +
	geom_indicator(aes(colour=sample_source, indicator=counts, group=sample_source)) +
	ggtitle("Hydroxyl Donor Hydrogen Bonds AHD Angle to Sidechain Acceptors\n(normalized for equal weight per unit distance)") +
	labs(x=expression(paste('Acceptor -- Hydrogen -- Donor (degrees)')),
              y="Feature Density")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id <- "OHdonor_chi_sidechain_acceptor_types"
dens <- estimate_density_1d_logspline(
  data = sidechain_geom,
  ids = c("sample_source"),
  variable = "chi")
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=acos(x)*360/pi, y=y, colour=sample_source)) +
	geom_indicator(aes(colour=sample_source, indicator=counts, group=sample_source)) +
	ggtitle("Hydroxyl Donor Hydrogen Bonds chi Torsion Angle to Sidechain Acceptors\n(normalized for equal weight per unit distance)") +
	labs(x=expression(paste('Acceptor Base -- Acceptor Torsion (degrees)')),
	     y="log(FeatureDensity + 1)")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
