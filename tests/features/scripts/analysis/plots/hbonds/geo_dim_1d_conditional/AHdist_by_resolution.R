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
id = "AHdist_by_resolution",
author = "Matthew O'Meara",
brief_description = "This measures the H-Bond A-H distance conditional on the resolution. Note that currently there is no features reporter for resolution so it must be included after the fact.",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

source("scripts/analysis/plots/hbonds/hbond_geo_dim_scales.R")

sele <-"
SELECT
  geom.AHdist,
  acc_site.HBChemType AS acc_chem_type,
  don_site.HBChemType AS don_chem_type,
  resolutions.resolution
FROM
  hbond_geom_coords AS geom,
  hbonds AS hbond,
  hbond_sites AS don_site,
  hbond_sites AS acc_site,
  resolutions
WHERE
  hbond.struct_id = resolutions.struct_id AND
  hbond.struct_id = geom.struct_id AND
  hbond.hbond_id =  geom.hbond_id AND
  hbond.struct_id = don_site.struct_id AND
  hbond.don_id = don_site.site_id AND
  hbond.struct_id = acc_site.struct_id AND
  hbond.acc_id = acc_site.site_id;"
f <- query_sample_sources(sample_sources, sele)

f$don_chem_type_name <- don_chem_type_name_linear(f$don_chem_type))
f$acc_chem_type_name <- acc_chem_type_name_linear(f$acc_chem_type))

f$resolution_q <-
  cut(f$resolution,
  quantile(f$resolution, seq(0,1,by=.25)))


dens <- estimate_density_1d(
  f, c("resolution_q", "don_chem_type_name", "acc_chem_type_name"),
  "AHdist", weight_fun = radial_3d_normalization)


plot_id <- "AHdist_by_resolution"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=resolution_q)) +
	geom_indicator(aes(indicator=counts, colour=resolution_q, group=resolution_q)) +
	facet_grid(don_chem_type_name ~ acc_chem_type_name) +
	ggtitle("Hydrogen Bonds A-H Distance by Chemical Type by Resolution Quantiles\nnormalized for equal weight per unit distance") +
	scale_y_continuous("FeatureDensity", limits=c(0,2.9), breaks=0:2) +
	scale_x_continuous(expression(paste('Acceptor -- Proton Distance (', ring(A), ')')), limits=c(1.4,2.7), breaks=c(1.6, 1.9, 2.2, 2.6))

save_plots(self, plot_id, sample_sources, output_dir, output_formats)



})) # end FeaturesAnalysis
