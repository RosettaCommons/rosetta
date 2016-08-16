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
id = "cosAHD_beta_regression",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

source("scripts/analysis/plots/hbonds/hbond_geo_dim_scales.R")

sele <-"
SELECT
  geom.cosAHD,
  acc_site.HBChemType AS acc_chem_type,
  don_site.HBChemType AS don_chem_type
FROM
  hbond_geom_coords AS geom,
  hbonds AS hbond,
  hbond_sites AS don_site,
  hbond_sites AS acc_site
WHERE
  hbond.struct_id = geom.struct_id AND hbond.hbond_id = geom.hbond_id AND
  hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
  hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id;"

f <- query_sample_sources(sample_sources, sele)


dens <- estimate_density_1d_reflect_boundary(
  f, c("sample_source"), "cosAHD", reflect_right=T, right_boundary=1)



plot_parts <- list(
  theme_bw(),
  scale_x_cosAHD,
  scale_y_continuous("FeatureDensity"))

plot_id <- "hbond_cosAHD"
ggplot(dens) + plot_parts +
  geom_line(aes(x=x, y=y, color=sample_source)) +
  geom_indicator(aes(indicator=counts, color=sample_source, group=sample_source)) +
  ggtitle("Hydrogen Bonds cosAHD") +
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
