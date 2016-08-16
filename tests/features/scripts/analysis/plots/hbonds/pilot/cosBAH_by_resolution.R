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
id = "cosBAH_by_resolution",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

source("scripts/analysis/plots/hbonds/hbond_geo_dim_scales.R")

sele <-"
SELECT
  geom.cosBAH,
  acc.HBChemType AS acc_chem_type, don.HBChemType AS don_chem_type,
  resolution.resolution
FROM
  hbond_geom_coords AS geom,
  hbonds AS hb,
  hbond_sites AS don, hbond_sites AS acc,
  resolutions AS resolution
WHERE
	geom.struct_id = hb.struct_id AND geom.hbond_id = hb.hbond_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	resolution.struct_id = hb.struct_id;"

f <- query_sample_sources(sample_sources, sele)

f$don_chem_type_name <- don_chem_type_name_linear(f$don_chem_type)
f$acc_chem_type_name <- acc_chem_type_name_linear(f$acc_chem_type)
f <- na.omit(f, method="r")

f$resolution_q <-
  cut(f$resolution, quantile(f$resolution, seq(0,1,by=.25)))

dens <- estimate_density_1d(
  data = f,
  ids = c("resolution_q"),
  variable = "cosBAH")

ddply(f, .(sample_source), function(sub_f){
	ss_id <- f$sample_source[1]
	ss <- sample_sources[sample_sources$sample_source == ss_id,]

	plot_id = "hbonds_cosBAH_chem_type_by_resolution_qantile"
	ggplot(data=dens) + theme_bw() +
		geom_line(aes(x=180-acos(x)*180/pi, y=y, colour=resolution_q)) +
		geom_indicator(aes(colour=resolution_q, indicator=counts, group=resolution_q)) +
		facet_grid(don_chem_type ~ acc_chem_type) +
		ggtitle("HBonds BAH Angle by Resolution Quantile\n(normalized for equal volume per unit distance)") +
		scale_x_continuous(paste('Base -- Acceptor -- Hydrogen (degrees)')) +
		scale_y_continuous("FeatureDensity", limits=c(0,2.9), breaks=0:2)

	save_plots(self, plot_id, ss, output_dir, output_formats)
})

})) # end FeaturesAnalysis
