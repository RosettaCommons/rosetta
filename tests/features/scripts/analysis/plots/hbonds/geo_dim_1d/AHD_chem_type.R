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
id = "AHD_chem_type",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

source("scripts/analysis/plots/hbonds/hbond_geo_dim_scales.R")

sele <-"
SELECT
  geom.cosAHD,
	don.HBChemType AS don_chem_type, acc.HBChemType AS acc_chem_type
FROM
	hbonds AS hb,
	hbond_geom_coords AS geom,
	hbond_sites AS don, hbond_sites AS acc,
	hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb
WHERE
	geom.struct_id = hb.struct_id AND geom.hbond_id = hb.hbond_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30 AND
	abs(don.resNum - acc.resNum) > 5;"

f <- query_sample_sources(sample_sources, sele)

f$don_chem_type_name <- don_chem_type_name_linear(f$don_chem_type)
f$acc_chem_type_name <- acc_chem_type_name_linear(f$acc_chem_type)
f <- na.omit(f, method="r")


f$AHD <- acos(f$cosAHD)

#AHD goes from 0 to pi/2, where 0 is linear. Since there is significant
#density at 0, to accurately model a discontinuity, reflect the data
#across the left boundary, in computing the density esitmation
dens <- estimate_density_1d_reflect_boundary(
	data=f,
	ids = c("sample_source", "acc_chem_type_name", "don_chem_type_name"),
	variable = "AHD",
	reflect_left=TRUE,
	left_boundary=0,
	weight_fun=conical_3d_normalization,
	adjust=.5)

plot_id = "hbond_AHD_chem_type"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=180-x*180/pi, y=y, colour=sample_source)) +
	geom_indicator(aes(colour=sample_source, indicator=counts, group=sample_source)) +
	facet_grid(don_chem_type_name ~ acc_chem_type_name) +
	ggtitle("HBond AHD Angle by Chemical Type, SeqSep > 5, B-Fact < 30\n(normalized for equal volume per unit distance)") +
	scale_y_continuous("Feature Density") +
	scale_x_continuous(
		"Acceptor -- Hydrogen -- Donor (degrees)", trans="reverse",
		limits=c(180, 120), breaks=c(180, 160, 140))

if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
