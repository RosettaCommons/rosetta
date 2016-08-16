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
id = "ADdist_backbone_backbone_by_squence_separation",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

source("scripts/analysis/plots/hbonds/hbond_geo_dim_scales.R")

sele <-"
SELECT
	don_atoms.base_x AS dx, don_atoms.base_y AS dy, don_atoms.base_z AS dz,
	acc_atoms.atm_x  AS ax, acc_atoms.atm_y  AS ay, acc_atoms.atm_z  AS az,
	CASE don.resNum - acc.resNum
		WHEN -1 THEN '-1' WHEN -2 THEN '-2' WHEN -3 THEN '-3' WHEN -4 THEN '-4'
		WHEN 2 THEN '2' WHEN 3 THEN '3' WHEN 4 THEN '4' WHEN 5 THEN '5'
		ELSE 'long' END AS seq_sep
FROM
	hbonds AS hb,
	hbond_sites AS don,	hbond_sites AS acc,
	hbond_site_atoms AS don_atoms, hbond_site_atoms AS acc_atoms
WHERE
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	acc.HBChemType == 'hbacc_PBA' AND
	don.HBChemType == 'hbdon_PBA' AND
	don_atoms.struct_id = hb.struct_id AND don_atoms.site_id = hb.don_id AND
	acc_atoms.struct_id = hb.struct_id AND acc_atoms.site_id = hb.acc_id;"
f <- query_sample_sources(sample_sources, sele)

# A-D distance is not stored directly in the features database,
# however it can be computed from the coordinates of the hydrogen
# bonding atoms.
f <- transform(f,
	ADdist = vector_distance(cbind(dx, dy, dz), cbind(ax, ay, az)))


# This shouldn't happend, but if it does get rid of them
f <- f[f$ADdist < 5,]


f$seq_sep <- factor(f$seq_sep,
	levels = c("-4", "-3", "-2", "-1", "2", "3", "4", "5", "long"),
	labels = c("-4", "-3", "-2", "-1", "2", "3", "4", "5", "long"))


dens <- estimate_density_1d(
	f, c("sample_source", "seq_sep"),
	"ADdist", weight_fun = radial_3d_normalization)

plot_id <- "hbond_ADdist_backbone_backbone_by_sequence_separation"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	facet_wrap( ~ seq_sep ) +
	ggtitle("Backbone-Backbone HBonds A-D Distance by Sequence Separation\n(DonRes - AccRes) normalized for equal weight per unit distance") +
	scale_x_ADdist +
	scale_y_continuous("FeatureDensity", limits=c(0,6), breaks=c(1,3,5))
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
