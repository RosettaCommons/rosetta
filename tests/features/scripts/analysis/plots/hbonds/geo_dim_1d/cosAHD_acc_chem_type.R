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
id = "cosAHD_acc_chem_type",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	geom.cosAHD,
	acc_site.HBChemType AS acc_chem_type
FROM
	hbond_geom_coords AS geom,
	hbonds AS hb,
	hbond_sites AS don_site,
	hbond_sites AS acc_site,
	hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb
WHERE
	hb.struct_id = geom.struct_id AND
	hb.hbond_id =  geom.hbond_id AND
	hb.struct_id = don_site.struct_id AND
	hb.don_id = don_site.site_id AND
	hb.struct_id = acc_site.struct_id AND
	hb.acc_id = acc_site.site_id AND
	don_pdb.struct_id = hb.struct_id AND
	don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND
	acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30 AND
	abs( don_site.resNum - acc_site.resNum ) > 5;"

f <- query_sample_sources(sample_sources, sele)

# Order the plots better and give more descriptive labels
f$acc_chem_type <- factor(f$acc_chem_type,
	levels =
		c("hbacc_CXA", "hbacc_AHX", "hbacc_IMD",
			"hbacc_CXL", "hbacc_HXL", "hbacc_IME",
			"hbacc_PBA"),
	labels =
		c("aCXA: n,q", "aAHX: y",   "aIMD: h",
			"aCXL: d,e", "aHXL: s,t", "aIME: h",
			"aPBA: bb"))

#coAHD goes from 0 to 1, where 1 is linear
#since there is significant density at 1,
#to accurately model a discontinuity, reflect
#the data across the right boundary, in computing the density esitmation
dens <- estimate_density_1d_reflect_boundary(
	data=f,
	ids = c("sample_source", "acc_chem_type"),
	variable = "cosAHD",
	reflect_right=TRUE,
	right_boundary=1,
	adjust=.3)

plot_id = "cosAHD_acc_chem_type"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=180-acos(x)*180/pi, y=y, colour=sample_source)) +
	geom_indicator(aes(colour=sample_source, indicator=counts, group=sample_source)) +
	facet_wrap( ~ acc_chem_type) +
	ggtitle("Hydrogen Bonds AHD Angle by Acceptor Chemical Type; SeqSep > 5; BFact < 30\n(normalized for equal volume per unit distance)") +
	scale_y_continuous("FeatureDensity", limits=c(0,20), breaks=c(0,5,10,15)) +
	scale_x_continuous("Acceptor -- Hydrogen -- Donor (degrees)", trans="reverse")+
	theme(legend.position=c(.58,.35)) +
	theme(legend.justification=c("left", "top"))

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
