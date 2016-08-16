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
id = "chi_BAH_eqpoldens_scatter_lrbbbbhbonds",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	geom.cosBAH, geom.chi,
	acc.HBChemType AS acc_chem_type, don.HBChemType AS don_chem_type
FROM
	hbond_geom_coords AS geom,
	hbonds AS hbond,
	hbond_sites AS don, hbond_sites AS acc
WHERE
  don.HBChemType == 'hbdon_PBA' AND acc.HBChemType == 'hbacc_PBA' AND
	hbond.struct_id = geom.struct_id AND hbond.hbond_id = geom.hbond_id AND
	hbond.struct_id = don.struct_id AND hbond.don_id = don.site_id AND
	hbond.struct_id = acc.struct_id AND hbond.acc_id = acc.site_id AND
	ABS(don.resNum - acc.resNum) > 5;";
f <- query_sample_sources(sample_sources, sele)

#equal area projection
f <- transform(f,
	capx = 2*sin(acos(cosBAH)/2)*cos(chi),
	capy = 2*sin(acos(cosBAH)/2)*sin(chi))

capx_limits <- c(-1.5,1.5)
capy_limits <- capx_limits


f_first <- f[ f$sample_source == levels(sample_sources$sample_source)[1], ]
f_second <- f[ f$sample_source == levels(sample_sources$sample_source)[2], ]

f_first <- ddply(
	f_first, .(sample_source), transform, counts = length(sample_source))

plot_id = "chi_BAH_eqpoldens_and_scatter_lrbb"
ggplot(data=f_first) + theme_bw() +
	theme(panel.background=element_rect(fill="#00007F", colour="#00007F")) +
	stat_density2d(
		aes(x=capx,y=capy, fill=..density..), geom="tile", contour=FALSE ) +
	polar_equal_area_grids_bw() +
	geom_indicator(aes(indicator=counts), color="white") +
	geom_point( data=f_second,aes(x=capx,y=capy),colour="white",size=2) +
	ggtitle(
		paste("Hydrogen Bonds chi vs BAH Angles with Sequence Separation > 5\n",
		"Backbone/Backbone Hydrogen Bonds, Equal Coordinate Projection\n",
		"Reference (density) vs Test (white circles)", sep="")) +
	scale_x_continuous(
		'2*sin(BAH/2) * cos(CHI)', limits=capx_limits, breaks=c(-1, 0, 1)) +
	scale_y_continuous(
		'2*sin(BAH/2) * sin(CHI)', limits=capy_limits, breaks=c(-1, 0, 1)) +
	coord_fixed(ratio = 1) +
	scale_fill_gradientn('Density', colours=jet.colors(10))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
