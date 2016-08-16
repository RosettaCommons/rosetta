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
id = "chi_BAH_eq_polar_density_srPBAtoPBA",
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
	hbond_sites AS don,	hbond_sites AS acc
WHERE
  don.HBChemType == 'hbdon_PBA' AND acc.HBChemType == 'hbacc_PBA' AND
	hbond.struct_id = geom.struct_id AND hbond.hbond_id = geom.hbond_id AND
	hbond.struct_id = don.struct_id AND hbond.don_id = don.site_id AND
	hbond.struct_id = acc.struct_id AND hbond.acc_id = acc.site_id AND
	ABS(don.resNum - acc.resNum) <= 5;";
f <- query_sample_sources(sample_sources, sele)

#equal area projection
f <- transform(f,
	capx = 2*sin(acos(cosBAH)/2)*cos(chi),
	capy = 2*sin(acos(cosBAH)/2)*sin(chi))

f <- ddply(
	f, .(sample_source, don_chem_type, acc_chem_type),
	transform, counts = length(sample_source))

capx_limits <- c(-1.5,1.5)
capy_limits <- capx_limits

plot_id = "chi_BAH_eq_polar_density_sr_bbbb"
d_ply(f, .(sample_source), function(sub_f){
	ss_id <- sub_f$sample_source[1]
	ss <- sample_sources[sample_sources$sample_source == ss_id, ]

	ggplot(data=sub_f) + theme_bw() +
		stat_density2d(
			aes(x=capx,y=capy, fill=..density..), geom="tile", contour=FALSE) +
		polar_equal_area_grids_bw(bgcolor="#00007F") +
		geom_indicator(aes(indicator=counts), color="white") +
		facet_grid(acc_chem_type ~ don_chem_type) +
		ggtitle(
			paste("Hydrogen Bonds chi vs BAH Angles with Sequence Separation <= 5\n",
			"Backbone/Backbone Hydrogen Bonds\n",
			"Equal Coordinate Projection   Sample Source: ", ss_id, sep="")) +
		scale_x_continuous(
			'2*sin(BAH/2) * cos(CHI)', limits=capx_limits, breaks=c(-1, 0, 1)) +
		scale_y_continuous(
			'2*sin(BAH/2) * sin(CHI)', limits=capy_limits, breaks=c(-1, 0, 1)) +
		coord_fixed(ratio = 1) +
		scale_fill_gradientn('Density', colours=jet.colors(10))
	save_plots(self, plot_id, ss, output_dir, output_formats)
})


})) # end FeaturesAnalysis
