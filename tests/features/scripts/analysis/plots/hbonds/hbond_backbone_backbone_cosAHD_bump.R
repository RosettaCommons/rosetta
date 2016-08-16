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
id = "hbond_backbone_backbone_cosAHD_bump",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <- "
SELECT
	geom.cosAHD AS cosAHD,
	geom.AHdist AS AHdist,
	geom.cosBAH AS cosBAH,
	don_site.resType AS don_res_type,
	acc_site.resType AS acc_res_type
FROM
	hbonds AS hbond,
	hbond_geom_coords AS geom,
	hbond_sites AS don_site,
	hbond_sites AS acc_site
WHERE
	hbond.struct_id = geom.struct_id AND
	hbond.hbond_id =  geom.hbond_id AND
	hbond.struct_id = don_site.struct_id AND
	hbond.don_id = don_site.site_id AND
	hbond.struct_id = acc_site.struct_id AND
	hbond.acc_id = acc_site.site_id AND
	don_site.HBChemType = 'hbdon_PBA' AND
	acc_site.HBChemType = 'hbacc_PBA' AND
	cosAHD > .63 AND
	cosAHD < .77;"

f <- query_sample_sources(sample_sources, sele)

plot_parts <- list(
	theme_bw(),
	geom_line(aes(x=acos(x)*360/pi, y=log(y+1), colour=sample_source)),
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)),
	geom_vline(xintercept = 91.05),
	scale_x_continuous("Acceptor -- Hydrogen -- Donor (Degrees)"),
	scale_y_continuous("log(FeatureDensity + 1)"))

plot_id <- "hbond_backbone_backbon_cosAHD_bump"
dens <- estimate_density_1d(
	f, c("sample_source"), "cosAHD", adjust=.1, histogram=TRUE)
ggplot(data=dens) + plot_parts +
	ggtitle("Backbone Backbone Hydrogen Bond AHD Angle\nNormalized for equal Weight per Unit Angle") +
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id <- "hbond_backbone_backbon_cosAHD_bump_by_don_res_type"
dens <- estimate_density_1d(
	f, c("sample_source", "don_res_type"), "cosAHD", adjust=.1)
ggplot(data=dens) + plot_parts +
	facet_wrap( ~ don_res_type) +
	ggtitle("Backbone Backbone Hydrogen Bond AHD Angle by donor and acceptor residue types\nnormalized for equal weight per unit distance")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id <- "hbond_backbone_backbon_cosAHD_bump_by_acc_res_type"
dens <- estimate_density_1d(
	f, c("sample_source", "acc_res_type"), "cosAHD", adjust=.1)
ggplot(data=dens) + plot_parts +
	facet_wrap(~acc_res_type) +
	ggtitle("Backbone Backbone Hydrogen Bond AHD Angle by donor and acceptor residue types\nnormalized for equal weight per unit distance")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


#dens <- estimate_density_1d(
#	f, c("sample_source"), "AHdist",
#	weight_fun=radial_3d_normalization,
#	adjust=.1)
#
#plot_id <- "hbond_backbone_backbon_cosAHD_bump_AHdist"
#ggplot(data=dens) + plot_parts +
#	facet_grid(don_res_type ~ acc_res_type) +
#	ggtitle("Backbone Backbone Hydrogen Bond AHdist by donor and acceptor residue types (restricted to cosAHD bump)\nnormalized for equal weight per unit distance")
#save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
