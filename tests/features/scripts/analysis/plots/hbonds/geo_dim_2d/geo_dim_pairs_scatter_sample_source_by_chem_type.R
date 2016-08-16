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
id = "geo_dim_pairs_scatter_chem_type_by_sample_source",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

source("scripts/analysis/plots/hbonds/hbond_geo_dim_scales.R")

sele <-"
SELECT
	geom.AHdist, geom.cosBAH, geom.cosAHD, geom.chi,
	acc.HBChemType AS acc_chem_type, don.HBChemType AS don_chem_type
FROM
	hbonds AS hb cross join
	hbond_geom_coords AS geom cross join
	hbond_sites AS don cross join
	hbond_sites AS acc cross join
	hbond_sites_pdb AS don_pdb cross join
	hbond_sites_pdb AS acc_pdb cross join
	residue_secondary_structure AS don_ss cross join
	dssp_codes AS don_ss_code cross join
	residue_secondary_structure AS acc_ss cross join
	dssp_codes AS acc_ss_code
WHERE
	geom.struct_id = hb.struct_id AND geom.hbond_id = hb.hbond_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30 AND
	geom.AHdist < 3.2 AND
	acc_ss.struct_id = acc.struct_id AND acc_ss.resNum = acc.resNum AND
	acc_ss_code.code = acc_ss.dssp AND
	don_ss.struct_id = don.struct_id AND don_ss.resNum = don.resNum AND
	don_ss_code.code = don_ss.dssp AND
	(don_ss.dssp != 'H' OR acc_ss.dssp != 'H' OR don.resNum - acc.resNum != 3);"

f <- query_sample_sources(sample_sources, sele)

f <- transform(f,
	don_chem_type_name = don_chem_type_name_linear(don_chem_type),
	acc_chem_type_name = acc_chem_type_name_linear(acc_chem_type),
	BAH_CHI_x = 2*sin(acos(cosBAH)/2)*cos(chi),
	BAH_CHI_y = 2*sin(acos(cosBAH)/2)*sin(chi))

f <- na.omit(f, method="r")

plot_parts <- list(
	theme_bw(),
	facet_wrap(~sample_source, ncol=2),
	geom_point(size=.8),
	stat_density2d(size=.2),
	coord_equal(ratio=1))


set_plot_title <- function(xdim, ydim, don_chem_type, acc_chem_type){
	ggtitle(
		paste(
			"Hydrogen Bonds ", xdim, " vs ", ydim, "  ",
			"don: ", don_chem_type, " acc: ", acc_chem_type, sep=""))
}

set_plot_id <- function(xdim, ydim, don_chem_type, acc_chem_type){
	paste("hbond_geo_dim_pairs_scatter", xdim, ydim, don_chem_type, acc_chem_type, sep="_")
}



plot_id <- "hbond_geom_dim_pairs_scatter_cosBAH_AHdist3"
ggplot(data=f, aes(x=cosBAH, y=AHdist^3)) +
	theme_bw() +
	geom_point(size=.4) +
	stat_density2d(size=.2) +
	scale_x_cosBAH + scale_y_AHdist3 +
	ggtitle("HBonds cosBAH vs AHdist^3") +
	facet_grid( ~ sample_source)
save_plots(self, plot_id, sample_sources, output_dir, alt_output_formats)

plot_id <- "hbond_geom_dim_pairs_scatter_cosAHD_AHdist3"
ggplot(data=f, aes(x=cosAHD, y=AHdist^3)) +
	theme_bw() +
	geom_point(size=.4) +
	stat_density2d(size=.2) +
	scale_x_cosAHD + scale_y_AHdist3 +
	ggtitle("HBonds cosAHD vs AHdist^3") +
	facet_grid( ~ sample_source)
save_plots(self, plot_id, sample_sources, output_dir, alt_output_formats)

plot_id <- "hbond_geom_dim_pairs_scatter_cosAHD_cosBAH"
ggplot(data=f, aes(x=cosAHD, y=cosBAH)) +
	theme_bw() +
	geom_point(size=.4) +
	stat_density2d(size=.2) +
	scale_x_cosAHD + scale_y_cosBAH +
	ggtitle("HBonds cosAHD vs cosBAH") +
	facet_wrap( ~ sample_source, ncol=ceiling(sqrt(nrow(sample_sources)))) +
	coord_equal(ratio=1)
alt_output_formats <- transform(output_formats, width=height)
save_plots(self, plot_id, sample_sources, output_dir, alt_output_formats)


d_ply(
	f[f$don_chem_type == "hbdon_PBA" & f$acc_chem_type == "hbacc_PBA",],
	.(seq_sep), function(sub_f){
	seq_sep <- as.character(sub_f$seq_sep[1])
	plot_id <- paste("hbond_geom_dim_pairs_scatter_cosAHD_AHdist3_by_don_chem_type", seq_sep, sep="_")
	narrow_output_formats <- transform(output_formats, width=height/1.8)
	ggplot(data=sub_f, aes(x=cosBAH, y=AHdist^3)) +
		theme_bw() +
		geom_point(size=.4) +
		stat_density2d(size=.2) +
		scale_x_cosBAH + scale_y_AHdist3 +
		ggtitle(paste("HBonds cosBAH vs AHdist^3, seq_sep:", seq_sep)) +
		facet_grid( ~ sample_source)
	save_plots(self, plot_id, sample_sources, output_dir, narrow_output_formats)
})

plot_each_chem_type <- function(sub_f){
	don_chem_type <- as.character(sub_f$don_chem_type_name[1])
	acc_chem_type <- as.character(sub_f$acc_chem_type_name[1])
	print(don_chem_type)
	print(acc_chem_type)

	plot_id <- set_plot_id("cosBAH", "AHdist", don_chem_type, acc_chem_type)
	ggplot(data=sub_f, aes(x=cosBAH, y=AHdist)) + plot_parts +
		scale_x_cosBAH + scale_y_AHdist +
		set_plot_title("cosBAH", "AHdist", don_chem_type, acc_chem_type)
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

#	plot_id <- set_plot_id("cosAHD", "AHdist", don_chem_type, acc_chem_type)
#	p <- ggplot(data=sub_f, aes(x=cosAHD, y=AHdist)) + plot_parts +
#		scale_x_cosAHD + scale_y_AHdist +
#		set_plot_title("cosAHD", "AHdist", don_chem_type, acc_chem_type)
#	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
#
#	plot_id <- set_plot_id("CHI", "AHdist", don_chem_type, acc_chem_type)
#	p <- ggplot(data=sub_f, aes(x=chi, y=AHdist)) + plot_parts +
#		scale_x_chi + scale_y_AHdist +
#		set_plot_title("CHI", "AHdist", don_chem_type, acc_chem_type)
#	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
#
#	plot_id <- set_plot_id("cosAHD", "cosBAH", don_chem_type, acc_chem_type)
#	p <- ggplot(data=sub_f, aes(x=cosAHD, y=cosBAH)) + plot_parts +
#		scale_x_cosAHD + scale_y_cosBAH +
#		set_plot_title("cosAHD", "cosBAH", don_chem_type, acc_chem_type)
#	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
#
##	plot_id <- set_plot_id("CHI", "cosBAH", don_chem_type, acc_chem_type)
##	p <- ggplot(data=sub_f, aes(x=chi, y=cosBAH)) + plot_parts +
##		scale_x_chi + scale_y_cosBAH +
##		set_plot_title("CHI", "cosBAH", don_chem_type, acc_chem_type)
##	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
#
#	plot_id <- set_plot_id("CHI", "cosAHD", don_chem_type, acc_chem_type)
#	p <- ggplot(data=sub_f, aes(x=chi, y=cosAHD)) + plot_parts +
#		scale_x_chi + scale_y_cosAHD +
#		set_plot_title("CHI", "cosAHD", don_chem_type, acc_chem_type)
#	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
#
##	plot_id <- set_plot_id("CHI", "BAH", don_chem_type, acc_chem_type)
##	p <- ggplot(data=sub_f, aes(x=BAH_CHI_x, y=BAH_CHI_y)) + plot_parts +
##		set_plot_title("CHI", "BAH", don_chem_type, acc_chem_type)
##	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

}

runtime <- system.time(d_ply(f, .(don_chem_type, acc_chem_type), .fun=plot_each_chem_type))
print(paste("Plot Generation Time: ", runtime, sep=""))



})) # end FeaturesAnalysis
