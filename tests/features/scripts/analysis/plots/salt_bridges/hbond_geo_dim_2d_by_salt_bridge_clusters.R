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
id = "geo_dim_pairs_scatter_chem_type",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

source("scripts/analysis/plots/hbonds/hbond_geo_dim_scales.R")
source("scripts/analysis/plots/salt_bridges/salt_bridge_clusters.R")

sele <-"
SELECT
	geom.AHdist, geom.cosBAH, geom.cosAHD, geom.chi,
	acc.HBChemType AS acc_chem_type, don.HBChemType AS don_chem_type
FROM
	hbond_geom_coords AS geom,
	hbonds AS hb,
	hbond_sites AS don, hbond_sites AS acc,
	hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb
WHERE
	geom.struct_id = hb.struct_id AND geom.hbond_id = hb.hbond_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30;"
f <- query_sample_sources(sample_sources, sele)

f$don_chem_type_name <- don_chem_type_name_linear(f$don_chem_type)
f$acc_chem_type_name <- acc_chem_type_name_linear(f$acc_chem_type)
f <- na.omit(f, method="r")

plot_parts <- list(
	theme_bw(),
	facet_grid(don_chem_type_name ~ acc_chem_type_name),
	geom_point(size=.4),
	stat_density2d(size=.2))

set_plot_title <- function(xdim, ydim, ss_id){
	ggtitle(
		paste("Hydrogen Bonds ", xdim, " vs ", ydim, "  ss_id: ", ss_id, sep=""))
}

set_plot_id <- function(xdim, ydim, ss_id){
	paste("hbond_geo_dim_pairs_scatter", xdim, ydim, "chem_type", ss_id,  sep="_")
}

plot_each_ss <- function(sub_f){
	ss_id <- sub_f$sample_source[1]
	ss <- sample_sources[sample_sources$sample_source == ss_id,]

	print(ss_id)
	print(ss)

	plot_id <- set_plot_id("cosBAH", "AHdist", ss_id)
	ggplot(data=sub_f, aes(x=cosBAH, y=AHdist)) + plot_parts +
		scale_x_cosBAH + scale_y_AHdist + set_plot_title("cosBAH", "AHdist", ss_id)
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

	plot_id <- set_plot_id("cosAHD", "AHdist", ss_id)
	p <- ggplot(data=sub_f, aes(x=cosAHD, y=AHdist)) + plot_parts +
		scale_x_cosAHD + scale_y_AHdist + set_plot_title("cosAHD", "AHdist", ss_id)
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

	plot_id <- set_plot_id("CHI", "AHdist", ss_id)
	p <- ggplot(data=sub_f, aes(x=chi, y=AHdist)) + plot_parts +
		scale_x_chi + scale_y_AHdist + set_plot_title("CHI", "AHdist", ss_id)
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

	plot_id <- set_plot_id("cosAHD", "cosBAH", ss_id)
	p <- ggplot(data=sub_f, aes(x=cosAHD, y=cosBAH)) + plot_parts +
		scale_x_cosAHD + scale_y_cosBAH + set_plot_title("cosAHD", "cosBAH", ss_id)
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

	plot_id <- set_plot_id("CHI", "cosAHD", ss_id)
	p <- ggplot(data=sub_f, aes(x=chi, y=cosBAH)) + plot_parts +
		scale_x_chi + scale_y_cosBAH + set_plot_title("CHI", "cosBAH", ss_id)
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

	plot_id <- set_plot_id("CHI", "cosAHD", ss_id)
	p <- ggplot(data=sub_f, aes(x=chi, y=cosAHD)) + plot_parts +
		scale_x_chi + scale_y_cosAHD + set_plot_title("CHI", "cosAHD", ss_id)
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
}

runtime <- system.time(d_ply(f, .(sample_source), .fun=plot_each_ss))
print(paste("Plot Generation Time: ", runtime, sep=""))



plot_parts <- list(
	theme_bw(),
	facet_wrap(~sample_source),
	geom_point(size=.4),
	stat_density2d(size=.2))

set_plot_title <- function(xdim, ydim, don_chem_type, acc_chem_type){
	ggtitle(
		paste(
			"Hydrogen Bonds ", xdim, " vs ", ydim, "  ",
			"don: ", don_chem_type, " acc: ", acc_chem_type, sep=""))
}

set_plot_id <- function(xdim, ydim, don_chem_type, acc_chem_type){
	paste("hbond_geo_dim_pairs_scatter", xdim, ydim, don_chem_type, acc_chem_type, sep="_")
}


plot_each_chem_type <- function(sub_f){
	don_chem_type <- as.character(sub_f$don_chem_type_name[1])
	acc_chem_type <- as.character(sub_f$acc_chem_type_name[1])

	plot_id <- set_plot_id("cosBAH", "AHdist", don_chem_type, acc_chem_type)
	ggplot(data=sub_f, aes(x=cosBAH, y=AHdist)) + plot_parts +
		scale_x_cosBAH + scale_y_AHdist +
		set_plot_title("cosBAH", "AHdist", don_chem_type, acc_chem_type)
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

	plot_id <- set_plot_id("cosAHD", "AHdist", don_chem_type, acc_chem_type)
	p <- ggplot(data=sub_f, aes(x=cosAHD, y=AHdist)) + plot_parts +
		scale_x_cosAHD + scale_y_AHdist +
		set_plot_title("cosAHD", "AHdist", don_chem_type, acc_chem_type)
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

	plot_id <- set_plot_id("CHI", "AHdist", don_chem_type, acc_chem_type)
	p <- ggplot(data=sub_f, aes(x=chi, y=AHdist)) + plot_parts +
		scale_x_chi + scale_y_AHdist +
		set_plot_title("CHI", "AHdist", don_chem_type, acc_chem_type)
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

	plot_id <- set_plot_id("cosAHD", "cosBAH", don_chem_type, acc_chem_type)
	p <- ggplot(data=sub_f, aes(x=cosAHD, y=cosBAH)) + plot_parts +
		scale_x_cosAHD + scale_y_cosBAH +
		set_plot_title("cosAHD", "cosBAH", don_chem_type, acc_chem_type)
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

	plot_id <- set_plot_id("CHI", "cosBAH", don_chem_type, acc_chem_type)
	p <- ggplot(data=sub_f, aes(x=chi, y=cosBAH)) + plot_parts +
		scale_x_chi + scale_y_cosBAH +
		set_plot_title("CHI", "cosBAH", don_chem_type, acc_chem_type)
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

	plot_id <- set_plot_id("CHI", "cosAHD", don_chem_type, acc_chem_type)
	p <- ggplot(data=sub_f, aes(x=chi, y=cosAHD)) + plot_parts +
		scale_x_chi + scale_y_cosAHD +
		set_plot_title("CHI", "cosAHD", don_chem_type, acc_chem_type)
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

}

runtime <- system.time(d_ply(f, .(don_chem_type, acc_chem_type), .fun=plot_each_chem_type))
print(paste("Plot Generation Time: ", runtime, sep=""))



})) # end FeaturesAnalysis
