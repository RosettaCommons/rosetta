# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

check_setup()
feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "hbond_geo_dim_2d_dARG_aCXL_by_salt_bridge_clusters",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

source("scripts/analysis/plots/hbonds/hbond_geo_dim_scales.R")
source("scripts/analysis/plots/salt_bridges/salt_bridge_clusters.R")

sele <-"
CREATE INDEX IF NOT EXISTS salt_bridges_struct_id_don_resNum_acc_id ON
	salt_bridges ( struct_id, don_resNum, acc_id );

SELECT
	geom.AHdist, geom.cosBAH, geom.cosAHD, geom.chi,
	acc.HBChemType AS acc_chem_type, don.HBChemType AS don_chem_type,
	sb.psi*180/3.141593 AS psi_degrees, sb.rho
FROM
	hbond_geom_coords AS geom,
	hbonds AS hb,
	hbond_sites AS don, hbond_sites AS acc,
	hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb,
	salt_bridges AS sb
WHERE
	geom.struct_id = hb.struct_id AND geom.hbond_id = hb.hbond_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30 AND
	sb.struct_id = don.struct_id AND sb.don_resNum = don.resNum AND
	sb.acc_id = acc.site_id AND
	(don.HBChemType == 'hbdon_GDE' OR don.HBChemType == 'hbdon_GDH') AND
	acc.HBChemType == 'hbacc_CXL';"
f <- query_sample_sources(sample_sources, sele)

#equal area projection
f <- transform(f,
	BAH_CHI_x = 2*sin(acos(cosBAH)/2)*cos(chi),
	BAH_CHI_y = 2*sin(acos(cosBAH)/2)*sin(chi))


f$don_chem_type_name <- don_chem_type_name_linear(f$don_chem_type)
f$acc_chem_type_name <- acc_chem_type_name_linear(f$acc_chem_type)
f <- na.omit(f, method="r")

f$cluster <- add_salt_bridge_clusters(f)
f <- na.omit(f, method="r")

f$is_bifurcated_a <- as.character(f$cluster) == "bifurcated_a"

plot_parts <- list(
	theme_bw(),
	facet_wrap(~sample_source, nrow=2),
	#geom_point(aes(colour=is_bifurcated_a), size=.8),
	stat_density2d(aes(colour=is_bifurcated_a), size=.4))

set_plot_title <- function(xdim, ydim){
	ggtitle(
		paste(
			"Hydrogen Bonds ", xdim, " vs ", ydim, " by Salt Bridge Clusters", sep=""))
}

set_plot_id <- function(xdim, ydim){
	paste("hbond_geo_dim_2d_dARG_aCXL_by_salt_bridge_clusters", xdim, ydim, sep="_")
}


generate_plots <- function(sub_f){

	plot_id <- set_plot_id("cosBAH", "AHdist")
	ggplot(data=sub_f, aes(x=cosBAH, y=AHdist)) + plot_parts +
		scale_x_cosBAH + scale_y_AHdist +
		set_plot_title("cosBAH", "AHdist")
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

	plot_id <- set_plot_id("cosAHD", "AHdist")
	p <- ggplot(data=sub_f, aes(x=cosAHD, y=AHdist)) + plot_parts +
		scale_x_cosAHD + scale_y_AHdist +
		set_plot_title("cosAHD", "AHdist")
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

	plot_id <- set_plot_id("CHI", "AHdist")
	p <- ggplot(data=sub_f, aes(x=chi, y=AHdist)) + plot_parts +
		scale_x_chi + scale_y_AHdist +
		set_plot_title("CHI", "AHdist")
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

	plot_id <- set_plot_id("cosAHD", "cosBAH")
	p <- ggplot(data=sub_f, aes(x=cosAHD, y=cosBAH)) + plot_parts +
		scale_x_cosAHD + scale_y_cosBAH +
		set_plot_title("cosAHD", "cosBAH")
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

        plot_id <- set_plot_id("CHI", "cosBAH")
	p <- ggplot(data=sub_f, aes(x=chi, y=cosBAH)) + plot_parts +
		scale_x_chi + scale_y_cosBAH +
		set_plot_title("CHI", "cosBAH")
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

	plot_id <- set_plot_id("CHI", "cosAHD")
	p <- ggplot(data=sub_f, aes(x=chi, y=cosAHD)) + plot_parts +
		scale_x_chi + scale_y_cosAHD +
		set_plot_title("CHI", "cosAHD")
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)


	plot_id <- set_plot_id("CHI", "BAH")
	p <- ggplot(data=sub_f, aes(x=BAH_CHI_x, y=BAH_CHI_y)) + plot_parts +
		set_plot_title("CHI", "BAH")
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
        

}
generate_plots(f)

})) # end FeaturesAnalysis
