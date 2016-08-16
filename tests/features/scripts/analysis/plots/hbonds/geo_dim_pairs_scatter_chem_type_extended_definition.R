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
id = "geo_dim_pairs_scatter_chem_type_extended_definition",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

source("scripts/analysis/plots/hbonds/hbond_geo_dim_scales.R")

sele <-"
SELECT
  geom.AHdist, geom.cosBAH, geom.cosAHD, geom.chi,
	don_atoms.base_x AS dx, don_atoms.base_y AS dy, don_atoms.base_z AS dz,
	acc_atoms.atm_x  AS ax, acc_atoms.atm_y  AS ay, acc_atoms.atm_z  AS az,
  acc_site.HBChemType AS acc_chem_type,
  don_site.HBChemType AS don_chem_type,
	hb.energy < 0 AS hb_exists
FROM
  hbond_geom_coords AS geom,
  hbonds AS hb,
  hbond_sites AS don_site,
  hbond_sites AS acc_site,
	hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb,
	hbond_site_atoms AS don_atoms, hbond_site_atoms AS acc_atoms
WHERE
  geom.struct_id = hb.struct_id AND
  geom.hbond_id =  hb.hbond_id AND
  don_site.struct_id = hb.struct_id AND
  don_site.site_id = hb.don_id AND
  acc_site.struct_id = hb.struct_id AND
  acc_site.site_id = hb.acc_id AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30 AND
	don_atoms.struct_id = hb.struct_id AND don_atoms.site_id = hb.don_id AND
	acc_atoms.struct_id = hb.struct_id AND acc_atoms.site_id = hb.acc_id;"
f <- query_sample_sources(sample_sources, sele)

f <- transform(f,
	ADdist = vector_distance(cbind(dx, dy, dz), cbind(ax, ay, az)))

# Filter to HBonds having A-D distance not greater than 3.5 Angstroms.
f <- f[f$ADdist <= 3.5,]


f$hb_exists <- factor(f$hb_exists)

# Order the plots better and give more descriptive labels
f$don_chem_type <- factor(f$don_chem_type,
	levels = c("hbdon_IMD", "hbdon_IME", "hbdon_GDE", "hbdon_GDH",
		"hbdon_AHX", "hbdon_HXL", "hbdon_IND", "hbdon_AMO", "hbdon_CXA", "hbdon_PBA"),
	labels = c("dIMD: h", "dIME: h", "dGDE: r", "dGDH: r",
		"dAHX: y", "dHXL: s,t", "dIND: w", "dAMO: k", "dCXA: n,q", "dPBA: bb"))

# Order the plots better and give more descriptive labels
f$acc_chem_type <- factor(f$acc_chem_type,
	levels = c("hbacc_IMD", "hbacc_IME", "hbacc_AHX", "hbacc_HXL",
		"hbacc_CXA", "hbacc_CXL", "hbacc_PBA"),
	labels = c("aIMD: h", "aIME: h", "aAHX: y", "aHXL: s,t",
		"aCXA: n,q", "aCXL: d,e", "aPBA: bb"))

plot_parts <- list(
	theme_bw(),
	facet_grid(don_chem_type ~ acc_chem_type),
	geom_point(size=.4),
	stat_density2d(size=.2, colour="black"))

set_plot_title <- function(xdim, ydim, ss_id){
	ggtitle(
		paste("Hydrogen Bonds ", xdim, " vs ", ydim, " BFact < 30; ss_id: ", ss_id, sep=""))
}

plot_each_ss <- function(sub_f){
	ss_id <- sub_f$sample_source[1]
	ss <- sample_sources[sample_sources$sample_source == ss_id,]

	plot_id <- "geo_dim_pairs_scatter_extended_definition_cosBAH_AHdist_chem_type"
	ggplot(data=sub_f, aes(x=cosBAH, y=AHdist, colour=hb_exists)) + plot_parts +
		scale_x_cosBAH + scale_y_AHdist + set_plot_title("cosBAH", "AHdist", ss_id)
	save_plots(self, plot_id, ss, output_dir, output_formats)

	plot_id <- "geo_dim_pairs_scatter_extended_definition_cosAHD_AHdist_chem_type"
	p <- ggplot(data=sub_f, aes(x=cosAHD, y=AHdist, colour=hb_exists)) + plot_parts +
		scale_x_cosAHD + scale_y_AHdist + set_plot_title("cosAHD", "AHdist", ss_id)
	save_plots(self, plot_id, ss, output_dir, output_formats)

	plot_id <- "geo_dim_pairs_scatter_extended_definition_chi_AHdist_chem_type"
	p <- ggplot(data=sub_f, aes(x=chi, y=AHdist, colour=hb_exists)) + plot_parts +
		scale_x_chi + scale_y_AHdist + set_plot_title("CHI", "AHdist", ss_id)
	save_plots(self, plot_id, ss, output_dir, output_formats)

	plot_id <- "geo_dim_pairs_scatter_extended_definition_cosAHD_cosBAH_chem_type"
	p <- ggplot(data=sub_f, aes(x=cosAHD, y=cosBAH, colour=hb_exists)) + plot_parts +
		scale_x_cosAHD + scale_y_cosBAH + set_plot_title("cosAHD", "cosBAH", ss_id)
	save_plots(self, plot_id, ss, output_dir, output_formats)

	plot_id <- "geo_dim_pairs_scatter_extended_definition_chi_cosBAH_chem_type"
	p <- ggplot(data=sub_f, aes(x=chi, y=cosBAH, colour=hb_exists)) + plot_parts +
		scale_x_chi + scale_y_cosBAH + set_plot_title("CHI", "cosBAH", ss_id)
	save_plots(self, plot_id, ss, output_dir, output_formats)

	plot_id <- "geo_dim_pairs_scatter_extended_definition_chi_cosAHD_chem_type"
	p <- ggplot(data=sub_f, aes(x=chi, y=cosAHD, colour=hb_exists)) + plot_parts +
		scale_x_chi + scale_y_cosAHD + set_plot_title("CHI", "cosAHD", ss_id)
	save_plots(self, plot_id, ss, output_dir, output_formats)
}

runtime <- system.time(d_ply(f, .(sample_source), .fun=plot_each_ss))
print(paste("Plot Generation Time: ", runtime, sep=""))

})) # end FeaturesAnalysis
