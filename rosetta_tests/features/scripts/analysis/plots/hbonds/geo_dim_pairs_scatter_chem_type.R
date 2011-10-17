# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

check_setup()

sele <-"
SELECT
  geom.AHdist, geom.cosBAH, geom.cosAHD, geom.chi,
  hblj.don_acc_atrE + hblj.don_acc_base_atrE +
    h_acc_atrE + h_acc_base_atrE AS lj_atrE,
  hblj.don_acc_repE + hblj.don_acc_base_repE +
    h_acc_repE + h_acc_base_repE AS lj_repE,
  hblj.don_acc_solv + hblj.don_acc_base_solv +
    h_acc_solv + h_acc_base_solv AS lj_solv,
  acc_site.HBChemType AS acc_chem_type,
  don_site.HBChemType AS don_chem_type
FROM
  hbond_geom_coords AS geom,
  hbonds AS hbond,
  hbond_sites AS don_site,
  hbond_sites AS acc_site,
  hbond_lennard_jones AS hblj
WHERE
  hbond.struct_id = geom.struct_id AND
  hbond.hbond_id =  geom.hbond_id AND
  hbond.struct_id = don_site.struct_id AND
  hbond.don_id = don_site.site_id AND
  hbond.struct_id = acc_site.struct_id AND
  hbond.acc_id = acc_site.site_id AND
  hbond.struct_id = hblj.struct_id AND
  hbond.hbond_id = hblj.hbond_id;"
f <- query_sample_sources(sample_sources, sele)

# This is deprecated please use the hbond_chem_types table for the lables instead
# Order the plots better and give more descriptive labels
f$don_chem_type <- factor(f$don_chem_type,
	levels = c("hbdon_IMD", "hbdon_IME", "hbdon_GDE", "hbdon_GDH",
		"hbdon_AHX", "hbdon_HXL", "hbdon_IND", "hbdon_AMO", "hbdon_CXA", "hbdon_PBA"),
	labels = c("dIMD: h", "dIME: h", "dGDE: r", "dGDH: r",
		"dAHX: y", "dHXL: s,t", "dIND: w", "dAMO: k", "dCXA: n,q", "dPBA: bb"))

# This is deprecated please use the hbond_chem_types table for the lables instead
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
	stat_density2d(size=.2))

set_plot_title <- function(xdim, ydim, ss_id){
	opts(title =
		paste("Hydrogen Bonds ", xdim, " vs ", ydim, "  ss_id: ", ss_id, sep=""))
}

plot_each_ss <- function(sub_f){
	ss_id <- sub_f$sample_source[1]
	ss <- sample_sources[sample_sources$sample_source == ss_id,]

	plot_id <- "geo_dim_pairs_scatter_cosBAH_AHdist_chem_type"
	ggplot(data=sub_f, aes(x=cosBAH, y=AHdist)) + plot_parts +
		scale_x_cosBAH + scale_y_AHdist + set_plot_title("cosBAH", "AHdist", ss_id)
	save_plots(plot_id, ss, output_dir, output_formats)

	plot_id <- "geo_dim_pairs_scatter_cosAHD_AHdist_chem_type"
	p <- ggplot(data=sub_f, aes(x=cosAHD, y=AHdist)) + plot_parts +
		scale_x_cosAHD + scale_y_AHdist + set_plot_title("cosAHD", "AHdist", ss_id)
	save_plots(plot_id, ss, output_dir, output_formats)

	plot_id <- "geo_dim_pairs_scatter_chi_AHdist_chem_type"
	p <- ggplot(data=sub_f, aes(x=chi, y=AHdist)) + plot_parts +
		scale_x_chi + scale_y_AHdist + set_plot_title("CHI", "AHdist", ss_id)
	save_plots(plot_id, ss, output_dir, output_formats)

	plot_id <- "geo_dim_pairs_scatter_cosAHD_cosBAH_chem_type"
	p <- ggplot(data=sub_f, aes(x=cosAHD, y=cosBAH)) + plot_parts +
		scale_x_cosAHD + scale_y_cosBAH + set_plot_title("cosAHD", "cosBAH", ss_id)
	save_plots(plot_id, ss, output_dir, output_formats)

	plot_id <- "geo_dim_pairs_scatter_chi_cosBAH_chem_type"
	p <- ggplot(data=sub_f, aes(x=chi, y=cosBAH)) + plot_parts +
		scale_x_chi + scale_y_cosBAH + set_plot_title("CHI", "cosBAH", ss_id)
	save_plots(plot_id, ss, output_dir, output_formats)

	plot_id <- "geo_dim_pairs_scatter_chi_cosAHD_chem_type"
	p <- ggplot(data=sub_f, aes(x=chi, y=cosAHD)) + plot_parts +
		scale_x_chi + scale_y_cosAHD + set_plot_title("CHI", "cosAHD", ss_id)
	save_plots(plot_id, ss, output_dir, output_formats)
}

runtime <- system.time(d_ply(f, .(sample_source), .fun=plot_each_ss))
print(paste("Plot Generation Time: ", runtime, sep=""))

