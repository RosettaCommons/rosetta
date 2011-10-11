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
  acc_site.HBChemType AS acc_chem_type,
  don_site.HBChemType AS don_chem_type
FROM
  hbond_geom_coords AS geom,
  hbonds AS hbond,
  hbond_sites AS don_site,
  hbond_sites AS acc_site
WHERE
  hbond.struct_id = geom.struct_id AND
  hbond.hbond_id =  geom.hbond_id AND
  hbond.struct_id = don_site.struct_id AND
  hbond.don_id = don_site.site_id AND
  hbond.struct_id = acc_site.struct_id AND
  hbond.acc_id = acc_site.site_id;"
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
	facet_grid(don_chem_type ~ acc_chem_type))

scale_x_AHdist <- scale_y_continuous(
	expression(paste('Acceptor -- Hydrogen Distance (', ring(A), ')')),
	limits=c(1.4, 3), breaks=c(1.4, 1.8, 2.2, 2.6, 3))

scale_y_AHdist <- scale_y_continuous(
	expression(paste('Acceptor -- Hydrogen Distance (', ring(A), ')')),
	limits=c(1.4, 3), breaks=c(1.4, 1.8, 2.2, 2.6, 3))

scale_x_cosBAH <- scale_x_continuous(
	"cos(Base -- Acceptor -- Hydrogen)",
	limit=c(-.4,1), breaks=c(-.4 -.2, 0, .2, .4, .6, .8, 1))

scale_y_cosBAH <- scale_y_continuous(
	"cos(Base -- Acceptor -- Hydrogen)",
	limit=c(-.4,1), breaks=c(-.4, -.2, 0, .2, .4, .6, .8, 1))

scale_x_cosAHD <- scale_x_continuous(
	"cos(Acceptor -- Hydrogen -- Donor)",
	limit=c(0,1), breaks=c(0, .2, .4, .6, .8, 1))

scale_y_cosAHD <- scale_y_continuous(
	"cos(Acceptor -- Hydrogen -- Donor)",
	limit=c(0,1), breaks=c(0, .2, .4, .6, .8, 1))

scale_x_chi <- scale_x_continuous(
	"Base -- Acceptor Torsion (Radians)",
	limit=c(0,2*pi), breaks=c(0, pi/3, pi*2/3, pi, pi*4/3, pi*5/3, 2*pi))

scale_y_chi <- scale_y_continuous(
	"Base -- Acceptor Torsion (Radians)",
	limit=c(0,2*pi), breaks=c(0, pi/3, pi*2/3, pi, pi*4/3, pi*5/3, 2*pi))

fill_scale_compression <- 100
plot_each_ss <- function(sub_f){
	ss_id <- sub_f$sample_source[1]
	ss <- sample_sources[sample_sources$sample_source == ss_id,]

	plot_id <- "geo_dim_pairs_density_cosBAH_AHdist_chem_type"
	p <- ggplot(data=sub_f) + plot_parts +
		stat_density2d(aes(x=cosBAH, y=AHdist, fill=log(..density..+100)), geom="tile", contour=FALSE) +
			scale_x_cosBAH + scale_y_AHdist +
			scale_fill_gradientn('Density', colour=jet.colors(10)) +
			opts(plot.background = theme_rect(colour = "#00007F")) +
		opts(title = paste("Hydrogen Bonds cosBAH vs AHdist  ss_id: ", ss_id, sep="_"))
	save_plots(plot_id, ss, output_dir, output_formats)

	plot_id <- "geo_dim_pairs_density_cosAHD_AHdist_chem_type"
	p <- ggplot(data=sub_f) + plot_parts +
		stat_density2d(aes(x=cosAHD, y=AHdist, fill=log(..density..+100)), geom="tile", contour=FALSE) +
			scale_x_cosAHD + scale_y_AHdist +
			scale_fill_gradientn('Density', colour=jet.colors(10)) +
			opts(plot.background = theme_rect(colour = "#00007F")) +
		opts(title = paste("Hydrogen Bonds cosAHD vs AHdist  ss_id: ", ss_id, sep="_"))
	save_plots(plot_id, ss, output_dir, output_formats)

	plot_id <- "geo_dim_pairs_density_chi_AHdist_chem_type"
	p <- ggplot(data=sub_f) + plot_parts +
		stat_density2d(aes(x=chi, y=AHdist, fill=log(..density..+100)), geom="tile", contour=FALSE) +
			scale_x_chi + scale_y_AHdist +
			scale_fill_gradientn('Density', colour=jet.colors(10)) +
			opts(plot.background = theme_rect(colour = "#00007F")) +
			opts(title = paste("Hydrogen Bonds CHI vs AHdist  ss_id: ", ss_id, sep="_"))
	save_plots(plot_id, ss, output_dir, output_formats)

	plot_id <- "geo_dim_pairs_density_cosAHD_cosBAH_chem_type"
	p <- ggplot(data=sub_f) + plot_parts +
		stat_density2d(aes(x=cosAHD, y=cosBAH, fill=log(..density..+100)), geom="tile", contour=FALSE) +
			scale_x_cosAHD + scale_y_cosBAH +
			scale_fill_gradientn('Density', colour=jet.colors(10)) +
			opts(plot.background = theme_rect(colour = "#00007F")) +
			opts(title = paste("Hydrogen Bonds cosAHD vs cosBAH  ss_id: ", ss_id, sep="_"))
	save_plots(plot_id, ss, output_dir, output_formats)

	plot_id <- "geo_dim_pairs_density_chi_cosBAH_chem_type"
	p <- ggplot(data=sub_f) + plot_parts +
		stat_density2d(aes(x=chi, y=cosBAH, fill=log(..density..+100)), geom="tile", contour=FALSE) +
			scale_x_chi + scale_y_cosBAH +
			scale_fill_gradientn('Density', colour=jet.colors(10)) +
			opts(plot.background = theme_rect(colour = "#00007F")) +
		opts(title = paste("Hydrogen Bonds CHI vs cosBAH  ss_id: ", ss_id, sep="_"))
	save_plots(plot_id, ss, output_dir, output_formats)

	plot_id <- "geo_dim_pairs_density_chi_cosAHD_chem_type"
	p <- ggplot(data=sub_f) + plot_parts +
		stat_density2d(aes(x=chi, y=cosAHD, fill=log(..density..+100)), geom="tile", contour=FALSE) +
			scale_x_chi + scale_y_cosAHD +
			scale_fill_gradientn('Density', colour=jet.colors(10)) +
			opts(plot.background = theme_rect(colour = "#00007F")) +
		opts(title = paste("Hydrogen Bonds CHI vs cosAHD  ss_id: ", ss_id, sep="_"))
	save_plots(plot_id, ss, output_dir, output_formats)
}

runtime <- system.time(d_ply(f, .(sample_source), .fun=plot_each_ss))
print(paste("Plot Generation Time: ", runtime, sep=""))

