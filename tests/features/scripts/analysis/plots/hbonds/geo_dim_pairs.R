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
id = "geo_dim_pairs",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

do_query <- function(dOH_clause, sample_source){
	sele <- paste("
SELECT
	geo.AHdist, geo.cosBAH, geo.cosAHD,
	don.HBChemType AS don_chem_type,
	acc.HBChemType AS acc_chem_type
FROM
	hbonds AS hb,
	hbond_sites AS don,
	hbond_sites AS acc,
	hbond_geom_coords AS geo
WHERE
	geo.struct_id = hb.struct_id AND geo.hbond_id = hb.hbond_id AND
	don.struct_id = hb.struct_id AND don.site_id  = hb.don_id AND
	acc.struct_id = hb.struct_id AND acc.site_id  = hb.acc_id AND
	(", dOH_clause, ");", sep="")

	f <- query_sample_sources(sample_source, sele)

	# give more descriptive plot labels
	f$acc_chem_type <- factor(f$acc_chem_type,
		levels = c("hbacc_IMD", "hbacc_IME", "hbacc_AHX", "hbacc_HXL",
			"hbacc_CXA", "hbacc_CXL", "hbacc_PBA"),
		labels = c("aIMD: h", "aIME: h", "aAHX: y", "aHXL: s,t",
			"aCXA: n,q", "aCXL: d,e", "aPBA: bb"))
	f
}

# compute the density for each of the quantile ranges
estimate_density <- function(f, geo_dim, quantile_ranges, ...){
	adply(quantile_ranges, 1, function(range){
			 sub <- f[f$AHdist >= range$xmin[1] & f$AHdist < range$xmax[1],]
			 dens <- estimate_density_1d_logspline(
							 			sub, c("acc_chem_type", "sample_source"), geo_dim)
			 dens$AHdist_range <- factor(paste("[", round(range$xmin,2), ", ", round(range$xmax,2), ")", sep=""))
			 dens
	})
}

plot_parts_scatter <- list(
  theme_bw(),
  aes(y=AHdist^2),
	geom_point(size=.7),
	stat_density2d(size=.2),
	facet_wrap( ~ acc_chem_type, ncol=1),
	labs(y=expression(paste('(Acceptor -- Proton Distance)^2 (', ring(A), ')'))),
	scale_x_continuous(limits=c(0,1)),
	scale_y_continuous(limits=c(2.7,9)))


plot_parts_window <- list(
  theme_bw(),
  aes(y=-log(y+1), colour=AHdist_range, group=AHdist_range),
	geom_line(),
	geom_indicator(aes(indicator=counts, colour=AHdist_range, group=AHdist_range)),
	facet_wrap( ~ acc_chem_type, ncol=2),
	labs(y="-log(FeatureDensity+1)"),
  scale_colour_brewer(palette="Spectral"))

d_ply(sample_sources, .(sample_source), function(ss){
	ss_id <- ss$sample_source[1]
	################## dOH_1 ####################################
	f <- do_query(
		"don.HBChemType='hbdon_AHX' AND
		(acc.HBChemType='hbacc_AHX' OR acc.HBChemType='hbacc_HXL' OR
		 acc.HBChemType='hbacc_CXA' OR acc.HBChemType='hbacc_CXL' OR
	   acc.HBChemType='hbacc_PBA')", ss)

	plot_id <- paste("geo_dim_scatter_dOH_1_cosBAH", ss_id, sep="_")
	ggplot(f, aes(x=cosBAH)) + plot_parts_scatter +
	  labs(x="cos(Base -- Acceptor -- Hydrogen)") +
		ggtitle(paste("Hydrogen Bonds cosBAH vs AHdist for AHX donors with AHX, HXL, CXA, CXL, or PBA acceptors\nss_id: ", ss_id, sep=""))
	save_plots(self, plot_id, ss, output_dir, output_formats)

	plot_id <- paste("geo_dim_scatter_dOH_1_cosAHD", ss_id, sep="_")
	ggplot(f, aes(x=cosAHD)) + plot_parts_scatter +
	  labs(x="cos(Acceptor -- Hydrogen -- Donor)") +
		ggtitle(paste("Hydrogen Bonds cosAHD vs AHdist for AHX donors with AHX, HXL, CXA, CXL, or PBA acceptors\nss_id: ", ss_id, sep=""))
	save_plots(self, plot_id, ss, output_dir, output_formats)

	quantile_ranges <-
	    data.frame(xmin=c(1.579960, 1.653311, 1.676390, 1.698214, 1.717584, 1.733008, 1.748605, 1.768427, 1.791709, 1.863508),
	               xmax=c(1.733008, 1.748605, 1.768427, 1.791709, 1.863508, 2.156245, 2.999787, 2.999787, 2.999787, 2.999787))


	plot_id <- paste("geo_dim_sliding_window_dOH_1_cosBAH", ss_id, sep="_")
	dens <- estimate_density(f, "cosBAH", quantile_ranges)
	ggplot(dens, aes(x=acos(x)*180/pi)) + plot_parts_window +
		labs(x='Base -- Acceptor -- Hydrogen (degrees)') +
		ggtitle(paste("Hydrogen Bonds cosBAH by AHdist windows for AHX donors with AHX, HXL, CXA, CXL, or PBA acceptors\nnormalized for equal weight per unit distance\nss_id: ", ss_id, sep=""))
	save_plots(self, plot_id, ss, output_dir, output_formats)

	plot_id <- paste("geo_dim_sliding_window_dOH_1_cosAHD", ss_id, sep="_")
	dens <- estimate_density(f, "cosAHD", quantile_ranges)
	ggplot(dens, aes(x=acos(x)*180/pi)) + plot_parts_window +
	  labs(x='Acceptor -- Hydrogen -- Donor (degrees)') +
		ggtitle(paste("Hydrogen Bonds cosAHD by AHdist windows for AHX donors with AHX, HXL, CXA, CXL, or PBA acceptors\nnormalized for equal weight per unit distanc\nss_id: ", ss_id, sep=""))
	save_plots(self, plot_id, ss, output_dir, output_formats)


	################## dOH_2 ####################################
	f <- do_query(
		"don.HBChemType='hbdon_HXL' AND
		(acc.HBChemType='hbacc_AHX' OR acc.HBChemType='hbacc_HXL' OR
		 acc.HBChemType='hbacc_CXA' OR acc.HBChemType='hbacc_CXL' OR
	   acc.HBChemType='hbacc_PBA')", ss)

	plot_id <- paste("geo_dim_scatter_dOH_2_cosBAH", ss_id, sep="_")
	ggplot(f, aes(x=cosBAH)) + plot_parts_scatter +
	  labs(x="cos(Base -- Acceptor -- Hydrogen)") +
		ggtitle(paste("Hydrogen Bonds cosBAH vs AHdist for HXL donors with AHX, HXL, CXA, CXL, or PBA acceptors\nss_id: ", ss_id, sep=""))
	save_plots(self, plot_id, ss, output_dir, output_formats)

	plot_id <- paste("geo_dim_scatter_dOH_2_cosAHD", ss_id, sep="_")
	ggplot(f, aes(x=cosAHD)) + plot_parts_scatter +
	  labs(x="cos(Acceptor -- Hydrogen -- Donor)") +
		ggtitle(paste("Hydrogen Bonds cosAHD vs AHdist for HXL donors with AHX, HXL, CXA, CXL, or PBA acceptors\nss_id: ", ss_id, sep=""))
	save_plots(self, plot_id, ss, output_dir, output_formats)

	quantile_ranges <-
	    data.frame(xmin=c(1.577235, 1.674924, 1.714609, 1.737923, 1.760393, 1.788916, 1.825436, 1.877424, 1.956450, 2.100728),
	               xmax=c(1.788916, 1.825436, 1.877424, 1.956450, 2.100728, 2.570788, 2.999943, 2.999943, 2.999943, 2.999943))

	plot_id <- paste("geo_dim_sliding_window_dOH_2_cosBAH", ss_id, sep="_")
	dens <- estimate_density(f, "cosBAH", quantile_ranges)
	ggplot(dens, aes(x=acos(x)*180/pi)) + plot_parts_window +
		labs(x='Base -- Acceptor -- Hydrogen (degrees)') +
	  ggtitle(paste("Hydrogen Bonds cosBAH by AHdist windows for HXL donors with AXL, HXL, CXA, CXL, or PBA acceptors\nnormalized for equal weight per unit distance\nss_id: ", ss_id, sep=""))
	save_plots(self, plot_id, ss, output_dir, output_formats)

	plot_id <- paste("geo_dim_sliding_window_dOH_2_cosAHD", ss_id, sep="_")
	dens <- estimate_density(f, "cosAHD", quantile_ranges)
	ggplot(dens, aes(x=acos(x)*180/pi)) + plot_parts_window +
	  labs(x='Acceptor -- Hydrogen -- Donor (degrees)') +
	  ggtitle(paste("Hydrogen Bonds cosAHD by AHdist windows for HXL donors with AXL, HXL, CXA, CXL, or PBA acceptors\nnormalized for equal weight per unit distance\nss_id: ", ss_id, sep=""))
	save_plots(self, plot_id, ss, output_dir, output_formats)


	################## dOH_3 ####################################
	f <- do_query(
		"(don.HBChemType='hbdon_AHX' OR don.HBChemType='hbdon_HXL') AND
		 (acc.HBChemType='hbacc_IMD' OR acc.HBChemType='hbacc_IME')", ss)

	plot_id <- paste("geo_dim_scatter_dOH_3_cosBAH", ss_id, sep="_")
	ggplot(f, aes(x=cosBAH)) + plot_parts_scatter +
	  labs(x="cos(Base -- Acceptor -- Hydrogen)") +
		ggtitle(paste("Hydrogen Bonds cosBAH vs AHdist for AHX or HXL donors with IMD or IME acceptors\nss_id: ", ss_id, sep=""))
	save_plots(self, plot_id, ss, output_dir, output_formats)

	plot_id <- paste("geo_dim_scatter_dOH_3_cosAHD", ss_id, sep="_")
	ggplot(f, aes(x=cosAHD)) + plot_parts_scatter +
	  labs(x="cos(Acceptor -- Hydrogen -- Donor)") +
		ggtitle(paste("Hydrogen Bonds cosAHD vs AHdist for AHX or HXL donors with IMD or IME acceptors\nss_id: ", ss_id, sep=""))
	save_plots(self, plot_id, ss, output_dir, output_formats)

	quantile_ranges <-
	    data.frame(xmin=c(1.834432, 1.901671, 1.921783, 1.936392, 1.948883, 1.965945, 1.991784, 2.065187, 2.282928, 2.351500),
	               xmax=c(1.965945, 1.991784, 2.065187, 2.282928, 2.351500, 2.448438, 2.944114, 2.944114, 2.944114, 2.944114))

	plot_id <- paste("geo_dim_sliding_window_dOH_3_cosBAH", ss_id, sep="_")
	dens <- estimate_density(f, "cosBAH", quantile_ranges)
	ggplot(dens, aes(x=acos(x)*180/pi)) + plot_parts_window +
		labs(x='Base -- Acceptor -- Hydrogen (degrees)') +
	  ggtitle(paste( "Hydrogen Bonds cosBAH by AHdist windows for AHX or HXL donors with IMD or IME acceptors\nss_id: ", ss_id, sep=""))
	save_plots(self, plot_id, ss, output_dir, output_formats)

	plot_id <-paste( "geo_dim_sliding_window_dOH_3_cosAHD", ss_id, sep="_")
	dens <- estimate_density(f, "cosAHD", quantile_ranges)
	ggplot(dens, aes(x=acos(x)*180/pi)) + plot_parts_window +
	  labs(x='Acceptor -- Hydrogen -- Donor (degrees)') +
	  ggtitle(paste("Hydrogen Bonds cosAHD by AHdist windows for AHX or HXL donors with IMD or IME acceptors\nss_id: ", ss_id, sep=""))
	save_plots(self, plot_id, ss, output_dir, output_formats)


	################## aOHsc ####################################
	f <- do_query(
		"(acc.HBChemType='hbacc_AHX' OR acc.HBChemType='hbacc_HXL') AND
		 (don.HBChemType!='hbdon_PBA' AND don.HBChemType!='hbdon_PBA')", ss)

	plot_id <- paste("geo_dim_scatter_aOHsc_cosBAH", ss_id, sep="_")
	ggplot(f, aes(x=cosBAH)) + plot_parts_scatter +
	  labs(x="cos(Base -- Acceptor -- Hydrogen)") +
		ggtitle(paste("Hydrogen Bonds cosBAH vs AHdist for AHX or HXL acceptors with SC donors\nss_id: ", ss_id, sep=""))
	save_plots(self, plot_id, ss, output_dir, output_formats)

	plot_id <- paste("geo_dim_scatter_aOHsc_cosAHD", ss_id, sep="_")
	ggplot(f, aes(x=cosAHD)) + plot_parts_scatter +
	  labs(x="cos(Acceptor -- Hydrogen -- Donor)") +
		ggtitle(paste("Hydrogen Bonds cosAHD vs AHdist for AHX or HXL acceptors with SC donors\nss_id: ", ss_id, sep=""))
	save_plots(self, plot_id, ss, output_dir, output_formats)
})


})) # end FeaturesAnalysis
