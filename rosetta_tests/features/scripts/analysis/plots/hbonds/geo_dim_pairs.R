# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

check_setup()

do_query <- function(dOH_clause){
	sele <- paste("
	SELECT AHdist, cosBAH, cosAHD,
         acc_chem_type,
         CASE acc_chem_type WHEN 'hbacc_PBA' THEN 'BB Acceptor'
                                             ELSE 'SC Acceptor' END as bb_acc
	FROM  (SELECT geo.AHdist, geo.cosBAH, geo.cosAHD,
   	     	      don.HBChemType AS don_chem_type,
                acc.HBChemType AS acc_chem_type
   		   FROM   hbonds AS hb,
   		          hbond_sites AS don,
     	          hbond_sites AS acc,
     	          hbond_geom_coords AS geo
   		   WHERE  geo.struct_id = hb.struct_id AND geo.hbond_id = hb.hbond_id AND
   		          don.struct_id = hb.struct_id AND don.site_id  = hb.don_id AND
   		          acc.struct_id = hb.struct_id AND acc.site_id  = hb.acc_id AND
   		          ", dOH_clause, ");", sep="")

	query_sample_sources(sample_sources, sele)
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
	facet_wrap( ~ acc_chem_type, ncol=1),
	labs(y=expression(paste('(Acceptor -- Proton Distance)^2 (', ring(A), ')'))))


plot_parts_window <- list(
  theme_bw(),
  aes(y=-log(y+1), colour=AHdist_range, group=AHdist_range),
	geom_line(),
	geom_indicator(aes(indicator=counts, colour=AHdist_range)),
	facet_wrap( ~ acc_chem_type, ncol=2),
	labs(y="-log(FeatureDensity+1)"),
  scale_colour_brewer(palette="Spectral"))


################## dOH_1 ####################################
f <- do_query(
	"don.HBChemType='hbdon_AHX' AND
	(acc.HBChemType='hbacc_AHX' OR acc.HBChemType='hbacc_HXL' OR
	 acc.HBChemType='hbacc_CXA' OR acc.HBChemType='hbacc_CXL' OR
   acc.HBChemType='hbacc_PBA')")

plot_id <- "geo_dim_scatter_dOH_1_cosBAH"
ggplot(f, aes(x=cosBAH)) + plot_parts_scatter +
  labs(x="cos(Base -- Acceptor -- Hydrogen)") +
	opts(title = "Hydrogen Bonds cosBAH vs AHdist for AHX donors with AHX, HXL, CXA, CXL, or PBA acceptors")
save_plots(plot_id, sample_sources, output_dir, output_formats)

plot_id <- "geo_dim_scatter_dOH_1_cosAHD"
ggplot(f, aes(x=cosAHD)) + plot_parts_scatter +
  labs(x="cos(Acceptor -- Hydrogen -- Donor)") +
	opts(title = "Hydrogen Bonds cosAHD vs AHdist for AHX donors with AHX, HXL, CXA, CXL, or PBA acceptors")
save_plots(plot_id, sample_sources, output_dir, output_formats)

quantile_ranges <-
    data.frame(xmin=c(1.579960, 1.653311, 1.676390, 1.698214, 1.717584, 1.733008, 1.748605, 1.768427, 1.791709, 1.863508),
               xmax=c(1.733008, 1.748605, 1.768427, 1.791709, 1.863508, 2.156245, 2.999787, 2.999787, 2.999787, 2.999787))



plot_id <- "geo_dim_sliding_window_dOH_1_cosBAH"
dens <- estimate_density(f, "cosBAH", quantile_ranges)
ggplot(dens, aes(x=acos(x)*180/pi)) + plot_parts_window +
	labs(x='Base -- Acceptor -- Hydrogen (degrees)') +
	opts(title = "Hydrogen Bonds cosBAH by AHdist windows for AHX donors with AHX, HXL, CXA, CXL, or PBA acceptors\nnormalized for equal weight per unit distance")
save_plots(plot_id, sample_sources, output_dir, output_formats)

plot_id <- "geo_dim_sliding_window_dOH_1_cosAHD"
dens <- estimate_density(f, "cosAHD", quantile_ranges)
ggplot(dens, aes(x=acos(x)*180/pi)) + plot_parts_window +
  labs(x='Acceptor -- Hydrogen -- Donor (degrees)') +
	opts(title = "Hydrogen Bonds cosAHD by AHdist windows for AHX donors with AHX, HXL, CXA, CXL, or PBA acceptors\nnormalized for equal weight per unit distance")
save_plots(plot_id, sample_sources, output_dir, output_formats)


################## dOH_2 ####################################
f <- do_query(
	"don.HBChemType='hbdon_HXL' AND
	(acc.HBChemType='hbacc_AHX' OR acc.HBChemType='hbacc_HXL' OR
	 acc.HBChemType='hbacc_CXA' OR acc.HBChemType='hbacc_CXL' OR
   acc.HBChemType='hbacc_PBA')")


plot_id <- "geo_dim_scatter_dOH_2_cosBAH"
ggplot(f, aes(x=cosBAH)) + plot_parts_scatter +
  labs(x="cos(Base -- Acceptor -- Hydrogen)") +
	opts(title = "Hydrogen Bonds cosBAH vs AHdist for HXL donors with AHX, HXL, CXA, CXL, or PBA acceptors")
save_plots(plot_id, sample_sources, output_dir, output_formats)

plot_id <- "geo_dim_scatter_dOH_2_cosAHD"
ggplot(f, aes(x=cosAHD)) + plot_parts_scatter +
  labs(x="cos(Acceptor -- Hydrogen -- Donor)") +
	opts(title = "Hydrogen Bonds cosAHD vs AHdist for HXL donors with AHX, HXL, CXA, CXL, or PBA acceptors")
save_plots(plot_id, sample_sources, output_dir, output_formats)

quantile_ranges <-
    data.frame(xmin=c(1.577235, 1.674924, 1.714609, 1.737923, 1.760393, 1.788916, 1.825436, 1.877424, 1.956450, 2.100728),
               xmax=c(1.788916, 1.825436, 1.877424, 1.956450, 2.100728, 2.570788, 2.999943, 2.999943, 2.999943, 2.999943))

plot_id <- "geo_dim_sliding_window_dOH_2_cosBAH"
dens <- estimate_density(f, "cosBAH", quantile_ranges)
ggplot(dens, aes(x=acos(x)*180/pi)) + plot_parts_window +
	labs(x='Base -- Acceptor -- Hydrogen (degrees)') +
  opts(title = "Hydrogen Bonds cosBAH by AHdist windows for HXL donors with AXL, HXL, CXA, CXL, or PBA acceptors\nnormalized for equal weight per unit distance")
save_plots(plot_id, sample_sources, output_dir, output_formats)

plot_id <- "geo_dim_sliding_window_dOH_2_cosAHD"
dens <- estimate_density(f, "cosAHD", quantile_ranges)
ggplot(dens, aes(x=acos(x)*180/pi)) + plot_parts_window +
  labs(x='Acceptor -- Hydrogen -- Donor (degrees)') +
  opts(title = "Hydrogen Bonds cosAHD by AHdist windows for HXL donors with AXL, HXL, CXA, CXL, or PBA acceptors\nnormalized for equal weight per unit distance")
save_plots(plot_id, sample_sources, output_dir, output_formats)


################## dOH_3 ####################################
f <- do_query(
	"(don.HBChemType='hbdon_AHX' OR don.HBChemType='hbdon_HXL') AND
	 (acc.HBChemType='hbacc_IMD' OR acc.HBChemType='hbacc_IME')")

plot_id <- "geo_dim_scatter_dOH_3_cosBAH"
ggplot(f, aes(x=cosBAH)) + plot_parts_scatter +
  labs(x="cos(Base -- Acceptor -- Hydrogen)") +
	opts(title = "Hydrogen Bonds cosBAH vs AHdist for AHX or HXL donors with IMD or IME acceptors")
save_plots(plot_id, sample_sources, output_dir, output_formats)

plot_id <- "geo_dim_scatter_dOH_3_cosAHD"
ggplot(f, aes(x=cosAHD)) + plot_parts_scatter +
  labs(x="cos(Acceptor -- Hydrogen -- Donor)") +
	opts(title = "Hydrogen Bonds cosAHD vs AHdist for AHX or HXL donors with IMD or IME acceptors")
save_plots(plot_id, sample_sources, output_dir, output_formats)


quantile_ranges <-
    data.frame(xmin=c(1.834432, 1.901671, 1.921783, 1.936392, 1.948883, 1.965945, 1.991784, 2.065187, 2.282928, 2.351500),
               xmax=c(1.965945, 1.991784, 2.065187, 2.282928, 2.351500, 2.448438, 2.944114, 2.944114, 2.944114, 2.944114))

plot_id <- "geo_dim_sliding_window_dOH_3_cosBAH"
dens <- estimate_density(f, "cosBAH", quantile_ranges)
ggplot(dens, aes(x=acos(x)*180/pi)) + plot_parts_window +
	labs(x='Base -- Acceptor -- Hydrogen (degrees)') +
  opts(title = "Hydrogen Bonds cosBAH by AHdist windows for AHX or HXL donors with IMD or IME acceptors\nnormalized for equal weight per unit distance")
save_plots(plot_id, sample_sources, output_dir, output_formats)

plot_id <- "geo_dim_sliding_window_dOH_3_cosAHD"
dens <- estimate_density(f, "cosAHD", quantile_ranges)
ggplot(dens, aes(x=acos(x)*180/pi)) + plot_parts_window +
  labs(x='Acceptor -- Hydrogen -- Donor (degrees)') +
  opts(title = "Hydrogen Bonds cosAHD by AHdist windows for AHX or HXL donors with IMD or IME acceptors\nnormalized for equal weight per unit distance")
save_plots(plot_id, sample_sources, output_dir, output_formats)
