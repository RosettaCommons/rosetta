# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

check_setup()

morse_fn <- function(x, D_a, a, r_0, min_e){
	D_a*(1+exp(-2*a*(x-r_0))-2*exp(-a*(x-r_0)))+min_e
}


estimate_densities <- function(dOH_clause){
	sele <- paste("
	SELECT AHdist,
	       REPLACE(don_chem_type, 'hbdon_', '') AS d,
	       REPLACE(don_chem_type, 'hbacc_', '') AS a
	FROM  (SELECT geo.AHdist,
   	     	      don.HBChemType AS don_chem_type,
                acc.HBChemType AS acc_chem_type
   		   FROM   hbonds AS hb,
   		          hbond_sites AS don,
     	          hbond_sites AS acc,
     	          hbond_geom_coords AS geo
   		   WHERE  geo.struct_id = hb.struct_id AND geo.hbond_id = hb.hbond_id AND
   		          don.struct_id = hb.struct_id AND don.site_id  = hb.don_id AND
   		          acc.struct_id = hb.struct_id AND acc.site_id  = hb.acc_id AND
                geo.AHdist < 2.6 AND
   		          ", dOH_clause, ");", sep="")

	f <- query_sample_sources(sample_sources, sele)

	cat("estimate_density\n")
	dens <- estimate_density_1d(f, c("sample_source"), "AHdist", radial_3d_normalization)
	dens$energy = -log(dens$y)

	cat("computting morse regression\n")
	ddply(dens, .(sample_source), function(d){
		d$morse <- 0
		tryCatch({
			m1 <- nls(energy ~ morse_fn(x, D_a, a, r_0, min_e), d,
				start=list(D_a=7, a=1.8, r_0=1.75, min_e=-1), algorithm="port",
				trace=TRUE, weights=y)
			d$morse <- predict(m1,dens$x)

			},
			error=function(e){
				cat("Unable to fit morse function to ", as.character(d$sample_source[1]), "\n")
				print(paste("failed with error: ", e, sep=""))
			})
		d
	})
}

plot_parts <- list(
  theme_bw(),
	geom_line(aes(x=x, y=energy)),
	geom_line(aes(x=x, y=morse), colour="blue"),
	facet_wrap( ~ sample_source, ncol=1),
	labs(x=expression(paste('Acceptor -- Proton Distance (', ring(A), ')')),
		y="-log(FeatureDensity)"),
	scale_y_continuous(limits=c(-3,5), breaks=((0:8)-3)),
	scale_x_continuous(limits=c(1.5,2.6), breaks=c(1.6, 1.9, 2.2, 2.5)))


plot_id <- "AHdist_dOH_1_fit_energy"
dens <- estimate_densities(
	"don.HBChemType='hbdon_AHX' AND
	(acc.HBChemType='hbacc_AHX' OR acc.HBChemType='hbacc_HXL' OR
	 acc.HBChemType='hbacc_CXA' OR acc.HBChemType='hbacc_CXL' OR
   acc.HBChemType='hbacc_PBA')")
ggplot(dens) + plot_parts +
	opts(title = "Hydrogen Bonds A-H Distance by Chemical Type for AHX donors with AHX, HXL, CXA, CXL, or PBA acceptors\nnormalized for equal weight per unit distance")
save_plots(plot_id, sample_sources, output_dir, output_formats)

plot_id <- "AHdist_dOH_2_fit_energy"
dens <- estimate_densities(
	"don.HBChemType='hbdon_HXL' AND
	(acc.HBChemType='hbacc_AHX' OR acc.HBChemType='hbacc_HXL' OR
	 acc.HBChemType='hbacc_CXA' OR acc.HBChemType='hbacc_CXL' OR
   acc.HBChemType='hbacc_PBA')")
ggplot(dens) + plot_parts +
	opts(title = "Hydrogen Bonds A-H Distance by Chemical Type for HXL donors with AXL, HXL, CXA, CXL, or PBA acceptors\nnormalized for equal weight per unit distance")
save_plots(plot_id, sample_sources, output_dir, output_formats)

plot_id <- "AHdist_dOH_3_fit_energy"
dens <- estimate_densities(
	"(don.HBChemType='hbdon_AHX' OR don.HBChemType='hbdon_HXL') AND
	 (acc.HBChemType='hbacc_IMD' OR acc.HBChemType='hbacc_IME')")
ggplot(dens) + plot_parts +
	opts(title = "Hydrogen Bonds A-H Distance by Chemical Type for AHX or HXL donors with IMD or IME acceptors\nnormalized for equal weight per unit distance")
save_plots(plot_id, sample_sources, output_dir, output_formats)
