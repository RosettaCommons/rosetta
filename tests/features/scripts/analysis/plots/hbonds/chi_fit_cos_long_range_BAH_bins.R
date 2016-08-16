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
id = "chi_fit_cos_long_range_BAH_bins",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

d_ply(sample_sources, .variables=("sample_source"), function(sample_source){
	ss <- sample_source[1,"sample_source"]

	sele <-"
	SELECT
		geom.chi, geom.cosBAH,
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
		hbond.acc_id = acc_site.site_id AND
		(acc_site.HBChemType='hbacc_PBA' OR
			acc_site.HBChemType='hbacc_CXL' OR
			acc_site.HBChemType='hbacc_CXA') AND
		ABS(don_site.resNum - acc_site.resNum) > 5;"

	f <- query_sample_sources(sample_source, sele)
	f$chi <- (f$chi*180/pi)

	# Order the plots better and give more descriptive labels
	f$don_chem_type <- factor(f$don_chem_type,
		levels=c(
			"hbdon_IMD",
	  	"hbdon_GDE",
	    "hbdon_HXL",
	    "hbdon_PBA",
	    "hbdon_IME",
	    "hbdon_GDH",
	    "hbdon_AHX",
	    "hbdon_IND",
	    "hbdon_CXA",
	    "hbdon_AMO"),
		labels = c(
			"Donor IMD: his",
			"Donor GDE: arg",
			"Donor HXL: ser, thr",
			"Donor PBA: backbone",
			"Donor IME: his",
			"Donor GDH: arg",
			"Donor AHX: tyr",
			"Donor IND: trp",
			"Donor CXA: gln, asn",
			"Donor AMO: lys"))


	# when the BAH angle is straigth (i.e. close to 0), the chi angle is less interesting
	f$BAH_deg <- acos(f$cosBAH)*180/pi
  f$BAH_bins <- cut(f$BAH_deg, seq(0,90, length=9))


	sc_f <- f[f$acc_chem_type != 'hbacc_PBA',]
	dens <- estimate_density_1d_wrap(
		data = sc_f,
		ids = c("BAH_bins", "don_chem_type"),
		variable = "chi")

	# This is a prototype for a multiplicative factor of the hbond energy
	# to account for the chi angle preference.
	cos_potential <- function(chi) -1*(cos( 2 * chi * pi/180 ) + 2 ) / 3

	# Transform the potential energy to probability via the Boltzmann
	# transformation where the partition function is over the chi angle
	potential <- data.frame( x=seq(0,360,length.out=200))
	potential$y_unnormalized <- exp(-1*cos_potential(potential$x))
	potential$y <- exp(-1*cos_potential(potential$x))/(12.5803 *180/pi)

	plot_id = "chi_chem_type_sp2_long_range_BAH_bins"
	ggplot(data=dens) + theme_bw() +
		geom_line(data=potential, aes(x=x, y=y), size=1.2, colour="darkgray") +
		geom_line(aes(x=x, y=y, colour=BAH_bins)) +
		geom_indicator(aes(indicator=counts, colour=BAH_bins, group=BAH_bins)) +
		facet_wrap( ~ don_chem_type) +
		ggtitle(paste("Hydrogen Bonds CHI Angle for Sidechain sp2 Acceptors with Sequence Separation > 5\n(normalized for equal volume per Unit Distance) Sample Source: ", ss, sep="")) +
		scale_x_continuous('Acceptor Base -- Acceptor Torsion (degrees)', breaks=c(90,270)) +
		scale_y_continuous('Feature Density') +
		theme(legend.position=c(.58,.35)) +
		theme(legend.justification=c("left", "top"))
	save_plots(self, plot_id, sample_source, output_dir, output_formats)



	plot_id = "chi_chem_type_sp2_by_acc_long_range_BAH_bins"
	dens <- estimate_density_1d_wrap(
		data = f,
		ids = c("BAH_bins", "acc_chem_type", "don_chem_type"),
		variable = "chi")
	ggplot(data=dens) + theme_bw() +
		geom_line(aes(x=x, y=y, colour=BAH_bins)) +
		geom_indicator(aes(indicator=counts, colour=BAH_bins, group=BAH_bins)) +
		facet_grid(don_chem_type ~ acc_chem_type) +
		ggtitle(paste("Hydrogen Bonds CHI Angle by Acceptor types with sequence separation > 5\n(normalized for equal volume per unit distance) Sample Source: ", ss, sep="")) +
		scale_x_continuous('Acceptor Base -- Acceptor Torsion (degrees)', breaks=c(90,270)) +
		scale_y_continuous('Feature Density', breaks=c(90,270))
	save_plots(self, plot_id, sample_source, output_dir, output_formats)
})


})) # end FeaturesAnalysis
