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
id = "chi_kl_divergence",
filename = "scripts/analysis/statistics/hbonds/chi_kl_divergence.R",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self){

sele <-"
SELECT
	geom.chi,
	acc_site.HBChemType AS acc_chem_type,
	don_site.HBChemType AS don_chem_type
FROM
	hbond_geom_coords AS geom,
	hbonds AS hbond,
	hbond_sites AS don_site,
	hbond_sites AS acc_site
WHERE
	hbond.struct_id = geom.struct_id AND hbond.hbond_id   = geom.hbond_id AND
	hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
	hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
	(acc_site.HBChemType='hbacc_PBA' OR
		acc_site.HBChemType='hbacc_CXL' OR
		acc_site.HBChemType='hbacc_CXA') AND
	ABS(don_site.resNum - acc_site.resNum) > 5;"

f <- query_sample_sources(sample_sources, sele)
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

sc_f <- f[f$acc_chem_type != 'hbacc_PBA',]
dens <- estimate_density_1d_wrap(
	data = sc_f,
	ids = c("sample_source", "don_chem_type"),
	variable = "chi")

# This is a prototype for a multiplicative factor of the hbond energy
# to account for the chi angle preference.
cos_potential <- function(chi) -1*(cos( 2 * chi * pi/180 ) + 2 ) / 3

# Transform the potential energy to probability via the Boltzmann
# transformation where the partition function is over the chi angle
potential <- data.frame( x=seq(0,360,length.out=200))
potential$y_unmormalized <- exp(-1*cos_potential(potential$x))
potential$y <- exp(-1*cos_potential(potential$x))/(12.5803 *180/pi)

comp <-
  ldply(levels(dens$sample_source), function(ss1){
  ldply(levels(dens$sample_source), function(ss2){
    ldply(levels(dens$don_chem_type), function(dct){
      data.frame( sample_source1 = ss1, sample_source2 = ss2, don_chem_type = dct,
                 kl_div = 
      sum(dens[dens$sample_source == ss1 & dens$don_chem_type == dct, "y"] *
      	(log(dens[dens$sample_source == ss1 & dens$don_chem_type == dct, "y"]) -
	log(dens[dens$sample_source == ss2 & dens$don_chem_type == dct, "y"]))))
    })
  })
})

plot_id = "chi_kl_divergence_by_don_chem_type"
ggplot(comp) + theme_bw() +
  geom_bar(aes(x=don_chem_type, y=kl_div)) +
  facet_grid(sample_source1 ~ sample_source2) +
  opts(title = "Hydrogen Bonds CHI Angle KL Divergence with sequence separation > 5\n(normalized for equal volume per unit distance)") +
  coord_flip()
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

comp <-
  ldply(levels(dens$sample_source), function(ss1){
  ldply(levels(dens$sample_source), function(ss2){
    ldply(levels(dens$don_chem_type), function(dct){
      data.frame( sample_source1 = ss1, sample_source2 = ss2, don_chem_type = dct,
                 cross_entropy = 
      sum(dens[dens$sample_source == ss1 & dens$don_chem_type == dct, "y"] *
	log(dens[dens$sample_source == ss2 & dens$don_chem_type == dct, "y"])))
    })
  })
})

plot_id = "chi_cross_entropy_by_don_chem_type"
ggplot(comp) + theme_bw() +
  geom_bar(aes(x=don_chem_type, y=cross_entropy)) +
  facet_grid(sample_source1 ~ sample_source2) +
  opts(title = "Hydrogen Bonds CHI Angle Cross Entropy with sequence separation > 5\n(normalized for equal volume per unit distance)") +
  coord_flip()
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
