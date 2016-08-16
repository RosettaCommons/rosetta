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
id = "chi_fit_cos_lr_bb_ssAnddon",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	geom.chi,
	acc_site.HBChemType AS acc_chem_type,
	don_site.HBChemType AS don_chem_type,
	don_ss.dssp AS don_ss, acc_ss.dssp AS acc_ss,
  hbond.donRank AS donRank,
  hbond.accRank AS accRank
FROM
	hbond_geom_coords AS geom,
	hbonds AS hbond,
	hbond_sites AS don_site,
	hbond_sites AS acc_site,
	residue_secondary_structure AS don_ss,
OB	residue_secondary_structure AS acc_ss
WHERE
	hbond.struct_id = geom.struct_id AND hbond.hbond_id   = geom.hbond_id AND
	hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
	hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
	don_ss.struct_id = don_site.struct_id AND don_ss.resNum = don_site.resNum AND
	acc_ss.struct_id = acc_site.struct_id AND acc_ss.resNum = acc_site.resNum AND
	acc_site.HBChemType  = 'hbacc_PBA' AND
	don_site.HBChemType != 'hbdon_PBA' AND
	ABS(don_site.resNum - acc_site.resNum) > 5;"

f <- query_sample_sources(sample_sources, sele)
f$chi <- (f$chi*180/pi)

f$acc_ss <- factor(f$acc_ss,
	levels = c(
		"H", "E", "T",
		"G", "B", "S",
		"I", " "),
	labels = c(
		'H: alpha-Helix',
		'E: beta-Sheet',
		'T: HB Turn',
		'G: 3/10 Helix',
		'B: isolated beta-Bridge',
		'S: Bend',
		'I: pi-Helix',
		'Loop or Irregular'))

plot_id = "chi_long_range_bbloopacc_by_dontype"

dens <- estimate_density_1d_wrap(
	f, c("sample_source", "acc_ss", "don_chem_type"), "chi")
ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	facet_grid(don_chem_type ~ acc_ss) +
	ggtitle("Hydrogen Bonds CHI Angle for Backbone Acceptors and Sidechain Donors with Sequence Separation > 5\nBy Acceptor DSSP Secondary Structure Type and Donor Type") +
	scale_x_continuous('Acceptor Base -- Acceptor Torsion (degrees)', breaks=c(90,270)) +
	scale_y_continuous('Feature Density') +
	theme(legend.position=c(.7, .35)) +
	theme(legend.justification=c("left", "top"))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)



})) # end FeaturesAnalysis
