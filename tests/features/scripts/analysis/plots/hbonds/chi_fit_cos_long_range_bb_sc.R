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
id = "chi_fit_cos_long_range_bb_sc",
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
	residue_secondary_structure AS acc_ss
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

# Change chi from radians to degrees
f$chi <- (f$chi*180/pi)

# The coding for loop and irregular is " ", change this to "C"
levels(f$acc_ss)[levels(f$acc_ss) == " "] <- "C"

# Add descriptive names to the DSSP codes
# Reorder factor to lign up the facet plots nicely
f$acc_ss_name <- factor(f$acc_ss,
	levels = c(
		"H", "E", "T",
		"G", "B", "S",
		"I", "C"),
	labels = c(
		'H: alpha-Helix',
		'E: beta-Sheet',
		'T: HB Turn',
		'G: 3/10 Helix',
		'B: isolated beta-Bridge',
		'S: Bend',
		'I: pi-Helix',
		'C: Loop or Irregular'))

# Make Rank a categorical variable
f$accRank <- factor(f$accRank)
f$donRank <- factor(f$donRank)

plot_id = "chi_long_range_acc_ss"
dens <- estimate_density_1d_wrap(
	f, c("sample_source", "acc_ss_name"), "chi")
ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	facet_wrap(~acc_ss_name, ncol=3) +
	ggtitle("Hydrogen Bonds CHI Angle for Backbone Acceptors and Sidechain Donors with Sequence Separation > 5\nBy Acceptor DSSP Secondary Structure Type") +
	scale_x_continuous('Acceptor Base -- Acceptor Torsion (degrees)', breaks=c(90,270)) +
	scale_y_continuous('Feature Density') +
	theme(legend.position=c(.7, .35)) +
	theme(legend.justification=c("left", "top"))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


plot_id = "chi_long_range_acc_ss_by_rank"
l_ply(levels(f$sample_source), function(ss){
	sub_f <- f[f$sample_source == ss,]
	dens <- estimate_density_1d_wrap(
		sub_f, c("acc_ss_name", "accRank"), "chi")
	ggplot(data=dens) + theme_bw() +
		geom_line(aes(x=x, y=y, colour=accRank)) +
		geom_indicator(aes(indicator=counts, colour=accRank, group=accRank)) +
		facet_wrap(~acc_ss_name, ncol=3) +
		ggtitle(paste("Hydrogen Bonds CHI Angle for Backbone Acceptors and Sidechain Donors with Sequence Separation > 5\nBy Acceptor DSSP Secondary Structure Type   Sample Source: ", ss, sep="")) +
		scale_x_continuous('Acceptor Base -- Acceptor Torsion (degrees)', breaks=c(90,270)) +
		scale_y_continuous('Feature Density') +
		theme(legend.position=c(.7, .35)) +
		theme(legend.justification=c("left", "top"))
	save_plots(self, plot_id, subset(sample_sources, ss %in% sample_source), output_dir, output_formats)

})

dens <- estimate_density_1d_wrap(
	f, c("sample_source", "acc_ss", "acc_ss_name", "don_chem_type"), "chi")
ylim <- range(dens$y)

l_ply(levels(dens$acc_ss), function(ss){
	plot_id = paste("chi_long_range_by_don_chem_type_acc_ss_",ss, sep="")
	sub_dens <- dens[ss == dens$acc_ss,]
	ss_name <- sub_dens[1,"acc_ss_name"]
	ggplot(data=sub_dens) + theme_bw() +
		geom_line(aes(x=x, y=y, colour=sample_source)) +
		geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
		facet_wrap(~don_chem_type, ncol=3) +
		ggtitle(paste("Hydrogen Bonds CHI Angle for Backbone Acceptors and Sidechain Donors with Sequence Separation > 5\nBy Acceptor DSSP Secondary Structure Type: ", ss_name, sep="")) +
		scale_x_continuous('Acceptor Base -- Acceptor Torsion (degrees)', breaks=c(90,270)) +
		scale_y_continuous('Feature Density', limits=ylim) +
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})


})) # end FeaturesAnalysis
