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
id = "bifurcated_at_donor",
author = "Matthew O'Meara",
brief_description = "",
long_description = "The rank of a hydrogen bond at donor site or acceptor site is rank of the relative Rosetta HBond energy of the hydrogen bond at the site.",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
CREATE TEMPORARY TABLE hbs AS SELECT
	don.struct_id AS struct_id, don.site_id AS don_site_id,
	geom.AHDist AS AHdist,
	acc.HBChemType AS acc_chem_type, don.HBChemType AS don_chem_type,
	hbond.donRank AS don_rank,
	don_ss.dssp AS don_ss
FROM
	hbonds AS hbond,
	hbond_sites AS acc,	hbond_sites AS don,
	hbond_geom_coords AS geom,
	hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb,
	residue_secondary_structure AS don_ss
WHERE
	hbond.struct_id = acc.struct_id AND hbond.acc_id = acc.site_id AND
	hbond.struct_id = don.struct_id AND hbond.don_id = don.site_id AND
	hbond.struct_id = geom.struct_id AND hbond.hbond_id = geom.hbond_id AND
	don.struct_id = don_pdb.struct_id AND don.site_id = don_pdb.site_id AND
	acc.struct_id = acc_pdb.struct_id AND acc.site_id = acc_pdb.site_id AND
	don_pdb.heavy_atom_temperature < 30 AND acc_pdb.heavy_atom_temperature < 30 AND
	don_ss.struct_id = don.struct_id AND don_ss.resNum = don.resNum;

CREATE INDEX hbs_struct_id_don_site_id ON hbs(struct_id, don_site_id);

SELECT
	hb1.AHdist AS AHdist1, hb2.AHdist AS AHdist2,
	hb1.acc_chem_type AS acc_chem_type1, hb2.acc_chem_type AS acc_chem_type2,
	hb1.don_chem_type AS don_chem_type,
	hb1.don_ss
FROM
	hbs AS hb1, hbs AS hb2
WHERE
	hb1.struct_id = hb2.struct_id AND hb1.don_site_id = hb2.don_site_id AND
	hb1.don_rank > hb2.don_rank;"

f <- query_sample_sources(sample_sources, sele)

# This is deprecated please use the hbond_chem_types table for the lables instead
# Order the plots better and give more descriptive labels
f$acc_chem_type_name1 <- factor(f$acc_chem_type1,
	levels = c("hbacc_IMD", "hbacc_IME", "hbacc_AHX", "hbacc_HXL",
		"hbacc_CXA", "hbacc_CXL", "hbacc_PBA"),
	labels = c("aIMD: h", "aIME: h", "aAHX: y", "aHXL: s,t",
		"aCXA: n,q", "aCXL: d,e", "aPBA: bb"))

# This is deprecated please use the hbond_chem_types table for the lables instead
# Order the plots better and give more descriptive labels
f$acc_chem_type_name2 <- factor(f$acc_chem_type2,
	levels = c("hbacc_IMD", "hbacc_IME", "hbacc_AHX", "hbacc_HXL",
		"hbacc_CXA", "hbacc_CXL", "hbacc_PBA"),
	labels = c("aIMD: h", "aIME: h", "aAHX: y", "aHXL: s,t",
		"aCXA: n,q", "aCXL: d,e", "aPBA: bb"))

# This is deprecated please use the hbond_chem_types table for the lables instead
# Order the plots better and give more descriptive labels
f$don_chem_type_name <- factor(f$don_chem_type,
	levels = c("hbdon_IMD", "hbdon_IME", "hbdon_GDE", "hbdon_GDH",
		"hbdon_AHX", "hbdon_HXL", "hbdon_IND", "hbdon_AMO", "hbdon_CXA", "hbdon_PBA"),
	labels = c("dIMD: h", "dIME: h", "dGDE: r", "dGDH: r",
		"dAHX: y", "dHXL: s,t", "dIND: w", "dAMO: k", "dCXA: n,q", "dPBA: bb"))

f <- na.omit(f, method="r")
f$don_ss_name <- factor(f$don_ss,
	levels = c("H", "E", "T", "G", "B", "S", "I", "C"),
	labels = c('H: a-Helix', 'E: b-Sheet', 'T: HB Turn', 'G: 3/10 Helix',
		'B: b-Bridge', 'S: Bend', 'I: pi-Helix','C: Irregular'))


plot_parts <- list(
	theme_bw(),
	geom_point(aes(x=AHdist1, y=AHdist2), size=.4),
#	stat_density2d(aes(x=AHdist1, y=AHdist2), size=.2),
	geom_indicator(aes(indicator=counts)), # TODO: fix colour
	facet_grid(acc_chem_type_name1 ~ acc_chem_type_name2),
	scale_x_continuous(
		expression(
			paste('Acceptor -- Hydrogen Distance (', ring(A), ') Worse Energy HBond')),
		limits=c(1.5, 3), breaks=c(1.5, 2.0, 2.5)),
	scale_y_continuous(
		expression(
			paste('Acceptor -- Hydrogen Distance (', ring(A), ') Better Energy HBond')),
		limits=c(1.5, 3), breaks=c(1.5, 2.0, 2.5)))

d_ply(f, .(sample_source), function(sub_f){
	ss_id <- sub_f$sample_source[1]

	plot_id <- "hbond_bifurcated_at_donor_AHdist"
	sub_f <- ddply(sub_f, .(acc_chem_type_name1, acc_chem_type_name2),
		transform, counts = length(sample_source))
	ggplot(sub_f) + plot_parts +
		ggtitle(paste("HBond Pairs Bifurcated at Donor; bFact < 30  ss_id:",ss_id))
	save_plots(self, 
		plot_id,
		sample_sources[sample_sources$sample_source == ss_id,],
		output_dir,
		output_formats)

	d_ply(sub_f, .(don_chem_type), function(sub_sub_f){
		don_chem_type <- sub_sub_f$don_chem_type[1]
		don_chem_type_name <- sub_sub_f$don_chem_type_name[1]

		plot_id <- paste("hbond_bifurcated_at_",don_chem_type, "_AHdist", sep="")
		sub_sub_f <- ddply(sub_sub_f, .(acc_chem_type_name1, acc_chem_type_name2),
			transform, counts = length(sample_source))
		ggplot(sub_sub_f) + plot_parts +
			ggtitle(paste(
				"HBond Pairs Bifurcated at ", don_chem_type_name, ";",
				" bFact < 30    ss_id:",ss_id))
		save_plots(self, 
			plot_id,
			sample_sources[sample_sources$sample_source == ss_id,],
			output_dir, output_formats)
	})

	plot_id <- paste("hbond_backbone_bifurcated_AHdist", sep="")
	sub_sub_f <- subset(sub_f,
		don_chem_type = "hbdon_PBA",
		acc_chem_type_name1 = "hbacc_PBA", acc_chem_type_name2 = "hbacc_PBA")
	sub_sub_f <- ddply(sub_sub_f, .(don_ss),
		transform, counts = length(sample_source))
	ggplot(sub_sub_f) + plot_parts + facet_wrap(~don_ss) +
		ggtitle(paste(
			"HBond Backbone Donor Bifurcated to 2 Backbone acceptors;",
			" bFact < 30    ss_id:",ss_id))
	save_plots(self, 
		plot_id,
		sample_sources[sample_sources$sample_source == ss_id,],
		output_dir, output_formats)

})


})) # end FeaturesAnalysis
