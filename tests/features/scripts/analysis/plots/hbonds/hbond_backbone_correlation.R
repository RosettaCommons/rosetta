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
id = "hbond_backbone_correlation",
author = "Matthew O'Meara",
brief_description = "This script addresses the question: Are the hydrogen bond geometries on either side of the backbone correlated?",
long_description = "
 Organization of joins to form query
rsd
	acc_site
		acc_hb
			acc_geom
		acc_env
	don_site
		don_hb
			don_geom
		don_env
",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <- "
CREATE INDEX IF NOT EXISTS hbond_sites_struct_id_resNum ON
	hbond_sites(struct_id, resNum);
CREATE INDEX IF NOT EXISTS hbonds_struct_id_acc_id ON
	hbonds(struct_id, acc_id);
CREATE INDEX IF NOT EXISTS hbonds_struct_id_don_id ON
	hbonds(struct_id, don_id);

CREATE TEMPORARY TABLE bb_hbs AS SELECT
	don.struct_id, don.resNum,
	acc.site_id AS acc_id, don.site_id AS don_id
FROM
	hbond_sites AS acc, hbond_sites AS don
WHERE
	acc.struct_id = don.struct_id AND acc.resNum = don.resNum AND
	acc.HBChemType = 'hbacc_PBA' AND don.HBChemType = 'hbdon_PBA';

SELECT
	acc_geo.AHdist AS acc_AHdist, don_geo.AHdist AS don_AHdist,
	acc_hb.energy AS acc_energy, don_hb.energy AS don_energy,
	res_ss.dssp AS dssp
FROM
	bb_hbs AS t,
	hbonds AS acc_hb, hbonds AS don_hb,
	hbond_geom_coords AS acc_geo, hbond_geom_coords AS don_geo,
	residue_secondary_structure AS res_ss
WHERE
	acc_hb.struct_id  = t.struct_id AND acc_hb.acc_id = t.acc_id AND
	don_hb.struct_id  = t.struct_id AND don_hb.don_id = t.don_id AND
	acc_geo.struct_id = t.struct_id AND acc_geo.hbond_id = acc_hb.hbond_id AND
	don_geo.struct_id = t.struct_id AND don_geo.hbond_id = don_hb.hbond_id AND
	res_ss.struct_id  = t.struct_id AND res_ss.resNum = t.resNum;"

f <- query_sample_sources(sample_sources, sele)

f <- ddply(f, c("sample_source", "dssp"), transform,
					 counts = length(sample_source),
					 correlation = cor(acc_AHdist, don_AHdist))

f$dssp <- factor(f$dssp,
	levels = c("H", "E", "T", "G", "B", "S", "I", " "),
	labels = c('H: a-Helix', 'E: b-Sheet', 'T: HB Turn', 'G: 3/10 Helix',
		'B: b-Bridge', 'S: Bend', 'I: pi-Helix','C: Irregular'))

f <- na.omit(f, method="r")

d_ply(f, .(sample_source), function(sub_f){
	ss_id = sub_f[1,"sample_source"]
	plot_id <- "hbond_backbone_correlation"
	ggplot(data=sub_f) + theme_bw() +
		geom_point(aes(acc_AHdist, don_AHdist), size=.3) +
		stat_density2d(aes(x=acc_AHdist, y=don_AHdist), size=.2) +
		geom_indicator(aes(indicator=counts)) +
		geom_indicator(aes(indicator=paste("cor:",round(correlation, 3)), xpos=.3)) +
		facet_wrap(~dssp) +
		coord_equal(ratio=1) +
		ggtitle(paste("HBond Backbone Amide Correlation by DSSP of Residue\nSample Source:", ss_id)) +
		labs(x=expression(paste('Backbone is Acceptor: Acceptor -- Proton Distance (', ring(A), ')')),
			y=expression(paste('Backbone is Donor: Acceptor -- Proton Distance (', ring(A), ')')))

	save_plots(self, 
		plot_id, sample_sources[sample_sources$sample_source == ss_id,],
		output_dir, output_formats)
})


})) # end FeaturesAnalysis
