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
id = "residue_pair_distances_vs_neighbors",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueFeatures", "ResidueTypeFeatures", "ResiduePairFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	rp.actcoord_dist,
	CASE
		WHEN rp.res1_10A_neighbors <= 16 THEN 'exposed'
		ELSE 'buried' END AS res1_burial,
	CASE
		WHEN rp.res2_10A_neighbors <= 16 THEN 'exposed'
		ELSE 'buried' END AS res2_burial,
	CASE WHEN res1.name3 < res2.name3 THEN res1.name3 ELSE res2.name3 END AS res1_type,
	CASE WHEN res1.name3 < res2.name3 THEN res2.name3 ELSE res1.name3 END AS res2_type
FROM
	residue_pairs AS rp,
	residues AS res1,
	residues AS res2,
	residue_pdb_confidence AS res1_conf,
	residue_pdb_confidence AS res2_conf
WHERE
	rp.actcoord_dist > 2 AND -- only difulfide bonds are less than 2A
	rp.actcoord_dist < 7.5 AND
	res1.struct_id = rp.struct_id AND res1.resNum = rp.resNum1 AND
	res2.struct_id = rp.struct_id AND res2.resNum = rp.resNum2 AND
	res1_conf.struct_id = rp.struct_id AND res1_conf.residue_number = rp.resNum1 AND
	res2_conf.struct_id = rp.struct_id AND res2_conf.residue_number = rp.resNum2 AND
	res1_conf.max_temperature < 30 AND
	res2_conf.max_temperature < 30 AND
	res1.name3 in ('ALA', 'ARG', 'ASN', 'ASP', 'GLU', 'GLN', 'PRO', 'SER', 'TYR', 'HIS', 'ILE', 'LEU', 'VAL', 'GLY', 'LYS', 'PHE', 'THR', 'TRP', 'MET') AND
	res2.name3 in ('ALA', 'ARG', 'ASN', 'ASP', 'GLU', 'GLN', 'PRO', 'SER', 'TYR', 'HIS', 'ILE', 'LEU', 'VAL', 'GLY', 'LYS', 'PHE', 'THR', 'TRP', 'MET');"


f <-  query_sample_sources(sample_sources, sele)

dens <- estimate_density_1d(
	f[f$res1_burial == 'buried' & f$res2_burial=='buried',],
	c("sample_source", "res1_type", "res2_type"),
	"actcoord_dist", weight_fun=radial_3d_normalization)

plot_id <- "residue_pair_distances_buried"
p <- ggplot(data=dens) +
	theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source, group=sample_source)) +
	geom_indicator(aes(indicator=counts, group=sample_source, colour=sample_source), size=3) +
	facet_grid(res1_type ~ res2_type) +
	ggtitle("Residue Pair Distances When Both Are Buried; B-Fact < 30") +
	scale_x_continuous("Distance between Action Coordinates", range=c(2, 7.5)) +
	scale_y_continuous("Feature Density")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
