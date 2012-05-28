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
id = "residue_pair_distances_vs_neighbors",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueFeatures", "ResidueTypeFeatures", "ResiduePairFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
DROP TABLE IF EXISTS pair_term_res_types;
CREATE TABLE pair_term_res_types (
	res_type TEXT,
	PRIMARY KEY (res_type));

INSERT INTO pair_term_res_types
SELECT
	res_type.name3 AS res_type
FROM
	residue_type AS res_type,
	residue_type_property AS polar_or_aromatic,
	residue_type_property AS protein
WHERE
	res_type.residue_type_set_name == 'fa_standard' AND

	polar_or_aromatic.residue_type_set_name = 'fa_standard' AND
	polar_or_aromatic.residue_type_name = res_type.name AND
	(polar_or_aromatic.property == 'POLAR' OR
	polar_or_aromatic.property == 'AROMATIC') AND

	protein.residue_type_set_name == 'fa_standard' AND
	protein.residue_type_name == res_type.name AND
	protein.property == 'PROTEIN' AND
  -- I'm not sure how to restrict the residue type so we don't have to enumerate. Perhaps we need to add a field to the residue type set about begin a canonical amino acid?
	res_type.name3 in ('ARG', 'ASN', 'ASP', 'GLN', 'HIS', 'LYS', 'PHE', 'SER', 'THR', 'TRP', 'TYR')
GROUP BY
	res_type;

SELECT
	CASE
  	WHEN rp.actcoord_dist < 4.5 THEN 4.5
  	WHEN rp.actcoord_dist < 6   THEN 6
  	WHEN rp.actcoord_dist < 7.5 THEN 7.5
    ELSE NULL END AS r12_bin,
	CASE
		WHEN rp.res1_10A_neighbors <= 16 THEN 'exposed'
		ELSE 'buried' END AS res1_burial,
	CASE
		WHEN rp.res2_10A_neighbors <= 16 THEN 'exposed'
		ELSE 'buried' END AS res2_burial,
  res1.name3 AS res1_type,
  res2.name3 AS res2_type,
  count(*) AS count
FROM
  residue_pairs AS rp,
  residues AS res1,
  residues AS res2,
	residue_pdb_confidence AS res1_conf,
	residue_pdb_confidence AS res2_conf,
	pair_term_res_types AS res1_valid_res_type,
	pair_term_res_types AS res2_valid_res_type
WHERE
  rp.actcoord_dist < 7.5 AND
  res1.struct_id = rp.struct_id AND res1.resNum = rp.resNum1 AND
  res2.struct_id = rp.struct_id AND res2.resNum = rp.resNum2 AND
  res1_conf.struct_id = rp.struct_id AND res1_conf.residue_number = rp.resNum1 AND
  res2_conf.struct_id = rp.struct_id AND res2_conf.residue_number = rp.resNum2 AND
	res1_conf.max_temperature < 30 AND
	res2_conf.max_temperature < 30 AND
	res1.name3 = res1_valid_res_type.res_type AND
	res2.name3 = res2_valid_res_type.res_type
GROUP BY
	r12_bin,
  res1_burial,
	res2_burial,
	res1_type,
	res2_type;"

f <-  query_sample_sources(sample_sources, sele)

f$r12_bin <- factor(f$r12_bin)

plot_id <- "residue_pair_distances_buried"
p <- ggplot(data=f[f$res1_burial == 'buried' & f$res2_burial=='buried',]) +
	theme_bw() +
	geom_line(aes(x=r12_bin, y=log(count), colour=sample_source, group=sample_source), size=.2) +
	geom_point(aes(x=r12_bin, y=log(count), colour=sample_source), size=2) +
	facet_grid(res1_type ~ res2_type) +
	opts(title = "Residue Pair Distances When Both Are Buried; B-Fact < 30") +
	scale_x_discrete("Distance between Action Coordinates") +
	scale_y_continuous("Log(Number of Neighbors)")
if(nrow(sample_sources) <= 3){
	p <- p + opts(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
