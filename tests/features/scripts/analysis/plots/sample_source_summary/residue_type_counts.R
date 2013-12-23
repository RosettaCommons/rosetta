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
id = "residue_type_counts",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueFeatures", "ResidueTypesFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <- "
SELECT
	res.name3 AS res_type,
	CASE
		WHEN res_protein.property IS NULL THEN 'not_protein'
		ELSE 'protein' END AS is_protein,
  count(res.res_type) AS count
FROM
  residues AS res,
	residue_pdb_confidence AS res_conf
	LEFT JOIN residue_type_property AS res_protein ON
		res_protein.residue_type_set_name = 'fa_standard' AND
		res_protein.residue_type_name = res.res_type AND
		res_protein.property = 'PROTEIN'
WHERE
	res_conf.struct_id = res.struct_id AND
	res_conf.residue_number = res.resNum
GROUP BY
	res.name3;"

f <-  query_sample_sources(sample_sources, sele)

plot_id <- "residue_type_counts"
p <- ggplot(data=f) + theme_bw() +
	geom_bar(aes(x=res_type, y=count, fill=sample_source), stat="identity", position="dodge") +
	coord_flip() +
	ggtitle("Residue Type Counts (B-Fact < 30)") +
	labs(x = "Residue Type", y = "Count")

save_plots(self, plot_id, sample_sources, output_dir, output_formats)




table_id <- "residue_types_counts"
f_c <- cast(f, res_type ~ sample_source, value="count")
save_tables(
	self,
	f_c,
	table_id,
	sample_sources, output_dir, output_formats,
	caption="Residue Type Counts", caption.placement="top")

table_id <- "protein_residue_types_counts"
f_c <- cast(f[as.character(f$is_protein) == 'protein',], res_type ~ sample_source, value="count")
save_tables(
	self,
	f_c,
	table_id,
	sample_sources, output_dir, output_formats,
	caption="Protein Residue Type Counts", caption.placement="top")

table_id <- "non_protein_residue_types_counts"
f_c <- cast(f[as.character(f$is_protein) != 'protein',], res_type ~ sample_source, value="count")
save_tables(
	self,
	f_c,
	table_id,
	sample_sources, output_dir, output_formats,
	caption="Non-Protein Residue Type Counts", caption.placement="top")



})) # end FeaturesAnalysis
