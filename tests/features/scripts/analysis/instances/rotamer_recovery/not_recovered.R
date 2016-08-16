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
id = "not_recovered",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("StructureFeatures", "ResidueFeatures", "RotaermRecoveryFeatures", "PdbDataFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	ref_struct.tag, '' AS chain, res.resNum,
	'CA' AS CA, 'C' AS C, 'CB' AS CB,
	res.name3 AS res_type,
	new_rr.divergence - ref_rr.divergence AS rel_div
FROM
	ref.structures AS ref_struct, new.structures AS new_struct,
	ref.rotamer_recovery AS ref_rr, new.rotamer_recovery AS new_rr,
	ref.residues AS res,
	ref.residue_pdb_confidence AS b,
	ref.residue_burial AS bur
WHERE
	ref_struct.tag = new_struct.tag AND
	res.struct_id = ref_struct.struct_id AND
	ref_rr.struct_id = ref_struct.struct_id AND ref_rr.resNum = res.resNum AND
	new_rr.struct_id = new_struct.struct_id AND new_rr.resNum = res.resNum AND
	b.struct_id = ref_struct.struct_id AND b.residue_number = res.resNum AND
	b.max_temperature < 20 AND
	bur.struct_id = res.struct_id AND bur.resNum = res.resNum AND
	bur.sasa_r140 = 0 AND
	res.res_type = 'ARG';"

f <-  query_sample_sources_against_ref(sample_sources, sele)

plot_id <- "rotamer_recovery_LYS_relative_divergence"
ref_ss_id <- f$ref_sample_source[1]
p <- ggplot(data=f) + theme_bw() +
	geom_histogram(aes(x=rel_div, fill=new_sample_source)) +
	ggtitle(paste("Lysine Relative Rotamer Recovery against ", ref_ss_id, ", 0 Sasa and B-Factor < 20", sep="")) +
	labs(x="New Recovery - Ref Recovery", y="log(FeatureDensity + 1)")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

n_examples <- 15

f <- f[order(f$rel_div, decreasing=T),]
f$id <- 1:nrow(f)
g <- melt(f[f$id <= n_examples,],
	id.vars=c("sample_source", "tag", "id", "chain", "resNum"),
	measure.vars=c("CA", "C", "CB"),
	variable_name = "atom")

instances_id <- "rotamer_recovery_LYS_max_diff_divergence"
prepare_feature_instances(instances_id, sample_sources, g, output_dir)

})) # end FeaturesAnalysis
