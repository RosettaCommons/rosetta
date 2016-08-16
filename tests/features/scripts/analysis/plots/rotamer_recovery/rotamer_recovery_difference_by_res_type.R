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
id = "rotamer_recovery_by_res_type",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("StructureFeatures", "ResidueFeatures", "RotamerRecoveryFeatures", "PdbDataFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

	
sele <-"
SELECT
	new_rr.divergence - ref_rr.divergence AS rel_div,
	res.name3 AS res_type,
	CASE bur.sasa_r140 WHEN 0 THEN 'buried' ELSE 'exposed' END AS burial
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
	res.res_type != 'ALA' AND res.res_type != 'GLY';"

f <- query_sample_sources_against_ref(sample_sources, sele)

f$res_type <- factor(f$res_type,
	levels = c(
		"SER", "VAL", "HIS", "ASN", "GLN",
		"THR", "ILE", "PHE", "ASP", "GLU",
		"CYS", "LEU", "TYR", "ARG", "LYS",
		"PRO", "MET", "TRP"))

f <- na.omit(f, method="r")

f <- ddply(f, .(sample_source, res_type, burial),
	transform, mean = round(mean(rel_div), 4))

f <- ddply(f, .(sample_source, res_type, burial),
	transform, counts = length(sample_source))


plot_id <- "rotamer_recovery_differences_by_res_type_and_burial"
d_ply(f, .(sample_source), function(sub_f) {
	ss_id <- sub_f$sample_source[1]

	p <- ggplot(data=sub_f) + theme_bw() +
		stat_bin(aes(x=rel_div, y=log(..count.. + 1), colour=burial),
			breaks=c(-4.5, -3.5, -2.5, -1.5, -.5, .5, 1.5, 2.5, 3.5, 4.5), geom="line") +
		geom_indicator(aes(indicator=counts, colour=burial, group=burial)) +
		geom_indicator(aes(indicator=mean, color=burial, group=burial), xpos="left") +
		facet_wrap( ~ res_type ) +
		ggtitle("Rotamer Recovery by Residue Type and Burial, B-Factor < 20") +
		labs(x="New Rotamer Recovery - Ref Recovery Score", y="log(counts + 1)") +
		theme(legend.position=c(.8, .25)) +
		theme(legend.justification=c("left", "top"))
	save_plots(self, plot_id, sample_sources[sample_sources$sample_source == ss_id,], output_dir, output_formats)
})

})) # end FeaturesAnalysis
