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
id = "HTFY_by_dssp_lower_chi_tight_kernel",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueFeatures", "ProteinResidueConformationFeatures", "PdbDataFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <-"
SELECT
	res.name3 AS res_type,
	res_dofs.chi2 AS chi_angle,
	dssp_code.code AS dssp_code,
	dssp_code.label AS dssp_label,
	CASE
-- figure out which is gauche+, gauche- and trans
		WHEN 0 < res_dofs.chi1 AND res_dofs.chi1 < 120 THEN 'chi1 60'
		WHEN -120 < res_dofs.chi1 AND res_dofs.chi1 < 0 THEN 'chi1 -60'
	  ELSE 'chi1 180' END AS chi1_bin
FROM
	residues AS res,
	residue_pdb_confidence AS res_conf,
	residue_secondary_structure AS res_ss,
	dssp_codes AS dssp_code,
	protein_residue_conformation AS res_dofs
WHERE
	res_conf.struct_id = res.struct_id AND res_conf.residue_number = res.resNum AND
	res_conf.max_temperature < 30 AND
	res_ss.struct_id = res.struct_id AND res_ss.resNum == res.resNum AND
	dssp_code.code = res_ss.dssp AND
	res_dofs.struct_id = res.struct_id AND res_dofs.seqpos == res.resNum AND
	(res.name3 == 'HIS' OR res.name3 == 'TRP' OR
	res.name3 = 'PHE' OR res.name3 = 'TYR');"

f <- query_sample_sources(sample_sources, sele)

f$res_type <- factor(f$res_type,
	levels=c("HIS", "TRP", "PHE", "TYR"),
 	labels=c("HIS", "TRP", "PHE", "TYR"))

dens <- estimate_density_1d_wrap(
	f, c("sample_source", "res_type", "chi1_bin", "dssp_code", "dssp_label"), "chi_angle", xlim=c(-180, 180), adjust=.2)


d_ply(dens, .(dssp_code, dssp_label), function(sub_dens){
	dssp_code <- as.character(sub_dens$dssp_code[1])
	dssp_label <- as.character(sub_dens$dssp_label[1])
	plot_id <- paste("rotamer_HTFY_", dssp_code, "_chi2_by_chi1_tight_kernel", sep="")
	p <- ggplot(data=sub_dens) +
		theme_bw() +
		geom_line(aes(x=x, y=y, colour=sample_source)) +
		geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
		facet_grid(chi1_bin ~ res_type) +
#		ggtitle(("HIS/TRP/PHE/TYR ", dssp_label, " chi2 by chi1, BFact < 30", sep="")) +
		scale_x_continuous("Dihedral Angle") +
		scale_y_continuous("FeatureDensity", limits=c(0,.07))
	if(nrow(sample_sources) <= 3){
		p <- p + theme(legend.position="bottom", legend.direction="horizontal")
	}

	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
})

})) # end FeaturesAnalysis
