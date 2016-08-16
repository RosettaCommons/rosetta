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
id = "DENQ_by_secondary_structure",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueFeatures", "ProteinResidueConformationFeatures", "PdbDataFeatures", "ResidueSecondaryStructureFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <-"
SELECT
	res.name3 AS res_type,
	res_dofs.chi2 AS last_chi_angle,
	dssp_code.label AS dssp
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
	(res.name3 == 'ASN' OR res.name3 == 'ASP')
UNION
SELECT
	res.name3 AS res_type,
	res_dofs.chi3 AS last_chi_angle,
	dssp_code.label AS dssp
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
	(res.name3 == 'GLN' OR res.name3 == 'GLU');"

f <- query_sample_sources(sample_sources, sele)

f$res_type <- factor(f$res_type,
	levels=c("ASN", "ASP", "GLN", "GLU"),
 	labels=c("ASN", "ASP", "GLN", "GLU"))

dens <- estimate_density_1d(
	f, c("sample_source", "res_type", "dssp"), "last_chi_angle", xlim=c(-180, 180))


plot_id <-"rotamer_denq_by_secondary_structure"
p <- ggplot(data=dens) +
	theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	facet_grid(dssp ~ res_type) +
#	ggtitle(("ASN/ASP chi2 and GLN/GLU chi3 by Secondary Structure BFact < 30", sep="")) +
	scale_x_continuous("Dihedral Angle") +
	scale_y_continuous("FeatureDensity", limits=c(0,.04))
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
