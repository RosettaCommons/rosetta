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
id = "ASN_conditional_on_short_range_backbone_hbond_formation",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueFeatures", "ProteinResidueConformationFeatures", "PdbDataFeatures", "HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <-"
CREATE TEMPORARY TABLE asn_sr_bb_hb AS
SELECT DISTINCT
	res.struct_id,
	res.resNum,
	res.name3,
	CASE don_bb.resNum - acc_asn.resNum
		WHEN 1 THEN 1 ELSE -1 END AS sr_bb_hb
FROM
	hbonds AS hb,
	hbond_sites AS acc_asn,
	hbond_sites AS don_bb,
	residues AS res
WHERE

	acc_asn.site_id = hb.acc_id AND
	acc_asn.struct_id = hb.struct_id AND
	acc_asn.HBChemType = 'hbacc_CXA' AND
	res.struct_id = acc_asn.struct_id AND
	res.resNum = acc_asn.resNum AND
	res.name3 = 'ASN' AND
	don_bb.site_id = hb.don_id AND
	don_bb.struct_id = hb.struct_id AND
	don_bb.HBChemType = 'hbdon_PBA' AND
	(don_bb.resNum - 1 = acc_asn.resNum OR
	don_bb.resNum + 1 = acc_asn.resNum);

SELECT
	res.name3 AS res_type,
	res_dofs.chi1, res_dofs.chi2,
	CASE sr_asn.sr_bb_hb
		WHEN 1 THEN 1 WHEN -1 THEN -1 ELSE 'other' END AS sr_bb_hb

FROM
	residues AS res,
	residue_pdb_confidence AS res_conf,
	protein_residue_conformation AS res_dofs
	LEFT JOIN asn_sr_bb_hb AS sr_asn ON
		res.struct_id = sr_asn.struct_id AND
		res.resNum = sr_asn.resNum
WHERE
	res.name3 = 'ASN' AND
	res_conf.struct_id = res.struct_id AND res_conf.residue_number = res.resNum AND
	res_conf.max_temperature < 30 AND
	res_dofs.struct_id = res.struct_id AND res_dofs.seqpos == res.resNum;"

f <- query_sample_sources(sample_sources, sele)

m_f <- melt(f, measure.vars=c("chi1", "chi2"), variable_name="chi_angle")

dens <- estimate_density_1d_wrap(
	m_f, c("sample_source", "chi_angle", "sr_bb_hb"), "value", xlim=c(-180, 180), adjust=.35)


plot_id <-"rotamer_ASN_conditional_on_short_range_backbone_hbond_formation"
p <- ggplot(data=dens) +
	theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	facet_grid(sr_bb_hb ~ chi_angle) +
	ggtitle(("ASN chi angles by presense of +1/-1 SC-BB hbonds; BFact < 30", sep="")) +
	scale_x_continuous("Dihedral Angle") +
	scale_y_continuous("Feature Density") 
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
