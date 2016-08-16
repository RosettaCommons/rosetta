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
id = "bidentite_hbond_asymmetry",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("StructureFeatures", "ResidueFeatures", "SaltBridgeFeatures", "PdbDataFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <- "
DROP TABLE IF EXISTS CXL_ARG_salt_bridges;
CREATE TABLE IF NOT EXISTS CXL_ARG_salt_bridges AS SELECT
	sb.*,
	sb_alt.struct_id IS NOT NULL AS bidentite
FROM
	(SELECT
		sb.struct_id,
		sb.don_resNum AS don_residue_number,
		acc.resNum AS acc_residue_number,
		sb.acc_id AS acc_id,
		sb.psi * 180/3.14159 AS psi,
		acc_env.sasa_r140 AS sasa
	FROM
		salt_bridges AS sb,
		hbond_sites AS acc,
		residues AS don,
		hbond_site_environment AS acc_env
	WHERE
		acc.struct_id = sb.struct_id AND
		acc.site_id = sb.acc_id AND
		acc.HBChemType = 'hbacc_CXL' AND

		don.struct_id = sb.struct_id AND
		don.resNum = sb.don_resNum AND
		don.name3 = 'ARG' AND

		ABS(acc.resNum - don.resNum) > 4 AND

		acc_env.struct_id = acc.struct_id AND
		acc_env.site_id = acc.site_id
	) AS sb LEFT JOIN
	(SELECT
		sb_alt.struct_id,
		sb_alt.don_resNum AS don_residue_number,
		sb_alt.acc_id AS acc_id,
		acc.resNum AS acc_residue_number
	FROM
		salt_bridges AS sb_alt,
		hbond_sites AS acc
	WHERE
		acc.struct_id = sb_alt.struct_id AND
		acc.site_id = sb_alt.acc_id) AS sb_alt
	ON
		(sb.struct_id = sb_alt.struct_id AND
		sb.don_residue_number = sb_alt.don_residue_number AND
		sb.acc_residue_number = sb_alt.acc_residue_number AND
		sb.acc_id != sb_alt.acc_id);

CREATE INDEX IF NOT EXISTS CXL_ARG_salt_bridges_struct_id_don_residue_number_acc_residue_number ON
	CXL_ARG_salt_bridges (struct_id, don_residue_number, acc_residue_number);

CREATE INDEX IF NOT EXISTS hbonds_struct_id_acc_id ON
	hbonds (struct_id, acc_id);

SELECT
	hb_coords1.AHdist AS AHdist1,
	hb_coords2.AHdist AS AHdist2,
	sb.sasa AS acc_sasa,
	sb.bidentite
FROM
	CXL_ARG_salt_bridges AS sb,
	hbonds AS hb1,
	hbonds AS hb2,
	hbond_geom_coords AS hb_coords1,
	hbond_geom_coords AS hb_coords2,
	hbond_sites AS don1,
	hbond_sites AS don2
WHERE
	hb1.struct_id = sb.struct_id AND
	hb2.struct_id = sb.struct_id AND
	hb_coords1.struct_id = sb.struct_id AND
	hb_coords2.struct_id = sb.struct_id AND
	don1.struct_id = sb.struct_id AND
	don2.struct_id = sb.struct_id AND

	hb1.hbond_id != hb2.hbond_id AND
	hb1.acc_id = sb.acc_id AND
	hb2.acc_id = sb.acc_id AND

	don1.resNum = sb.don_residue_number AND
	hb1.don_id = don1.site_id AND

	don2.resNum = sb.don_residue_number AND
	hb2.don_id = don2.site_id AND

	hb_coords1.hbond_id = hb1.hbond_id AND
	hb_coords2.hbond_id = hb2.hbond_id;"

f <- query_sample_sources(sample_sources, sele)
f$close_AHdist <- pmin(f$AHdist1, f$AHdist2)
f$far_AHdist <- pmax(f$AHdist1, f$AHdist2)

f$asymmetry <- abs(f$AHdist1 - f$AHdist2)
f$is_symmetric <- as.factor(ifelse(f$asymmetry == 0, "symmetric", "asymmetric"))

f$is_buried <- as.factor(ifelse(f$acc_sasa == 0, "buried", "exposed"))

f$is_bidentite <- as.factor(ifelse(f$bidentite == 1, "bidentite", "not bidentite"))


counts <- ddply(f, .(is_bidentite, is_symmetric, sample_source), function(df){
	data.frame(count = nrow(df))})
plot_id <- "ASP_GLU_to_ARG_bifurcated_symmetry_scatter"
p <- ggplot(data=f) + theme_bw() +
	geom_point(aes(x=close_AHdist, y=far_AHdist, colour=is_symmetric), size=.7) +
	geom_indicator(data=counts, aes(colour=is_symmetric, group=is_symmetric, indicator=count), ypos=.08) +
	ggtitle("ASP/GLU -> ARG Bifurcated Salt Bridge Asymmetry; SeqSep > 4") +
	scale_x_continuous("Close AHdist") +
	scale_y_continuous("Far AHdist") +
	facet_grid(is_bidentite ~sample_source)
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
