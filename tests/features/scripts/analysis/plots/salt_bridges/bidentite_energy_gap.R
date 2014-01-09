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
id = "bidentite_energy_gap",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("StructureFeatures", "ResidueFeatures", "SaltBridgeFeatures", "PdbDataFeatures", "ResidueScoreFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

problem_region <- data.frame(xmin=-119, xmax=-116)


sele <-paste("
CREATE TABLE IF NOT EXISTS CXL_ARG_salt_bridges AS SELECT
	struct.tag AS tag,
	sb.struct_id AS struct_id,
	sb.don_resNum AS don_residue_number,
	sb.acc_id AS acc_id,
	acc.resNum AS acc_residue_number,
	sb.psi * 180/3.14159 AS psi
FROM
	structures AS struct,
	salt_bridges AS sb,
	residues AS don,
	hbond_sites AS acc
WHERE
	sb.struct_id = struct.struct_id AND
	don.struct_id = sb.struct_id AND
	don.resNum = sb.don_resNum AND
	don.name3 = 'ARG' AND
	acc.struct_id = sb.struct_id AND
	acc.site_id = sb.acc_id AND
	acc.HBChemType = 'hbacc_CXL' AND
	ABS(acc.resNum - don.resNum) > 4;

DROP TABLE IF EXISTS total_residue_scores;
CREATE TABLE total_residue_scores (
	struct_id INTEGER,
	resNum INTEGER,
	score_value);

DROP INDEX IF EXISTS total_residue_scores_struct_id_resNum;
CREATE INDEX total_residue_scores_struct_id_resNum ON total_residue_scores (
	struct_id, resNum);

CREATE INDEX IF NOT EXISTS residue_scores_2b_struct_id_resNum1 ON residue_scores_2b (struct_id, resNum1);
CREATE INDEX IF NOT EXISTS residue_scores_2b_struct_id_resNum2 ON residue_scores_2b (struct_id, resNum2);

INSERT INTO total_residue_scores SELECT
	e2b.struct_id, e2b.resNum2 AS resNum,
	SUM((e2b.score_value)/2) AS score_value
FROM
	CXL_ARG_salt_bridges AS sb,
	residue_scores_2b AS e2b
WHERE
	e2b.struct_id = sb.struct_id AND e2b.resNum2 = sb.acc_residue_number
GROUP BY
	e2b.struct_id, e2b.resNum2;

INSERT INTO total_residue_scores SELECT
	e2a.struct_id, e2a.resNum2 AS resNum,
	SUM((e2a.score_value)/2) AS score_value
FROM
	CXL_ARG_salt_bridges AS sb,
	residue_scores_2b AS e2a
WHERE
	e2a.struct_id = sb.struct_id AND e2a.resNum1 = sb.acc_residue_number
GROUP BY
	e2a.struct_id, e2a.resNum2;

INSERT INTO total_residue_scores SELECT
	e2b.struct_id, e2b.resNum2 AS resNum,
	SUM((e2b.score_value)/2) AS score_value
FROM
	CXL_ARG_salt_bridges AS sb,
	residue_scores_2b AS e2b
WHERE
	e2b.struct_id = sb.struct_id AND e2b.resNum2 = sb.don_residue_number
GROUP BY
	e2b.struct_id, e2b.resNum2;

INSERT INTO total_residue_scores SELECT
	e2a.struct_id, e2a.resNum2 AS resNum,
	SUM((e2a.score_value)/2) AS score_value
FROM
	CXL_ARG_salt_bridges AS sb,
	residue_scores_2b AS e2a
WHERE
	e2a.struct_id = sb.struct_id AND e2a.resNum1 = sb.don_residue_number
GROUP BY
	e2a.struct_id, e2a.resNum2;

DROP TABLE summed_total_residue_scores;
CREATE TABLE summed_total_residue_scores AS SELECT
	struct_id, resNum, SUM(score_value) AS score_value
FROM
	total_residue_scores
GROUP BY
	struct_id, resNum;



CREATE TABLE IF NOT EXISTS CXL_ARG_salt_bridge_energies AS SELECT
	sb.struct_id AS struct_id,
	sb.acc_residue_number AS acc_residue_number,
	sb.don_residue_number AS don_residue_number,
	fa_atr.score_value AS fa_atr,
	fa_rep.score_value AS fa_rep,
	fa_sol.score_value AS fa_sol,
	fa_pair.score_value AS fa_pair,
	hbond_sc.score_value AS hbond_sc,
	.8*fa_atr.score_value + .44*fa_rep.score_value +
	.65*fa_sol.score_value + .49*fa_pair.score_value +
	1.1*hbond_sc.score_value  AS e_sum
FROM
	CXL_ARG_salt_bridges AS sb,
	residue_scores_2b AS fa_atr,
	residue_scores_2b AS fa_rep,
	residue_scores_2b AS fa_sol,
	residue_scores_2b AS fa_pair,
	residue_scores_2b AS hbond_sc
WHERE
	((fa_atr.struct_id = sb.struct_id AND
	fa_atr.resNum1 = sb.don_residue_number AND
	fa_atr.resNum2 = sb.acc_residue_number AND
	fa_atr.score_type == 'fa_atr') OR 
	(fa_atr.struct_id = sb.struct_id AND
	fa_atr.resNum1 = sb.acc_residue_number AND
	fa_atr.resNum2 = sb.don_residue_number AND
	fa_atr.score_type == 'fa_atr')) AND
	((fa_rep.struct_id = sb.struct_id AND
	fa_rep.resNum1 = sb.don_residue_number AND
	fa_rep.resNum2 = sb.acc_residue_number AND
	fa_rep.score_type == 'fa_rep') OR 
	(fa_rep.struct_id = sb.struct_id AND
	fa_rep.resNum1 = sb.acc_residue_number AND
	fa_rep.resNum2 = sb.don_residue_number AND
	fa_rep.score_type == 'fa_rep')) AND
	((fa_sol.struct_id = sb.struct_id AND
	fa_sol.resNum1 = sb.don_residue_number AND
	fa_sol.resNum2 = sb.acc_residue_number AND
	fa_sol.score_type == 'fa_sol') OR 
	(fa_sol.struct_id = sb.struct_id AND
	fa_sol.resNum1 = sb.acc_residue_number AND
	fa_sol.resNum2 = sb.don_residue_number AND
	fa_sol.score_type == 'fa_sol')) AND
	((fa_pair.struct_id = sb.struct_id AND
	fa_pair.resNum1 = sb.don_residue_number AND
	fa_pair.resNum2 = sb.acc_residue_number AND
	fa_pair.score_type == 'fa_pair') OR 
	(fa_pair.struct_id = sb.struct_id AND
	fa_pair.resNum1 = sb.acc_residue_number AND
	fa_pair.resNum2 = sb.don_residue_number AND
	fa_pair.score_type == 'fa_pair')) AND
	((hbond_sc.struct_id = sb.struct_id AND
	hbond_sc.resNum1 = sb.don_residue_number AND
	hbond_sc.resNum2 = sb.acc_residue_number AND
	hbond_sc.score_type == 'hbond_sc') OR 
	(hbond_sc.struct_id = sb.struct_id AND
	hbond_sc.resNum1 = sb.acc_residue_number AND
	hbond_sc.resNum2 = sb.don_residue_number AND
	hbond_sc.score_type == 'hbond_sc'))
GROUP BY
	sb.struct_id,
	sb.don_residue_number,
	sb.acc_residue_number;
	

CREATE INDEX IF NOT EXISTS hbonds_struct_id_acc_id ON
	hbonds (struct_id, acc_id);

SELECT
	sb.tag,
	sb.struct_id AS struct_id,
	sb.don_residue_number AS don_residue_number,
	sb.acc_residue_number AS acc_residue_number,
	sb.psi,
	CASE WHEN
		sb.psi < ", problem_region$xmax[1], " AND
		", problem_region$xmin, " < sb.psi THEN 1
	ELSE 0 END AS in_bifurcated_region,
	sb_e.fa_atr AS fa_atr,
	sb_e.fa_rep AS fa_rep,
	sb_e.fa_sol AS fa_sol,
	sb_e.fa_pair AS fa_pair,
	sb_e.hbond_sc AS hbond_sc,
	sb_e.e_sum AS sb_interaction_energy,
	don_e.score_value + acc_e.score_value AS res_pair_energy
FROM
	CXL_ARG_salt_bridges AS sb,
	hbonds AS hb1,
	hbonds AS hb2,
	hbond_sites AS don1,
	hbond_sites AS don2,
	CXL_ARG_salt_bridge_energies AS sb_e,
	summed_total_residue_scores AS don_e,
	summed_total_residue_scores AS acc_e
WHERE
	-- hb1 and hb2 are hydrogen bonds that are bifurcated at the acceptor
	-- between an ASP or GLU hbond site to ARG with seqsep > 4
	-- the ASP -> ARG forms a salt bridge which is evaluated
	hb1.struct_id = sb.struct_id AND hb1.acc_id = sb.acc_id AND
	hb2.struct_id = sb.struct_id AND hb2.acc_id = sb.acc_id AND
	don1.struct_id = hb1.struct_id AND don1.site_id = hb1.don_id AND
	don2.struct_id = hb2.struct_id AND don2.site_id = hb2.don_id AND
	don1.site_id != don2.site_id AND
	sb_e.struct_id = sb.struct_id AND
	sb_e.don_residue_number = sb.don_residue_number AND
	sb_e.acc_residue_number = sb.acc_residue_number AND
	don_e.struct_id = sb.struct_id AND
	don_e.resNum = sb.don_residue_number AND
	acc_e.struct_id = sb.struct_id AND
	acc_e.resNum = sb.acc_residue_number;", sep="")
f_bifurcated <- query_sample_sources(sample_sources, sele)
f_bifurcated <- na.omit(f_bifurcated, method="r")
f_bifurcated$salt_bridge_type <- factor("bifurcated")

sele <- paste("
SELECT
	sb1.tag,
	sb1.struct_id AS struct_id,
	sb1.don_residue_number AS don_residue_number,
	sb1.acc_residue_number AS acc_residue_number,
	sb1.psi AS psi,
	sb2.psi AS sb2_psi,
	CASE WHEN
		sb1.psi < ", problem_region$xmax[1], " AND
		", problem_region$xmin, " < sb1.psi THEN 1
	ELSE 0 END AS in_bifurcated_region,
	CASE WHEN
		sb2.psi < ", problem_region$xmax[1], " AND
		", problem_region$xmin, " < sb2.psi THEN 1
	ELSE 0 END AS sb2_in_bifurcated_region,
	sb1_e.fa_atr AS fa_atr,
	sb1_e.fa_rep AS fa_rep,
	sb1_e.fa_sol AS fa_sol,
	sb1_e.fa_pair AS fa_pair,
	sb1_e.hbond_sc AS hbond_sc,
	sb1_e.e_sum,
	sb2_e.fa_atr AS sb2_fa_atr,
	sb2_e.fa_rep AS sb2_fa_rep,
	sb2_e.fa_sol AS sb2_fa_sol,
	sb2_e.fa_pair AS sb2_fa_pair,
	sb2_e.hbond_sc AS sb2_hbond_sc,
	sb1_e.e_sum AS sb_interaction_energy,
	don_e.score_value + acc_e.score_value AS res_pair_energy
FROM
	CXL_ARG_salt_bridges AS sb1, CXL_ARG_salt_bridges AS sb2,
	hbonds AS hb1, hbonds AS hb2,
	hbond_sites AS don1, hbond_sites AS don2,
	CXL_ARG_salt_bridge_energies AS sb1_e,
	CXL_ARG_salt_bridge_energies AS sb2_e,
	summed_total_residue_scores AS don_e,
	summed_total_residue_scores AS acc_e
WHERE
	sb2.struct_id = sb1.struct_id AND
	sb1.acc_id != sb2.acc_id AND
	sb2.don_residue_number = sb1.don_residue_number AND
	sb2.acc_residue_number = sb1.acc_residue_number AND
	hb1.struct_id = sb1.struct_id AND hb1.acc_id = sb1.acc_id AND
	hb2.struct_id = sb2.struct_id AND hb2.acc_id = sb2.acc_id AND
	don1.struct_id = hb1.struct_id AND don1.site_id = hb1.don_id AND
	don2.struct_id = hb2.struct_id AND don2.site_id = hb2.don_id AND
	don1.resNum = sb1.don_residue_number AND
	don2.resNum = sb2.don_residue_number AND
	don1.site_id != don2.site_id AND
	sb1_e.don_residue_number = sb1.don_residue_number AND
	sb1_e.acc_residue_number = sb1.acc_residue_number AND
	sb2_e.struct_id = sb2.struct_id AND
	sb2_e.don_residue_number = sb2.don_residue_number AND
	sb2_e.acc_residue_number = sb2.acc_residue_number AND
	don_e.struct_id = sb1.struct_id AND
	don_e.resNum = sb1.don_residue_number AND
	acc_e.struct_id = sb1.struct_id AND
	acc_e.resNum = sb1.acc_residue_number;", sep="")
f_bidentite <- query_sample_sources(sample_sources, sele)
f_bidentite$salt_bridge_type <- factor("bidentite")
f_bidentite <- na.omit(f_bidentite, method="r")


plot_parts <- list(
	theme_bw(),
	geom_line(aes(x=x, y=y, colour=sample_source)),
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)),
	scale_y_continuous("FeatureDensity"))



plot_id <- "salt_bridge_bifurcated_ASP_GLU_to_ARG_psi"
dens <- estimate_density_1d_wrap(
	f_bifurcated,
	c("sample_source"),
	"psi",
	xlim=c(-180, 180),
	adjust=.2)

p <- ggplot(data=dens) + plot_parts +
	ggtitle("Bifurcated ASP/GLU -> ARG Salt Bridge PSI; SeqSep > 4") +
	geom_rect(data=problem_region, aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf), alpha=.2) +
	scale_x_continuous("Angle around Donor (degrees)")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)



plot_id <- "salt_bridge_bidentite_ASP_GLU_to_ARG_psi"
dens <- estimate_density_1d_wrap(
	f_bidentite,
	c("sample_source"),
	"psi",
	xlim=c(-180, 180),
	adjust=.2)

p <- ggplot(data=dens) + plot_parts +
	ggtitle("ASP/GLU -> ARG Bidentite Salt Bridges; SeqSep > 4") +
	geom_rect(data=problem_region, aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf), alpha=.2) +
	scale_x_continuous("Angle around Donor (degrees)")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)



plot_id <- "salt_bridge_bifurcated_ASP_GLU_to_ARG_psi_by_sb_interaction_energy"
p <- ggplot(data=f_bifurcated) + theme_bw() +
	geom_point(aes(x=psi, y=sb_interaction_energy), size=.8) +
	ggtitle("ASP/GLU -> ARG Bifurcated Salt Bridges; SeqSep > 4") +
	scale_x_continuous("Angle around Donor (degrees)") +
	scale_y_continuous("Weighted Interaction Energy", limits=c(-8, .5)) +
	facet_wrap(~sample_source, ncol=1)
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id <- "salt_bridge_bidentite_ASP_GLU_to_ARG_psi_by_sb_interaction_energy"
p <- ggplot(data=f_bidentite) + theme_bw() +
	geom_point(aes(x=psi, y=sb_interaction_energy), size=.8) +
	ggtitle("ASP/GLU -> ARG Bidentite Salt Bridges; SeqSep > 4") +
	scale_x_continuous("Angle around Donor (degrees)") +
	scale_y_continuous("Weighted Interaction Energy", limits=c(-8, .5)) +
	facet_wrap(~sample_source, ncol=1)
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)



plot_id <- "salt_bridge_bifurcated_ASP_GLU_to_ARG_psi_by_res_pair_energy"
p <- ggplot(data=f_bifurcated) + theme_bw() +
	geom_point(aes(x=psi, y=res_pair_energy), size=.8) +
	ggtitle("ASP/GLU -> ARG Bifurcated Salt Bridges; SeqSep > 4") +
	scale_x_continuous("Angle around Donor (degrees)") +
	scale_y_continuous("Weighted Total Energy for Residue Pair", limits=c(-60, .5)) +
	facet_wrap(~sample_source, ncol=1)
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id <- "salt_bridge_bidentite_ASP_GLU_to_ARG_psi_by_res_pair_energy"
p <- ggplot(data=f_bidentite) + theme_bw() +
	geom_point(aes(x=psi, y=res_pair_energy), size=.8) +
	ggtitle("ASP/GLU -> ARG Bidentite Salt Bridges; SeqSep > 4") +
	scale_x_continuous("Angle around Donor (degrees)") +
	scale_y_continuous("Weighted Total Energy for Residue Pair", limits=c(-8, .5)) +
	facet_wrap(~sample_source, ncol=1)
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
