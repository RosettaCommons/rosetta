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
id = "relative_score_in_environment",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <- "

DROP TABLE IF EXISTS total_residue_scores;
CREATE TEMPORARY TABLE total_residue_scores (
	struct_id INTEGER,
	resNum INTEGER,
	score_value);

CREATE INDEX total_residue_scores_struct_id_resNum ON total_residue_scores (
	struct_id, resNum);

INSERT INTO total_residue_scores
SELECT
	struct_id, resNum,
	SUM(score_value) AS score_value
FROM
	residue_scores_1b
GROUP BY
	struct_id,
	resNum;

INSERT INTO total_residue_scores
SELECT
	struct_id, resNum1 AS resNum,
	SUM(score_value)/2 AS total
FROM
	residue_scores_2b
GROUP BY
	struct_id,
	resNum1;

INSERT INTO total_residue_scores
SELECT
	struct_id, resNum2 AS resNum,
	SUM(score_value)/2 AS total
FROM
	residue_scores_2b
GROUP BY
	struct_id,
	resNum2;

DROP TABLE IF EXISTS residue_vDW_atr;
CREATE TABLE residue_vDW_atr(
	struct_id INTEGER,
	resNum INTEGER,
	vDW_atr);

CREATE INDEX residue_vDW_atr_struct_id_resNum ON residue_vDW_atr (
	struct_id, resNum);

INSERT INTO residue_vDW_atr
SELECT
	struct_id, resNum1 AS resNum,
	SUM(score_value) AS vDW_atr
FROM
	residue_scores_2b
WHERE
	score_type = 'fa_atr'
GROUP BY
	struct_id,
	resNum1;

INSERT INTO residue_vDW_atr
SELECT
	struct_id, resNum2 AS resNum,
	SUM(score_value) AS vDW_atr
FROM
	residue_scores_2b
WHERE
	score_type = 'fa_atr'
GROUP BY
	struct_id,
	resNum2;


SELECT
	res.name3 AS res_type,
	SUM(vDW_atr.vDW_atr) AS vDW_atr,
	SUM(total.score_value) AS total
FROM
	residues AS res,
	total_residue_scores AS total,
	residue_vDW_atr AS vDW_atr
WHERE
	total.struct_id = res.struct_id AND
	total.resNum = res.resNum AND
	vDW_atr.struct_id = res.struct_id AND
	vDW_atr.resNum = res.resNum
GROUP BY
	res.struct_id,
	res.resNum;"


f <- query_sample_sources(sample_sources, sele)

# remove residues that are clashing
clean_f <- f[f$total < 5,]
clean_f <- clean_f[clean_f$res_type != "CA",]

sub_f <- sample_rows(clean_f, 8000)

plot_id <- "average_residue_score_in_environment"
p <- ggplot() + theme_bw() +
	geom_point(data=sub_f, aes(x=vDW_atr, y=total), size=.5) +
	stat_smooth(data=sub_f, aes(x=vDW_atr, y=total), method="lm") +
	facet_wrap(~res_type, scales="free_x") +
	ggtitle("Total residue score as a function of Van der Waals attraction ") +
	labs(x="Residue Total Van der Waals Attration (Rosetta Energy Units)") +
	scale_y_continuous("Residue Total Energy (Rosetta Energy Units)")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

mean_score_model = lm(f, total ~ vDW_atr)

f$diff_avg_score <- predict(mean_score_model)


#write.csv(f, file="diff_avg_score.csv", sep=", ", row.names=F)



})) # end FeaturesAnalysis
