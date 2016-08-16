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
id = "hbond_scores",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueFeatures", "ResidueScoreFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <- "
SELECT
	res1.name3 AS res1_type, res2.name3 AS res2_type,
	s2b.score_type AS sc_type, s2b.score_value AS score
FROM
	residue_scores_2b AS s2b,
	residues AS res1, residues AS res2
WHERE

	res1.struct_id = s2b.struct_id AND res1.resNum = s2b.resNum1 AND
	res2.struct_id = s2b.struct_id AND res2.resNum = s2b.resNum2 AND
	(s2b.score_type = 'hbond_bb_sc' OR s2b.score_type = 'hbond_lr_bb_sc' OR
	s2b.score_type = 'hbond_sc' OR s2b.score_type = 'hbond_sr_bb_sc') AND
	s2b.score_value > -4;"


f <- query_sample_sources(sample_sources, sele)

dens <- estimate_density_1d(f, c("sample_source", "sc_type"), "score")

plot_id <- "hbond_scores"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x,y, color=sample_source)) +
	geom_indicator(aes(indicator=counts, color=sample_source, group=sample_source)) +
	ggtitle("Unweighted Rosetta HBond Scores") +
	facet_wrap( ~ sc_type ) +
	labs(x="Unweighted Rosetta Energy", y = "Feature Density")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
