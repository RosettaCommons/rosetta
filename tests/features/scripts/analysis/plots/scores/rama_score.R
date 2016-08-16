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
id = "rama_score",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueScoresFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele_1b <-"
SELECT
	score_type,
	score_value
FROM
	residue_scores_1b
WHERE
	score_type == 'rama' AND
	-2 < score_value AND score_value < 4;"

scores <- query_sample_sources(sample_sources, sele_1b)
scores$score_type <- factor(scores$score_type)

dens <- estimate_density_1d(
  data = scores,
  ids = c("sample_source", "score_type"),
  variable = "score_value")

d_ply(dens, .(score_type), function(sub_dens){
	score_type <- sub_dens$score_type[1]
	plot_id <- paste("component_scores_", score_type, sep="")
	p <- ggplot(data=sub_dens) + theme_bw() +
		geom_line(aes(x=x, y=y, colour=sample_source)) +
		geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
		ggtitle(paste("Rosetta ", score_type, " Scores", sep="")) +
		labs(x="Rosetta Energy Units") +
		scale_y_continuous("FeatureDensity") +
		theme(legend.position=c(.58,.15)) +
		theme(legend.justification=c("right", "top"))
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})


})) # end FeaturesAnalysis
