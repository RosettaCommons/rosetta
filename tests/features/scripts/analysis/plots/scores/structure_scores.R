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
id = "structure_scores",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("StructureScoreFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

# TODO join score_type_id with the score_types table
sele <-"
SELECT
	score_value
FROM
	structure_scores
WHERE
	(score_type_id = 262 OR score_type_id = 245) AND
	score_value < 20;"

f <- query_sample_sources(sample_sources, sele)

print(summary(f))

dens <- estimate_density_1d(
  data = f,
  ids = c("sample_source"),
  variable = "score_value")

plot_id <- "structure_scores"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source), size=1.4) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	ggtitle("Rosetta Structure Scores Scores") +
	labs(x="Rosetta Energy Units") +
	scale_y_continuous("FeatureDensity", breaks=c(0, .3, .6))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
