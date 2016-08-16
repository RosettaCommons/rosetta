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
id = "zero_score_cusp",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
  score_type,
  score_value
FROM
  residue_scores_1b
WHERE
  abs(score_value) < .1;"

residue_scores_1b <- query_sample_sources(sample_sources, sele)

sele <-"
SELECT
  score_type,
  score_value
FROM
  residue_scores_2b
WHERE
  abs(score_value) < .1;"

residue_scores_2b <- query_sample_sources(sample_sources, sele)

scores <- rbind( residue_scores_1b, residue_scores_2b )

scores$score_type <- factor(scores$score_type)

dens <- estimate_density_1d(
  data = scores,
  ids = c("sample_source", "score_type"),
  variable = "score_value")

plot_id <- "zero_score_cusp"
p <- ggplot(data=dens, aes(x=x, y=log(y+1), colour=sample_source, indicator=counts))
p <- p + geom_line()
p <- p + geom_indicator(aes(group=sample_source))
p <- p + facet_wrap( ~ score_type, ncol = 3)
p <- p + ggtitle("Residue Level Scores")
p <- p + labs(x=expression(paste('Acceptor -- Donor Distance (', ring(A), ')')),
              y="log(FeatureDensity + 1)")
p <- p + theme_bw()

save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
