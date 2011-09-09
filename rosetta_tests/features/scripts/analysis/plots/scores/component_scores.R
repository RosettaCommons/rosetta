# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


check_setup()

sele_1b <-"
SELECT score_type, score_value
FROM   residue_scores_1b
WHERE  -2 < score_value AND score_value < 5
ORDER BY RANDOM()
LIMIT 10000;"

sele_2b <-"
SELECT score_type, score_value
FROM   residue_scores_2b
WHERE  -2 < score_value AND score_value < 5
ORDER BY RANDOM()
LIMIT  100000;"

scores <- rbind(
  query_sample_sources(sample_sources, sele_1b),
  query_sample_sources(sample_sources, sele_2b))
scores$score_type <- factor(scores$score_type)

dens <- estimate_density_1d(
  data = scores,
  ids = c("sample_source", "score_type"),
  variable = "score_value")

plot_id <- "component_scores"
p <- ggplot(data=dens,
  aes(x=x, y=log(y+1), colour=sample_source, indicator=counts))
p <- p + geom_line() + geom_indicator()
p <- p + facet_wrap( ~ score_type, ncol = 4)
p <- p + opts(title = "Rosetta Component Scores")
p <- p + labs(x="Rosetta Energy Units", y="log(FeatureDensity + 1)")
p <- p + theme_bw()
p <- p + scale_y_continuous(breaks=c(0, .3, .6))

save_plots(plot_id, sample_sources, output_dir, output_formats)
