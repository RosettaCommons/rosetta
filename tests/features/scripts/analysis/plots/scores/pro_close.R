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
id = "pro_close",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueFeatuers"),
run=function(self, sample_sources, output_dir, output_formats){



sele <-"
SELECT
  score_value
FROM
  residue_scores_2b
WHERE
  score_type = 'pro_close' AND
  score_value < 3.2;"

f <- query_sample_sources(sample_sources, sele)

dens <- estimate_density_1d(
  data = f,
  ids = c("sample_source"),
  variable = "score_value")

plot_id <- "pro_close"
p <- ggplot(data=dens, aes(x=x, y=log(y+1), colour=sample_source, indicator=counts))
p <- p + geom_line()
p <- p + geom_indicator(aes(group=sample_source))
p <- p + ggtitle("Proline Closure")
p <- p + theme_bw()
p <- p + scale_y_continuous("Score", limits=c(0,1.1))
p <- p + scale_x_continuous("log(FeatureDensity + 1)", limits=c(0,3.2))

save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
