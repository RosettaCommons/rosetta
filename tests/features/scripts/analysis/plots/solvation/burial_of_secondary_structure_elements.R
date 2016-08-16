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
id = "burial_of_secondary_structure_elements",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("APSAFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <-"
SELECT DISTINCT
  rp.res1_10A_neighbors AS nbrs,
  ss.secondary_struct AS secondary_struct
FROM
  residue_pairs AS rp,
  apsa AS ss
WHERE
  rp.struct_id = ss.struct_id AND
  rp.resNum1 = ss.resNum;"

f <-  query_sample_sources(sample_sources, sele)

dens <- estimate_density_1d(
  data = f,
  ids = c("sample_source", "secondary_struct"),
  variable = "nbrs")

plot_id <- "burial_of_secondary_structure_elements"
p <- ggplot(data=dens, aes(x=log(x+1), y=log(y+1), color=secondary_struct, indicator=counts))
p <- p + geom_line()
p <- p + geom_indicator(aes(group=secondary_struct))
p <- p + facet_wrap( ~ sample_source )
p <- p + ggtitle("Burial by Secondary Structure")
p <- p + labs(x="Number of Neighbors",
              y="log(FeatureDensity + 1)")
p <- p + theme_bw()
p <- p + theme(axis.text.y=theme_blank())

save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
