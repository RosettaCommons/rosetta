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
id = "residue_pair_distance_vs_num_neighbors",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueFeatures", "ResiduePairFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
  rp.actcoord_dist AS dist,
  rp.res1_10A_neighbors + rp.res2_10A_neighbors AS nbrs,
  res1.name3 AS res1_type,
  res2.name3 AS res2_type
FROM
  residue_pairs AS rp,
  residues AS res1,
  residues AS res2
WHERE
  res1.resNum = rp.resNum1 AND
  res2.resNum = rp.resNum2 AND
  res1.struct_id = rp.struct_id AND
  res2.struct_id = rp.struct_id
ORDER BY
  RANDOM()
LIMIT
  100000;"

f <-  query_sample_sources(sample_sources, sele)

plot_id <- "residue_pair_distances_vs_neighbors"
p <- ggplot(data=f, aes(x=dist, y=nbrs, color=sample_source))
p <- p + geom_point(aes( size = .5)
p <- p + facet_wrap( ~ res1_type )
p <- p + ggtitle("Residue Pair Distances vs Burial")
p <- p + labs(x="Distance between Action Coordinates",
              y="Sum of Number of Neighbors")
p <- p + theme_bw()
p <- p + theme(axis.text.y=theme_blank())

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


										
})) # end FeaturesAnalysis
