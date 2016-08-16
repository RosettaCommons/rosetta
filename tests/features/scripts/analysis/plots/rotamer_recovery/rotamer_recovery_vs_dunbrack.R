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
id = "rotamer_recovery_vs_dunbrack",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueFeatures", "ResidueScoresFeatures", "RotamerRecoveryFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <- "
SELECT
  rr.divergence AS divergence,
  sc.score_value,
  res.name3 AS res_type
FROM
  rotamer_recovery AS rr,
  residue_scores_1b AS sc,
  residues AS res
WHERE
  rr.resNum = res.resNum AND
  rr.struct_id = res.struct_id AND
  sc.resNum = res.resNum AND
  sc.struct_id = res.struct_id AND
  sc.score_type = 'fa_dun' AND
  sc.score_value < 30;"

f <-  query_sample_sources(sample_sources, sele)

plot_id <- "rotamer_recovery_vs_dunbrack"
p <- ggplot(
  data=f,
  aes(x=score_value, y=divergence))
p <- p + geom_point(size=.3)
p <- p + facet_wrap( ~ res_type )
p <- p + ggtitle(paste("Rotamer Recovery vs Dunbrack Score by Residue Type"))
p <- p + labs(x="Score",
              y="AutomorphicRMSD")
p <- p + theme_bw()

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
