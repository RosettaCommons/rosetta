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
id = "rotamer_recovery_by_secondary_structure",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("RotamerRecoveryFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <-"
SELECT
  rr.divergence AS divergence,
  ss.secondary_struct AS secondary_struct
FROM
  rotamer_recovery AS rr,
  apsa AS ss
WHERE
  ss.struct_id  = rr.struct_id AND
  ss.resNum = rr.resNum;"

f <- query_sample_sources(sample_sources, sele)

dens <- estimate_density_1d(
  data = f,
  ids = c("secondary_struct"),
  variable = "divergence")

plot_id <- "rotamer_recovery_by_secondary_structure"
p <- ggplot(data=dens, aes(x=log(x+1), y=log(y+1), color=secondary_struct, indicator=counts))
p <- p + geom_line()
p <- p + geom_indicator(aes(group=secondary_struct))
p <- p + ggtitle("Rotamer Recovery by Secondary Structure")
p <- p + labs(x="<- better      log(Automorphic RMSD + 1)      worse ->",
              y="log(FeatureDensity + 1)")
p <- p + theme_bw()

save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
