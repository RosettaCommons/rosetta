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
id = "number_of_neighbors",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueFeatures", "ResidueBurialFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


# for now find # neighbors through the residue pair table
sele <-"
SELECT
  res.name3 as res_type,
  bur.ten_a_neighbors as nbrs
FROM
  residues as res,
  residue_burial as bur
WHERE
  bur.struct_id == res.struct_id AND bur.resNum == res.resNum;"

f <-  query_sample_sources(sample_sources, sele)

dens <- estimate_density_1d(
  data = f,
  ids = c("sample_source", "res_type"),
  variable = "nbrs")

plot_id <- "number_of_neighbors"
p <- ggplot(data=dens, aes(x=x, y=log(y+1), color=sample_source, indicator=counts))
p <- p + geom_line()
p <- p + geom_indicator(aes(group=sample_source))
p <- p + facet_wrap( ~ res_type )
p <- p + ggtitle("Residue Burial by Residue Type")
p <- p + labs(x="Number of Neighbors in 10 Angstroms",
              y="log(FeatureDensity + 1)")
p <- p + theme_bw()
p <- p + theme(axis.text.y=theme_blank())

save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
