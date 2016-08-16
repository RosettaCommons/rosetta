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
id = "radius_of_gyration",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies  = c("PoseConformationFeatures", "RadiusOfGyrationFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

# if you didn't extract PoseConformationFeatures you could compute
# total_residue from the residues table, but it'd be slower.
sele <-"
SELECT
  pc.total_residue,
  rg.radius_of_gyration
FROM
  pose_conformations AS pc,
  radius_of_gyration AS rg
WHERE
  pc.struct_id = rg.struct_id;"

f <- query_sample_sources(sample_sources, sele)

plot_id <- "radius_of_gyration_scatter"
ggplot(data=f) + theme_bw() +
  geom_point(aes(x=total_residue, y=radius_of_gyration, colour=sample_source)) +
  ggtitle("Radius of Gyration by Number of Residues") +
  labs(y="Radius of Gyration", x="Number of Residues")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

f <- transform(f, normed_radius_of_gyration = radius_of_gyration / total_residue)

#filter out non-compacted structures
f <- f[f$normed_radius_of_gyration < .3,]
dens <- estimate_density_1d(f, c("sample_source"), "normed_radius_of_gyration")

plot_id <- "radius_of_gyration"
ggplot(data=dens) + theme_bw() +
  geom_line(aes(x=x, y=y, color=sample_source)) +
  geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
  ggtitle("Normalized Radius of Gyration") +
  labs(x="Radius of Gyration / Number of Residues", y="Feature Density")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
