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
id = "buried_unsatisfied",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
  site.HBChemType as chem_type,
  COUNT( CASE WHEN site_env.sasa_r100 = 0 THEN 1 ELSE NULL END ) AS small_probe,
  COUNT( CASE WHEN site_env.sasa_r140 = 0 THEN 1 ELSE NULL END ) AS medium_probe,
  COUNT( CASE WHEN site_env.sasa_r200 = 0 THEN 1 ELSE NULL END ) AS large_probe
FROM
  hbond_sites AS site,
  hbond_site_environment AS site_env
WHERE
  site.struct_id = site_env.struct_id AND site.site_id = site_env.site_id AND
  site_env.num_hbonds = 0
GROUP BY
  site.HBChemType;"
f <- query_sample_sources(sample_sources, sele)

print(summary(f))

f <- melt(f, id.vars=c("sample_source", "chem_type"), variable_name="probe_radius")
names(f) <- c("sample_source", "chem_type", "probe_radius", "counts")

plot_id <- "hb_unsat"
ggplot(f) + theme_bw() +
  geom_bar(aes(x=sample_source, y=log(counts), fill=sample_source)) +
  facet_grid(chem_type ~ probe_radius) +
  ggtitle("Number of Buried Unsatisfied Hydrogen Bond Sites for Different Probe Radii") +
  coord_flip() +
	theme(strip.text.y=theme_text(size=7)) +
	theme(legend.position="none")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
