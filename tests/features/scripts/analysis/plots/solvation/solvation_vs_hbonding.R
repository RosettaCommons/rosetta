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
id = "solvation_vs_hbonding",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures", "GeometricSolvationFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	geo_sol.geometric_solvation * hb_site_env.hbond_energy AS gs_hb,
	hb_site.HBChemType AS hb_chem_type
FROM
	hbond_site_environment AS hb_site_env,
	hbond_sites AS hb_site,
	geometric_solvation AS geo_sol
WHERE
  hb_site_env.struct_id = hb_site.struct_id AND
	hb_site_env.site_id = hb_site.site_id AND
	geo_sol.struct_id  = hb_site.struct_id AND
	geo_sol.hbond_site_id = hb_site.site_id AND
  gs_hb > -6;"

f <- query_sample_sources(sample_sources, sele)


plot_id <- "solvation_vs_hbonding"

dens <- estimate_density_1d(
  f, c("sample_source", "hb_chem_type"), "gs_hb")

ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=-log(1-x), y=log(y+1), colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	facet_wrap( ~ hb_chem_type) +
	ggtitle("Geometric Solvation score * HBonding score") +
	labs(x="(Geometric Solvation)*(Hydrogen Bonding)",
	     y="log(FeatureDensity + 1)")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
