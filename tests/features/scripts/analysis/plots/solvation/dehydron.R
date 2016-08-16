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
id = "dehydron",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


### BY ENERGY ####

sele <-"
SELECT
  hbond.energy,
  dehydron.wrapping_count
FROM
  hbonds AS hbond,
  hbond_dehydrons AS dehydron,
  hbond_sites AS don_site,
  hbond_sites AS acc_site
WHERE
  hbond.struct_id = dehydron.struct_id AND
  hbond.hbond_id = dehydron.hbond_id AND
  hbond.struct_id = don_site.struct_id AND
  hbond.don_id = don_site.site_id AND
  hbond.struct_id = acc_site.struct_id AND
  hbond.acc_id = acc_site.site_id AND
  don_site.HBChemType = 'hbdon_PBA' AND
  acc_site.HBChemType = 'hbacc_PBA' AND
  energy < -.8;"

all_geom <- query_sample_sources(sample_sources, sele)

all_geom$wrapping_count_quantile <-
  cut(all_geom$wrapping_count,
  quantile(all_geom$wrapping_count, seq(0,1,by=.25)))

dens <- estimate_density_1d(
  data = all_geom,
  ids = c("sample_source", "wrapping_count_quantile"),
  variable = "energy")

plot_id <- "hbond_dehydron_energy_qantiles"
ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=wrapping_count_quantile)) +
	geom_indicator(aes(indicator=counts, colour=wrapping_count_quantile, group=wrapping_count_quantile)) +
	facet_wrap(~sample_source, ncol=1) +
	ggtitle("Backbone Backbone Hydrogen Bond Energy by Dehydron Wrapping Count") +
	labs(x="Rosetta Predicted Energy(lower is better geometry)",
	     y="FeatureDensity") +
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


#print(cov(all_geom$energy, all_geom$wrapping_count))

})) # end FeaturesAnalysis
