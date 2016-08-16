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
id = "rotamer_recovery_by_HBChemType",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("RotamerRecoveryFeatures", "HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <-"
SELECT
  rr.divergence AS divergence,
  hb_site_types.label AS hb_chem_type
FROM
  rotamer_recovery AS rr,
  hbond_sites AS hb_site,
  hbond_chem_types AS hb_site_types
WHERE
  hb_site.resNum = rr.resNum AND
  hb_site.struct_id = rr.struct_id AND
  hb_site.HBChemType = hb_site_types.chem_type
ORDER BY RANDOM()
LIMIT 10000;"

f <- query_sample_sources(sample_sources, sele)

exp <- ddply(f, .(sample_source, hb_chem_type),
	function(df) data.frame(mean=round(mean(df$divergence), 2)))


dens <- estimate_density_1d_reflect_boundary(
  data = f,
  ids = c("sample_source", "hb_chem_type"),
  variable = "divergence",
	reflect_left=TRUE, left_boundary=0)

dens <- merge(dens, exp)

plot_id <- "rotamer_recovery_by_HBChemType"
ggplot(data=dens) + theme_bw() +
  geom_line(aes(x=log(x+1), y=log(y+1), colour=sample_source)) +
	geom_indicator(aes(indicator=counts, color=sample_source, group=sample_source)) +
	geom_indicator(aes(indicator=mean, color=sample_source, group=sample_source), xpos="left") +
	facet_wrap( ~ hb_chem_type ) +
	ggtitle("Rotamer Recovery by Hydrogen Bond Chemical Type") +
	labs(x="<- better      log(Automorphic RMSD + 1)      worse ->",
			 y="log(FeatureDensity + 1)")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
