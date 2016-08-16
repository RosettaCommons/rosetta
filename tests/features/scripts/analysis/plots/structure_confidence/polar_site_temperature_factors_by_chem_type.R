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
id = "polar_site_temperature_factors_by_chem_type",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <- "
SELECT
  pdb.heavy_atom_temperature AS temperature,
	label.label AS chem_type
FROM
	hbond_sites AS site,
	hbond_chem_types AS label,
	hbond_sites_pdb AS pdb
WHERE
	site.struct_id = pdb.struct_id AND site.site_id = pdb.site_id AND
	site.HBChemType = label.chem_type;"

# Execute the SQL query on each sample source.
f <- query_sample_sources(sample_sources, sele)


plot_id <- "polar_site_temperature_factor_by_chem_type_ecdf"
f <- f[order(f$temperature),]
n_y_bins = 200
temperature_quantiles <- ddply(f, .(sample_source, chem_type), function(df){
	data.frame(
		x=quantile(df$temperature, seq(0,1,length.out=n_y_bins)),
		y=seq(0,1,length.out=n_y_bins),
		counts=nrow(df))
})

ggplot(data=temperature_quantiles) + theme_bw() +
  geom_line(aes(x=x, y=y, color=sample_source)) +
  facet_wrap( ~ chem_type ) +
  geom_indicator(aes(indicator=counts, color=sample_source, group=sample_source), xpos="left") +
  ggtitle("Cumulative distribution of H-Bond heavy atom temperature factor by chemical type") +
	scale_x_continuous("Temperature Factor", limit=c(0, 60)) +
	scale_y_continuous("Fraction of polar sites less than temperature factor cut off")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)





# PLOT DENSITY DISTRIBUTION OF BFACTOR VALUES FOR POLAR SITES

#plot_id <- "polar_site_temperature_factor_by_chem_type_density"
#dens <- estimate_density_1d(f, c("sample_source", "chem_type"), "temperature")
## Generate a lattice of density plots for each chemical type type
#p <- ggplot(dens) + theme_bw() +
#	geom_line(aes(x, log(y+1), colour=sample_source)) +
#	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
#	facet_wrap(~chem_type) +
#	ggtitle("Polar Site Heavy Atom Temperature Factor By Chemical Type") +
#	scale_x_continuous("Temperature Factor")
#	scale_y_continuous("log(FeatureDensity + 1)")
#if(nrow(sample_sources) <= 3){
#	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
#}
#
#save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
