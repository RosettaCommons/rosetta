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
id = "cosBAH_regression",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
  geom.cosBAH,
  acc_site.HBChemType AS acc_chem_type,
  don_site.HBChemType AS don_chem_type,
  CASE acc_site.HBChemType
		WHEN 'hbacc_IMD' THEN 'ring' WHEN 'hbacc_IME' THEN 'ring'
		WHEN 'hbacc_AHX' THEN 'sp3'  WHEN 'hbacc_HXL' THEN 'sp3'
		WHEN 'hbacc_CXA' THEN 'sp2'  WHEN 'hbacc_CXL' THEN 'sp2'
		WHEN 'hbacc_PBA' THEN 'bb_sp2'  END AS hybrid,
  CASE ABS(don_site.resNum - acc_site.resNum) > 5
		WHEN 0 THEN 'short_range' WHEN 1 THEN 'long_range' END AS seq_sep_geq_6
FROM
  hbond_geom_coords AS geom,
  hbonds AS hbond,
  hbond_sites AS don_site,
  hbond_sites AS acc_site
WHERE
  hbond.struct_id = geom.struct_id AND hbond.hbond_id = geom.hbond_id AND
  hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
  hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id;"

f <- query_sample_sources(sample_sources, sele)

#mirror data across "linear"
f <- rbind(f, transform(f, cosBAH = 2 - cosBAH))

plot_parts <- list(
	theme_bw(),
	geom_vline(x=1, colour="darkgrey", size=.9),
  scale_x_continuous("Acceptor Base -- Acceptor -- Hydrogen cos(degrees)"),
  scale_y_continuous("-log(FeatureDensity)"))

######## Join all types of hydrogen bonds together in a single distribution
dens <- estimate_density_1d(f, c("sample_source"), "cosBAH")
dens$neg_log_y = -log(dens$y)

plot_id <- "hbond_cosBAH_regression_unified"
ggplot(dens) + plot_parts +
  geom_line(aes(x=x, y=neg_log_y, color=sample_source)) +
  ggtitle("Hydrogen Bonds cosBAH\nnormalized for equal weight per unit distance") +
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


#####  Fit by hybrid and sequence separation####

dens <- estimate_density_1d(
	f, c("sample_source", "hybrid", "seq_sep_geq_6"), "cosBAH")
dens$neg_log_y = -log(dens$y)

plot_id <- "HBond_cosBAH_regression_by_range_hybrid"
ggplot(dens) + plot_parts +
  geom_line(aes(x=x, y=neg_log_y, colour=sample_source)) +
  geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
  facet_grid(seq_sep_geq_6 ~ hybrid) +
  ggtitle("Hydrogen Bonds cosBAH by Hybridization and Sequence Separation\nnormalized for equal weight per unit distance") +
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

#####  Fit by hybrid and sequence separation####

d_ply(f, c("sample_source"), function(s_f){
	ss <- s_f[1,"sample_source"]

	dens <- estimate_density_1d(
		s_f, c("hybrid", "seq_sep_geq_6", "don_chem_type"), "cosBAH")
	dens$neg_log_y = -log(dens$y)

	plot_id <- "HBond_cosBAH_regression_by_range_hybrid_and_don_chem_type"
	ggplot(dens) + plot_parts +
	  geom_line(aes(x=x, y=neg_log_y, colour=don_chem_type)) +
	  geom_indicator(aes(indicator=counts, colour=don_chem_type, group=don_chem_type)) +
	  facet_grid(seq_sep_geq_6 ~ hybrid) +
	  ggtitle(paste("Hydrogen Bonds cosBAH by Hybridization and Sequence Separation\nnormalized for equal weight per unit distance  Sample Source: ", ss, sep="")) +
		scale_y_continuous("-log(FeatureDensity)", limits=c(-2, 3))
	save_plots(self, plot_id, sample_sources[sample_sources$sample_source==ss,], output_dir, output_formats)
})



})) # end FeaturesAnalysis
