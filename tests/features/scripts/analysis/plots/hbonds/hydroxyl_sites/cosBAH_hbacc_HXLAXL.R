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
id = "cosBAH_hbacc_HXLAXL",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
  geom.cosBAH,
  acc_site.HBChemType AS acc_chem_type,
  don_site.HBChemType AS don_chem_type
FROM
  hbond_geom_coords AS geom,
  hbonds AS hbond,
  hbond_sites AS don_site,
  hbond_sites AS acc_site
WHERE
  geom.cosBAH < 0 AND
  (acc_site.HBChemType = 'hbacc_HXL' OR acc_site.HBChemType = 'hbacc_AHX') AND
	don_site.HBChemType != 'hbdon_PBA' AND
  hbond.struct_id = geom.struct_id AND
  hbond.hbond_id =  geom.hbond_id AND
  hbond.struct_id = don_site.struct_id AND
  hbond.don_id = don_site.site_id AND
  hbond.struct_id = acc_site.struct_id AND
  hbond.acc_id = acc_site.site_id AND
  ABS(don_site.resNum - acc_site.resNum) > 5;"

f <- query_sample_sources(sample_sources, sele)
print(summary(f))

dens <- estimate_density_1d(
  data = f, ids = c("sample_source" ), variable = "cosBAH")

plot_id = "cosBAH_scsc_sp3acc"
ggplot(data=dens) +
	geom_line(aes(x=(pi-acos(x))*180/pi, y=y, colour=sample_source),size=3) +
	geom_indicator(aes(colour=sample_source, indicator=counts, group=sample_source)) +
	ggtitle("Sidechain/Sidechain Hydrogen Bonds to SP3 Acceptors. BAH Angle by Chemical Type\n(normalized for equal volume per unit distance)") +
	scale_x_continuous(paste('Base -- Acceptor -- Hydrogen (degrees)')) +
	scale_y_continuous("FeatureDensity")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
