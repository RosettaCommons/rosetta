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
id = "hbond_AHdist_cosAHD_regression",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
  geom.AHdist, geom.cosAHD,
  acc_site.HBChemType AS acc_chem_type,
  don_site.HBChemType AS don_chem_type
FROM
  hbond_geom_coords AS geom,
  hbonds AS hbond,
  hbond_sites AS don_site,
  hbond_sites AS acc_site
WHERE
  hbond.struct_id = geom.struct_id AND
  hbond.hbond_id =  geom.hbond_id AND
  hbond.struct_id = don_site.struct_id AND
  hbond.don_id = don_site.site_id AND
  hbond.struct_id = acc_site.struct_id AND
  hbond.acc_id = acc_site.site_id AND
  hbond.donRank <= 1 AND hbond.accRank <= 1;"

f <- query_sample_sources(sample_sources, sele)

d <- estimate_density_2d(f, c("sample_source", "don_chem_type", "acc_chem_type"), "cosAHD", "AHdist", histogram=TRUE, n_pts=100)

#d_ply(f, c("sample_source", "don_chem_type", "acc_chem_type"), function(df){
#	d[d$sample_source==df$sample_source[1],
#	  d$don_chem_type==df$don_chem_type[1],
#    d$acc_chem_type==df$acc_chem_type[1],]$cov <- cov(f$cosAHd, f$AHdist)
#})

plot_id <- "hbond_AHdist_cosAHD_density"
ggplot(data=d, aes(x=x, y=y, fill=-log(z+1))) + theme_bw() +
  geom_tile() +
#  geom_indicator(aes(indicator=cov)) +
  facet_grid(don_chem_type ~ acc_chem_type) +
  scale_fill_gradient("AHdist", high="lightyellow", low="blue") +
  ggtitle("Hydrogen Bonds AHdist vs cosAHD Angle Fit with Beta Function\n(normalized for equal volume per unit distance)") +
  labs(x=expression(paste('cos(Acceptor -- Hydrogen -- Donor)')),
       y=expression(paste('(Acceptor -- Proton Distance)(', ring(A), ')')))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)



})) # end FeaturesAnalysis
