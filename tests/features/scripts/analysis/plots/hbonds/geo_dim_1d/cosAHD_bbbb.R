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
id = "cosBAH_bbbb",
author = "Matthew O'Meara",

brief_description = "Measure the Acceptor--Hydrogen--Donor angle for backbone-backbone hydrogen bond interactions by sequence separation.",

feature_reporter_dependencies = c("HBondFeatures"),

run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
  geom.cosAHD,
  CASE don_site.resNum - acc_site.resNum
    WHEN -1 THEN '-1' WHEN -2 THEN '-2' WHEN -3 THEN '-3' WHEN -4 THEN '-4'
    WHEN 1 THEN '1' WHEN 2 THEN '2' WHEN 3 THEN '3' WHEN 4 THEN '4'
    ELSE 'long' END AS seq_sep
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
	acc_site.HBChemType == 'hbacc_PBA' AND
	don_site.HBChemType == 'hbdon_PBA';"
f <- query_sample_sources(sample_sources, sele)

dens <- estimate_density_1d_reflect_boundary(
 data=f,
 ids = c("sample_source", "seq_sep"),
 variable = "cosAHD",
 reflect_right=TRUE,
 right_boundary=1)


plot_id <- "cosAHD_bbbb"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=180-180/pi*acos(x), y=y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	facet_wrap( ~ seq_sep ) +
	ggtitle("BB/BB Hydrogen Bonds A-H-D Angle by Sequence Separation\n(donres - accres) normalized for equal weight per unit distance") +
	scale_y_continuous("FeatureDensity" ) +
	scale_x_continuous("Acceptor -- Hydrogen -- Donor Angle (degrees)", limits=c(90,180), breaks=c(90,120,150,180))

if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
