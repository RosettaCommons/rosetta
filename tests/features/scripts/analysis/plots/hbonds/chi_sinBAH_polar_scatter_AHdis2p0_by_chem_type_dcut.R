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
id = "chi_sinBAH_polar_scatter_AHdis2p0_by_chem_type_dcut",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

d_ply(sample_sources, .variables=("sample_source"), function(sample_source){
  ss <- sample_source[1,"sample_source"]

  sele <-"
  SELECT
    geom.AHdist,
    geom.cosBAH,
    geom.chi,
    acc_site.HBChemType AS acc_chem_type,
    don_site.HBChemType AS don_chem_type
  FROM
    hbond_geom_coords AS geom,
    hbonds AS hbond,
    hbond_sites AS don_site,
    hbond_sites AS acc_site
  WHERE
    geom.AHdist < 2.0 AND
    hbond.struct_id = geom.struct_id AND hbond.hbond_id = geom.hbond_id AND
    hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
    hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
		acc_site.HBChemType = 'hbacc_CXL' AND don_site.HBChemType = 'hbdon_HXL' AND
		ABS( don_site.resNum - acc_site.resNum ) > 5;";

  f <- query_sample_sources(sample_source, sele)

  #equal area projection
  f <- transform(f,
    capx = 2*sin(acos(cosBAH)/2)*cos(chi),
    capy = 2*sin(acos(cosBAH)/2)*sin(chi))

  sub_f <- ddply(f, .variables=c("don_chem_type", "acc_chem_type"),
    function(df){sample_rows(df, 5000)})

  plot_id = "hbond_sinBAH_longrange_dAH2p0cut_colorbydist_eq_polar_scatter_by_chem_type"
  ggplot(data=sub_f) + theme_bw() +
    geom_point(aes(x=capx, y=capy, color=AHdist), size=1.0) +
    facet_grid(acc_chem_type ~ don_chem_type) +
    ggtitle(paste("Hydrogen Bonds chi vs sinBAH Angles by Chemical Type\nEqual Coordinate Projection   Sample Source: ", ss, sep="")) +
    scale_x_continuous('2*sin(BAH/2) * cos(CHI)', breaks=c(-1, 0, 1)) +
    scale_y_continuous('2*sin(BAH/2) * sin(CHI)', breaks=c(-1, 0, 1))
  save_plots(self, plot_id, sample_source, output_dir, output_formats)

  #orthographic projection
#  f <- transform(f,
#    capx = sin(acos(cosBAH))*cos(chi),
#    capy = sin(acos(cosBAH))*sin(chi))

#  sub_f <- ddply(f, .variables=c("don_chem_type", "acc_chem_type"),
#    function(df){sample_rows(df, 5000)})

#  plot_id = "hbond_sinBAH_ortho_polar_scatter_by_chem_type"
#  ggplot(data=sub_f) + theme_bw() +
#    geom_point(aes(x=capx, y=capy), size=.5) +
#    facet_grid(acc_chem_type ~ don_chem_type) +
#    ggtitle(paste("Hydrogen Bonds chi vs sinBAH Angles by Chemical Type\nOrthographic Projection   Sample Source: ", ss, sep="")) +
#    scale_x_continuous('sin(BAH) * cos(CHI)', breaks=c(-1, 0, 1)) +
#    scale_y_continuous('sin(BAH) * sin(CHI)', breaks=c(-1, 0, 1))
#  save_plots(self, plot_id, sample_source, output_dir, output_formats)

})


})) # end FeaturesAnalysis
