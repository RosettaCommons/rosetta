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
id = "ab_packing_angle_den",
author = "Jared Adolf-Bryfogle",
brief_description = "VL VH packing angle metrics",
feature_reporter_dependencies = c("AntibodyFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
  
  #First we run on all the interfaces in the database

  

  sele = "
  SELECT
    VL_VH_packing_angle,
    VL_VH_distance,
    VL_VH_opening_angle,
    VL_VH_opposite_opening_angle
  FROM
    ab_metrics
    "
  
  data = query_sample_sources(sample_sources, sele, char_as_factor=F)
 
  parts <- list(
    geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)),
    scale_y_continuous("Feature Density"),
    theme_bw())

  plot_field = function(p, plot_id, grid = NULL){
    
    if (! is.null(grid)){
      p <- p+ facet_wrap(facets=grid, ncol=3)
    }
    if(nrow(sample_sources) <= 3){
      p <- p + theme(legend.position="bottom", legend.direction="horizontal")
    }
    save_plots(self, plot_id, sample_sources, output_dir, output_formats)
  }
  
  
  #Packing Angle
  group = c("sample_source")
  dens <- estimate_density_1d(data, group, c("VL_VH_packing_angle"))
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    xlab("Degrees") +
    ggtitle("VL VH Packing Angle")
  plot_field(p, "vl_vh_packing_angle_den")
  
  #Distance
  group = c("sample_source")
  dens <- estimate_density_1d(data, group, c("VL_VH_distance"))
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    xlab("Angstroms") +
    ggtitle("VL VH Distance")
  plot_field(p, "vl_vh_distance_den")
  
  #Opening Angle
  group = c("sample_source")
  dens <- estimate_density_1d(data, group, c("VL_VH_opening_angle"))
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    xlab("Degrees") +
    ggtitle("VL VH Opening Angle")
  plot_field(p, "vl_vh_opening_angle_den")
  
  #Opposite Opening angle
  group = c("sample_source")
  dens <- estimate_density_1d(data, group, c("VL_VH_opposite_opening_angle"))
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    xlab("Degrees") +
    ggtitle("VL VH Opposite Opening Angle")
  plot_field(p, "vl_vh_opposite_opening_angle_den")
  
})) # end FeaturesAnalysis