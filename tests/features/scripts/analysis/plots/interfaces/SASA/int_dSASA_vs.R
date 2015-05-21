# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

check_setup()

feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "int_SASA-dSASA_vs",
author = "Jared Adolf-Bryfogle",
brief_description = "Graphs basic dSASA and SASA information",
feature_reporter_dependencies = c("InterfaceFeatures/AntibodyFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
  
  #First we run on all the interfaces in the database

  

  sele = "
  SELECT
    dSASA,
    dSASA_hphobic,
    dSASA_polar,
    dG,
    interface
  FROM
    interfaces"
  
  int_data = query_sample_sources(sample_sources, sele, char_as_factor=F)
 

  
  #dSASA sides
  sele = "
  SELECT
    dG,
    dSASA,
    dSASA_sc,
    dSASA - dSASA_sc as dSASA_bb,
    dhSASA,
    dhSASA_sc,
    dhSASA - dhSASA_sc as dhSASA_bb,
    dSASA-dhSASA as dpSASA,
    dhSASA_rel_by_charge,
    aromatic_dSASA_fraction,
    interface_nres,
    interface,
    side
  FROM
    interface_sides
  "
  plot_field = function(p, plot_id, grid = NULL){
    
    if (! is.null(grid)){
      p <- p+ facet_grid(facets=grid)
    }
    if(nrow(sample_sources) <= 3){
      p <- p + theme(legend.position="bottom", legend.direction="horizontal")
    }
    save_plots(self, plot_id, sample_sources, output_dir, output_formats)
  }
  
  data = query_sample_sources(sample_sources, sele, char_as_factor=F)
  #print(data)
  
  
  
  #ScatterPlots
  
  
  parts = list(
    geom_point(size=.75),
    #stat_smooth(color="grey"),
    stat_smooth(method=lm),
    geom_density2d(size=.5),
    #stat_density2d(aes(fill = ..level..), geom="polygon"),
    #stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE),
    theme_bw())
  
  fields = c("dSASA", "dSASA_bb", "dSASA_sc", "dhSASA", "dhSASA_bb", "dhSASA_sc", "dhSASA_rel_by_charge")
  for (f in fields){
    
    #dSASA vs dG
    field = f
    p <- ggplot(data=int_data[int_data$dG<=5000 & int_data$dG>-5000,], aes(x = field, y = dG, colour=sample_source)) + parts +
      ggtitle(paste(field,"vs dG")) +
      scale_x_continuous("SASA") +
      scale_y_continuous("REU")
    plot_field(p, paste(field, "vs_dG_by_all", sep="_"))
    plot_field(p, paste(field, "vs_dG_by_interface", sep="_"), grid=interface ~ .)
  
#    #dSASA vs aromatic dSASA fraction.  Shouldn't increase, but worth a plot 
#    p <- ggplot(data= data, aes(x = dSASA, y = aromatic_dSASA_fraction, colour=sample_source)) + parts +
#      ggtitle(paste(field,"vs Aromatic dSASA Fraction")) +
#      scale_x_continuous("SASA") +
#      scale_y_continuous("fraction", limit=c(0, 1.0))
#    plot_field(p, paste(field, "vs_aromatic_dSASA_fraction", sep="_"), grid= side ~ .)
    
#    #dSASA vs interface nres.  Again, shouldn't be interesting.  Mainly a control.
#    p <- ggplot(data = data, aes(y = interface_nres, x = dSASA, colour=sample_source)) + parts +
#      ggtitle(paste(field, "vs Interface nres")) +
#      scale_y_continuous("n") +
#      scale_x_continuous("SASA")
#    plot_field(p, paste("control", field, "vs_interface_nres", sep="_"), grid= side ~ .)
  
    #dSASA vs 'energy density' from Ben Strange's paper - Should be pretty much flat for natives

  }
  
    data$e_density = data$dG/data$dSASA
    field = c("dSASA")
    p <- ggplot(data = data, aes(x = dSASA, y = e_density, colour=sample_source)) + parts +
      ggtitle(paste(field, "vs Interface energy density")) +
      scale_x_continuous("dSASA") +
      scale_y_continuous("dG/dSASA")
    plot_field(p, paste("control", field, "vs_energy_density", sep="_"), grid=side ~ .)
  
#  #dhSASA vs dpSASA
#  p <- ggplot(data = data, aes(y = dpSASA, x = dhSASA, colour=sample_source)) + parts +
#    ggtitle("dhSASA vs dpSASA") +
#    scale_y_continuous("dSASA") +
#    scale_x_continuous("dSASA")
#  plot_field(p, "dhSASA_vs_dpSASA", grid=side ~ .)
  
  #Control -  Should be flat?
  p <- ggplot(data = data, aes(y = dSASA_bb, x = dSASA_sc, colour=sample_source)) + parts +
    ggtitle("dSASA_sc vs dSASA_bb") +
    scale_x_continuous("Sidechain dSASA") +
    scale_y_continuous("Backbone dSASA")
  plot_field(p, "control_dSASA_sc_vs_dSASA_bb", grid=side ~ .)
  
  p <- ggplot(data = data, aes(y = dhSASA_sc, x = dhSASA_bb, colour=sample_source)) + parts +
    ggtitle("dhSASA_sc vs dhSASA_bb") +
    scale_y_continuous("Sidechain dhSASA") +
    scale_x_continuous("Backbone dhSASA")
  plot_field(p, "control_dhSASA_sc_vs_dhSASA_bb", grid=side ~ .)
})) # end FeaturesAnalysis