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
id = "interface_dSASA",
author = "Jared Adolf-Bryfogle",
brief_description = "Graphs basic dSASA and SASA information",
feature_reporter_dependencies = c("InterfaceFeatures"),
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
  
  data = query_sample_sources(sample_sources, sele, char_as_factor=F)
 
  plot_parts <- list(
    geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)),
    scale_y_continuous("Feature Density"),
    theme_bw())
  
  plot_field = function(p, plot_id, grid = NULL, save=T){

    if (! is.null(grid)){
      p <- p+ facet_grid(facets=grid)
    }
    if(nrow(sample_sources) <= 3){
      p <- p + theme(legend.position="bottom", legend.direction="horizontal")
    }
    if(save){
      save_plots(self, plot_id, sample_sources, output_dir, output_formats)
    }
  }

  
  parts = list(plot_parts, scale_x_continuous("SASA"))
  fields = c("dSASA", "dSASA_hphobic", "dSASA_polar")
  for(field in fields){
    
    group = c("sample_source")
    dens <- estimate_density_1d(data,  group, field)
    p <- ggplot(data=dens, na.rm=T) + parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(paste("Buried", field, sep=" "))
    plot_field(p, paste(field, "den", sep="_"))
    
    group = c("sample_source", "interface")
    dens <- estimate_density_1d(data,  group, field)
    p <- ggplot(data=dens, na.rm=T) + parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(paste("Buried", field, sep=" "))
    plot_field(p, paste(field, "den_by_interface", sep="_"),grid=~interface)
  }


  
  #dSASA sides
  int_data = data
  
  sele = "
  SELECT
    dSASA,
    avg_per_residue_dSASA,
    avg_per_residue_SASA_int,
    avg_per_residue_SASA_sep,
    aromatic_dSASA_fraction,
    interface_nres,
    interface,
    side
  FROM
    interface_sides
  "
  print(sele)
  data = query_sample_sources(sample_sources, sele, char_as_factor=F)
  #print(data)
  
  
  field = "dSASA"
  parts = list(plot_parts, scale_x_continuous("SASA"))
  group = c("sample_source", "side")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle(paste("Buried", field, sep=" "))
  plot_field(p, paste(field, "den_sides", sep="_"), grid=side ~ .)
  
  group = c("sample_source", "interface", "side")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle(paste("Buried", field, sep=" "))
  plot_field(p, paste(field, "den_sides","by_interface", sep="_"), grid=side~interface)
  
  
  #AvgFields
  fields = c("avg_per_residue_dSASA",
            "avg_per_residue_SASA_int",
            "avg_per_residue_SASA_sep")

  parts = list(plot_parts, scale_x_continuous("SASA", limit=c(0, 100)))
  for(field in fields){
    fieldSP = unlist(strsplit(field, split="_"))
    group = c("sample_source", "side")
    dens <- estimate_density_1d(data,  group, field)
    p <- ggplot(data=dens, na.rm=T) + parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(paste(capwords(fieldSP), collapse=" "))
    plot_field(p, paste(field, "den_sides", sep="_"), grid=side ~ .)
    
    group=c("sample_source", "interface", "side")
    dens <- estimate_density_1d(data,  group, field)
    p <- ggplot(data=dens, na.rm=T) + parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(paste(capwords(fieldSP), collapse=" "))
    plot_field(p, paste(field, "den_sides","by_interface", sep="_"), grid=side~interface)
    
  }
  
  #Fractions
  field = "aromatic_dSASA_fraction"
  parts = list(plot_parts, scale_x_continuous("fraction", limit=c(0, 1.0)))
  group = c("sample_source", "side")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Aromatic dSASA Fraction")
  plot_field(p, paste(field, "den_sides", sep="_"), grid=side ~ .)
  
  group = c("sample_source", "interface", "side")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Aromatic dSASA Fraction")
  plot_field(p, paste(field, "den_sides","by_interface", sep="_"), grid=side~interface)
  
  #ScatterPlots
  parts = list(
    geom_point(size=1.5, pch="o"),
    #stat_smooth(color="grey"),
    stat_smooth(method=lm),
    geom_density2d(),
    #stat_density2d(aes(fill = ..level..), geom="polygon"),
    #stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE),
    theme_bw())
  
  #dSASA_hphobic vs dSASA_polar
  p <- ggplot(data=int_data, aes(x = dSASA_hphobic, y = dSASA_polar)) + parts +
    ggtitle("dSASA: Hydrophobic vs Polar") +
    scale_x_continuous("Hydrophobic SASA") +
    scale_y_continuous("Polar SASA")
  plot_field(p, "dSASA_hphobic_vs_polar", grid=sample_source ~ .)
  plot_field(p, "dSASA_hphobic_vs_polar_by_interface", grid=interface ~ sample_source)
  
  #dSASA vs dG
  p <- ggplot(data=int_data[int_data$dG<=5000 & int_data$dG>-5000,], aes(x = dSASA, y = dG)) + parts +
    ggtitle("dSASA vs dG") +
    scale_x_continuous("SASA") +
    scale_y_continuous("REU")
  plot_field(p, "dSASA_vs_dG", grid=sample_source ~ .)
  plot_field(p, "dSASA_vs_dG_by_interface", grid=interface ~ sample_source)
  
  #dSASA vs aromatic dSASA fraction.  Shouldn't increase, but worth a plot 
  p <- ggplot(data= data, aes(x = dSASA, y = aromatic_dSASA_fraction)) + parts +
    ggtitle("dSASA vs Aromatic dSASA Fraction") +
    scale_x_continuous("SASA") +
    scale_y_continuous("fraction", limit=c(0, 1.0))
  plot_field(p, "dSASA_vs_aromatic_dSASA_fraction", grid=side ~ sample_source)
  
  #aromatic dSASA fraction vs packstat
  
  #aromatic dSASA fraction vs sc_value
  
  #interface_nres vs dSASA.  Again, shouldn't be interesting.  Mainly a control.
  p <- ggplot(data = data, aes(x = interface_nres, y = dSASA)) + parts +
    ggtitle("Residues at Interface vs dSASA") +
    scale_x_continuous("n") +
    scale_y_continuous("SASA")
  plot_field(p, "interface_nres_vs_dSASA", grid=side ~ sample_source)
})) # end FeaturesAnalysis