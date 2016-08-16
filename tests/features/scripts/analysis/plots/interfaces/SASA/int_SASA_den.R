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
id = "int_SASA_den",
author = "Jared Adolf-Bryfogle",
brief_description = "Graphs basic dSASA and SASA information",
feature_reporter_dependencies = c("InterfaceFeatures/AntibodyFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
  
  #First we run on all the interfaces in the database

  

#  sele = "
#  SELECT
#    dSASA,
#    dSASA_hphobic,
#    dSASA_polar,
#    interface
#  FROM
#    interfaces"
#  
#  data = query_sample_sources(sample_sources, sele, char_as_factor=F)
# 

#  fields = c("dSASA", "dSASA_hphobic", "dSASA_polar")
#  for(field in fields){
#    
#    group = c("sample_source")
#    dens <- estimate_density_1d(data,  group, field)
#    p <- ggplot(data=dens, na.rm=T) + parts +
#      geom_line(aes(x, y, colour=sample_source), size=1.2) +
#      ggtitle(field)
#    plot_field(p, paste(field, "den_by_all", sep="_"))
#    
#    group = c("sample_source", "interface")
#    dens <- estimate_density_1d(data,  group, field)
#    p <- ggplot(data=dens, na.rm=T) + parts +
#      geom_line(aes(x, y, colour=sample_source), size=1.2) +
#      ggtitle(field)
#    plot_field(p, paste(field, "den_by_interface", sep="_"),grid=~interface)
#  }
#
#  
#  #dSASA sides
#  int_data = data
  
  plot_parts <- list(
    geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)),
    scale_y_continuous("Feature Density"),
    theme_bw())

 plot_field = function(p, plot_id, grid = NULL){
    
    if (! is.null(grid)){
      p <- p+ facet_grid(facets=grid)
    }
    if(nrow(sample_sources) <= 3){
      p <- p + theme(legend.position="bottom", legend.direction="horizontal")
    }
    save_plots(self, plot_id, sample_sources, output_dir, output_formats)
  }
  
  sele = "
  SELECT
    dSASA,
    dSASA_sc,
    dSASA - dSASA_sc as dSASA_bb,
    dhSASA,
    dhSASA_sc,
    dhSASA - dhSASA_sc as dhSASA_bb,
    dhSASA_rel_by_charge,
    aromatic_dSASA_fraction,
    interface,
    side
  FROM
    interface_sides
  ORDER BY dSASA DESC
  "
  
  #Polar fraction - from Ben Strange's Paper
  

  parts = list(plot_parts, xlab("SASA"))
  
  data = query_sample_sources(sample_sources, sele, char_as_factor=F)
  
  data$polar_fraction = (data$dSASA - data$dhSASA)/data$dSASA
  field = "polar_fraction"
  
  group = c("sample_source", "side")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    xlab("dSASA(Polar)/dSASA") +
    ggtitle("Polar dSASA Fraction")
  plot_field(p, "dSASA_polar_fraction_den_by_all", grid=side ~ .)
  
  group = c("sample_source", "interface", "side")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    xlab("dSASA(Polar)/dSASA") +
    ggtitle("Polar dSASA Fraction")
  plot_field(p, "dSASA_polar_fraction_den_by_interface", grid=side~interface)
  
  #print(data)
  
  #Backbone SASA may not be interesting, but I want I still want to know for now.
  fields = c("dSASA", "dSASA_bb", "dSASA_sc", "dhSASA", "dhSASA_bb", "dhSASA_sc", "dhSASA_rel_by_charge")
  for (field in fields){

    parts = list(plot_parts, scale_x_continuous("SASA"))
    group = c("sample_source", "side")
    dens <- estimate_density_1d(data,  group, field)
    p <- ggplot(data=dens, na.rm=T) + parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(paste("Buried", field, sep=" "))
    plot_field(p, paste(field, "den_sides_by_all", sep="_"), grid=side ~ .)
  
    group = c("sample_source", "interface", "side")
    dens <- estimate_density_1d(data,  group, field)
    p <- ggplot(data=dens, na.rm=T) + parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(paste("Buried", field, sep=" "))
    plot_field(p, paste(field, "den_sides","by_interface", sep="_"), grid=side~interface)
    
  
  }
  

#  Plotting all together - Might look like crap, but lets try it.
#  group = c("sample_source", "interface", "side")
#  dens_dsasa <- estimate_density_1d(data,  group, c("dSASA"))
#  dens_dsasa_bb <- estimate_density_1d(data, group, c("dSASA_bb"))
#  dens_dsasa_sc <- estimate_density_1d(data, group, c("dSASA_sc"))
#  
#  p <- ggplot(data=dens_dsasa, na.rm=T) + parts +
#    geom_line(aes(x, y, colour=sample_source), size=1.2) +
#    #geom_point(data=dens_dsasa, aes(x, y, colour=sample_source, size=.5, pch="o")) +
#    geom_line(data=dens_dsasa_bb, aes(x, y, colour=sample_source, linetype= "dotted"), size=1.2) +
#    geom_line(data=dens_dsasa_sc, aes(x, y, colour=sample_source, linetype= "dotdash"), size= 1.2) +
#    ggtitle("dSASA Density")
      
#  plot_field(p, paste("dSASA_all", "den_sides","by_interface", sep="_"), grid=side~interface)
  
  ####  Means  #########
  fields = c("dSASA", "dhSASA")
  for (field in fields){

  avgs <- ddply(data, .(sample_source, side, field), function(d2){
    data.frame(m = mean(d2[,field]), std_dev = sd(d2[,field]), m_top10 = mean(d2[1:10,field]), std_dev_top_10 = sd(d2[1:10,field]), top = d2[1,field])
  })
    
  p <- ggplot(data=avgs ) +
    geom_bar(position="dodge", stat='identity', aes(x = sample_source, y= m , fill=sample_source)) +
    theme_bw() +
    ggtitle(paste("Buried", field,"Average", sep=" "))+
    scale_x_discrete(labels = abbreviate)
  plot_field(p, "avg_sides_by_all", grid=side ~ .)
  
  #Average Top 10
  p <- ggplot(data=avgs ) + 
    geom_bar(position="dodge", stat='identity', aes(x = sample_source, y= m_top10 , fill=sample_source)) +
    theme_bw() +
    ggtitle(paste("Buried", field, "Average Best 10",sep=" ")) +
    scale_x_discrete(labels = abbreviate)
  plot_field(p, "avg_sides_top_10_by_all", grid=side ~ .)
  
  #Best
  p <- ggplot(data=avgs ) + 
    geom_bar(position="dodge", stat='identity', aes(x = sample_source, y= top , fill=sample_source)) +
    theme_bw() +
    ggtitle(paste("Buried", field, "top", sep=" ")) +
    scale_x_discrete(labels = abbreviate)
  plot_field(p, "sides_top_by_all", grid=side ~ .)
  
  avgs <- ddply(data, .(sample_source, side, field, interface), function(d2){
    data.frame(m = mean(d2[,field]), std_dev = sd(d2[,field]), m_top10 = mean(d2[1:10,field]), std_dev_top_10 = sd(d2[1:10,field]), top = d2[1,field])
  })
      
  p <- ggplot(data=avgs ) + 
    geom_bar(position="dodge", stat='identity', aes(x = sample_source, y= m , fill=sample_source)) +
    theme_bw() +
    ggtitle(paste("Buried", field,"Average", sep=" ")) +
    scale_x_discrete(labels = abbreviate)
  plot_field(p, "avg_sides_by_interface", grid=side ~ interface)
  
  #Average Top 10
  p <- ggplot(data=avgs ) + 
    geom_bar(position="dodge", stat='identity', aes(x = sample_source, y= m_top10 , fill=sample_source)) +
    theme_bw() +
    ggtitle(paste("Buried", field, "Average Best 10",sep=" ")) +
    scale_x_discrete(labels = abbreviate) +
    ylab(field)
  plot_field(p, "avg_sides_top_10_by_interface", grid=side ~ interface)
  
  #Best
  p <- ggplot(data=avgs ) + 
    geom_bar(position="dodge", stat='identity', aes(x = sample_source, y= top , fill=sample_source)) +
    theme_bw() +
    ggtitle(paste("Buried", field, "top", sep=" ")) +
    scale_x_discrete(labels = abbreviate) +
    ylab(field)
  plot_field(p, "sides_top_by_interface", grid=side ~ interface)
    
  }
  #Fractions
  field = "aromatic_dSASA_fraction"
  parts = list(plot_parts, scale_x_continuous("fraction", limit=c(0, 1.0)))
  group = c("sample_source", "side")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Aromatic dSASA Fraction") +
  plot_field(p, "dSASA_aromatic_fraction_den_by_all", grid=side ~ .)
  
  group = c("sample_source", "interface", "side")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Aromatic dSASA Fraction")
  plot_field(p, "dSASA_aromatic_fraction_den_by_interface", grid=side~interface)
  

  
})) # end FeaturesAnalysis