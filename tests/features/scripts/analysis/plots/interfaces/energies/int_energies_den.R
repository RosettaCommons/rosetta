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
id = "int_energies_den",
author = "Jared Adolf-Bryfogle",
brief_description = "Graphs basic interface energy information",
feature_reporter_dependencies = c("InterfaceFeatures/AntibodyFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
  
  sele <- "
  SELECT
    interfaces.dG as dG,
    interfaces.dG_cross as dG_cross,
    interfaces.hbond_E_fraction as hbond_E_fraction,
    interfaces.interface as interface,
    structure_scores.score_value as total_score
  FROM
    interfaces,
    score_types,
    structure_scores
  WHERE
    score_types.score_type_name='total_score' AND
    structure_scores.score_type_id = score_types.score_type_id AND
    structure_scores.struct_id = interfaces.struct_id 
  ORDER BY dG;
  "

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
  
  
  data = query_sample_sources(sample_sources, sele)
  #data_rm_out = data[data$dG<=5000 & data$dG>-5000,]#Remove high energy outliers
  
  #Basic Densities
  fields = c("dG", "dG_cross")
  for(field in fields){
    parts = list(plot_parts, scale_x_continuous("Rosetta Energy"))
    
    group = c("sample_source")
    dens <- estimate_density_1d(data,  group, field)
    p <- ggplot(data=dens, na.rm=T) + parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(field)
    plot_field(p, paste(field, "den_by_all", sep="_"), )
    
    group = c("sample_source", "interface")
    dens <- estimate_density_1d(data,  group, field)
    p <- ggplot(data=dens, na.rm=T) + parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(field)
    plot_field(p, paste(field, "den_by_interface", sep="_"), grid=interface ~ .)
  }
  
  
  field = "hbond_E_fraction"
  parts = list(plot_parts, xlab("fraction"))
  
  group = c("sample_source")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Hbond Interface Energy Fraction")
  plot_field(p, paste(field, "den_by_all", sep="_") )
  
  group = c("sample_source", "interface")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Hbond Interface Energy Fraction")
  plot_field(p, paste(field, "den_by_interface", sep="_"), grid=interface ~ .)
  
  #Side energies?
  
  
  #Averages:
  avgs <- ddply(data, .(sample_source, interface), function(d2){
    data.frame(m = mean(d2$dG), std_dev = sd(d2$dG), m_top10 = mean(d2[1:10,]$dG), std_dev_top_10 = sd(d2[1:10,]$dG), top = d2[1,]$dG)
  })
  p <- ggplot(data=avgs ) + 
    geom_bar(position="dodge", stat='identity', aes(x = sample_source, y= m , fill=sample_source)) +
    theme_bw() +
    ggtitle("Average Interface dG") +
    ylab("REU") +
    scale_x_discrete(labels = abbreviate)
  plot_field(p, "avg_dG_by_interface", grid=interface ~ .)
  
  #Average Top 10
  p <- ggplot(data=avgs ) + 
    geom_bar(position="dodge", stat='identity', aes(x = sample_source, y= m_top10 , fill=sample_source)) +
    theme_bw() +
    ggtitle("Average Best 10 Interface dG") +
    ylab("REU") +
    scale_x_discrete(labels = abbreviate)
  plot_field(p, "avg_dG_top_10_by_interface", grid=interface ~ .)
  
  #Best
  p <- ggplot(data=avgs ) + 
    geom_bar(position="dodge", stat='identity', aes(x = sample_source, y= top , fill=sample_source)) +
    theme_bw() +
    ggtitle("Top Interface dG") +
    ylab("REU") +
    scale_x_discrete(labels = abbreviate)
  plot_field(p, "dG_top_by_interface", grid=interface ~ .)
})) # end FeaturesAnalysis