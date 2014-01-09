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
id = "interface_energies",
author = "Jared Adolf-Bryfogle",
brief_description = "Graphs basic interface energy information",
feature_reporter_dependencies = c("InterfaceFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
  
  sele <- "
  SELECT
    dG,
    dG_cross,
    delta_unsatHbonds,
    hbond_E_fraction,
    dSASA,
    interface
  FROM
    interfaces
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
  data_rm_out = data[data$dG<=5000 & data$dG>-5000,]#Remove high energy outliers
  
  #Basic Densities
  fields = c("dG", "dG_cross")
  for(field in fields){
    parts = list(plot_parts, scale_x_continuous("Rosetta Energy"))
    
    group = c("sample_source")
    dens <- estimate_density_1d(data_rm_out,  group, field)
    p <- ggplot(data=dens, na.rm=T) + parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(field)
    plot_field(p, paste(field, "den", sep="_"), )
    
    group = c("sample_source", "interface")
    dens <- estimate_density_1d(data_rm_out,  group, field)
    p <- ggplot(data=dens, na.rm=T) + parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(field)
    plot_field(p, paste(field, "den_by_interface", sep="_"), grid=interface ~ .)
  }
  
  field = "delta_unsatHbonds"
  parts = list(plot_parts, scale_x_continuous("n"))
  
  group = c("sample_source")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Delta Unsaturated Hbonds")
  plot_field(p, paste(field, "den", sep="_"), )
  
  group = c("sample_source", "interface")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Delta Unsaturated Hbonds")
  plot_field(p, paste(field, "den_by_interface", sep="_"), grid=interface ~ .)
  
  field = "hbond_E_fraction"
  parts = list(plot_parts, scale_x_continuous("fraction", limit=c(-.5, .5)))
  
  group = c("sample_source")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Hbond Interface Energy Fraction")
  plot_field(p, paste(field, "den", sep="_") )
  
  group = c("sample_source", "interface")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Hbond Interface Energy Fraction")
  plot_field(p, paste(field, "den_by_interface", sep="_"), grid=interface ~ .)
  
  parts = list(
    geom_point(size=1.5, pch="o"),
    #stat_smooth(color="grey"),
    stat_smooth(method=lm),
    geom_density2d(),
    #stat_density2d(aes(fill = ..level..), geom="polygon"),
    #stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE),
    theme_bw())
  
  #dSASA vs dG
  p <- ggplot(data=data_rm_out, aes(y = dSASA, x = dG)) + parts +
    ggtitle("dG vs dSASA") +
    ylab("SASA") +
    xlab("REU")
  plot_field(p, "dG_vs_dSASA", grid=sample_source ~ .)
  plot_field(p, "dG_vs_dSASA_by_interface", grid=interface ~ sample_source)
  
  #dG vs dG_cross
  p <- ggplot(data = data_rm_out, aes(x=dG, y=dG_cross)) + parts +
    ggtitle("dG vs dG_cross") +
    xlab("REU") +
    ylab("REU")
  plot_field(p, "dG_vs_dG_cross", grid=sample_source ~ .)
  plot_field(p, "dG_vs_dG_cross_by_interface", grid=interface ~ sample_source)
  
  p <- ggplot(data = data_rm_out, aes(x=dG, y=dG_cross, color=dSASA)) + parts +
    geom_point() +
    ggtitle("dG vs dG_cross") +
    xlab("REU") +
    ylab("REU") +
    scale_fill_hue(l=40)
  plot_field(p, "dG_vs_dG_cross_col_by_dSASA", grid=sample_source ~ .)
  plot_field(p, "dG_vs_dG_cross_col_by_dSASA_by_interface", grid=interface ~ sample_source) 
  
  #dG_cross vs dSASA
  p <- ggplot(data = data_rm_out, aes(x=dG_cross, y=dSASA)) + parts +
    ggtitle("dG_cross vs dSASA") +
    xlab("REU") +
    ylab("SASA")
  plot_field(p, "dg_cross_vs_dSASA", grid=sample_source ~ .)
  plot_field(p, "dG_cross_vs_dSASA_by_interface", grid=interface ~ sample_source)
  
  #hbond fraction vs dunsat hbonds
  p <- ggplot(data = data_rm_out, aes(x=hbond_E_fraction, y=delta_unsatHbonds)) + parts +
    ggtitle("HbondEfraction vs Delta unsat Hbonds") +
    xlab("fraction") +
    ylab("n")
  plot_field(p, "HbondEfraction_vs_deltaUnsatHbonds", grid=sample_source ~ .)
  plot_field(p, "HbondEfraction_vs_deltaUnsatHbonds_by_interface", grid=interface ~ sample_source)
  
  #dG vs dunsat hbonds
  p <- ggplot(data = data_rm_out, aes(x=dG, y=delta_unsatHbonds)) + parts +
    ggtitle("dG vs Delta unsat Hbonds") +
    xlab("REU") +
    ylab("n")
  plot_field(p, "dG_vs_deltaUnsatHbonds", grid=sample_source ~ .)
  plot_field(p, "dG_vs_deltaUnsatHbonds_by_interface", grid=interface ~ sample_source)
  
  #dSASA vs dunsathbonds
  #hbond fraction vs dunsat hbonds
  p <- ggplot(data = data_rm_out, aes(x=dSASA, y=delta_unsatHbonds)) + parts +
    ggtitle("dSASA vs Delta unsat Hbonds") +
    xlab("SASA") +
    ylab("n")
  plot_field(p, "dSASA_vs_deltaUnsatHbonds", grid=sample_source ~ .)
  plot_field(p, "dSASA_vs_deltaUnsatHbonds_by_interface", grid=interface ~ sample_source)
  
  #dG_cross vs dunsathbonds
  #hbond fraction vs dunsat hbonds
  p <- ggplot(data = data_rm_out, aes(x=dG_cross, y=delta_unsatHbonds)) + parts +
    ggtitle("dG_cross vs Delta unsat Hbonds") +
    xlab("REU") +
    ylab("n")
  plot_field(p, "dG_cross_vs_deltaUnsatHbonds", grid=sample_source ~ .)
  plot_field(p, "dG_cross_vs_deltaUnsatHbonds_by_interface", grid=interface ~ sample_source)
  
  
})) # end FeaturesAnalysis