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
id = "interface_residues",
author = "Jared Adolf-Bryfogle",
brief_description = "Graphs all information of individual interface residues.
Should be same interface, same numbering scheme / decoy set for this to have any meaning.",
feature_reporter_dependencies = c("InterfaceFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

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
  
  sele <-"
  SELECT
    interface_residues.interface as interface,
    interface_residues.relative_dSASA_fraction as dSASA_fraction,
    interface_residues.dSASA as dSASA,
    interface_residues.dG as dG,
    interface_residues.energy_int as energy_int,
    interface_residues.energy_sep as energy_sep
  FROM
    interface_residues"
  
  #Density plots
  
  data = query_sample_sources(sample_sources, sele, char_as_factor=F)
  ##Overall plots for all residues: Add Side data once we have this.
  
  #Densities
  
  #dSASA
  field = "dSASA"
  group = c("sample_source")
  dens <- estimate_density_1d(data[data$dSASA > 0,],  group, field)
  p <- ggplot(data = dens, na.rm=T) + plot_parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle(paste("Residue", field, sep=" ")) +
    xlab(field)
  plot_field(p, paste("residue", field, "dens", sep="_"))
    
  group = c("sample_source", "interface")
  dens <- estimate_density_1d(data[data$dSASA > 0,],  group, field)
  p <- ggplot(data = dens, na.rm=T) + plot_parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle(paste("Residue", field, sep=" ")) +
    xlab(field)
  plot_field(p, paste("residue", field, "dens_by_interface", sep="_"), grid=interface ~ .)
  
  #dSASA fraction
  field = "dSASA_fraction"
  group = c("sample_source")
  dens <- estimate_density_1d(data[data$dSASA > 0,],  group, field)
  p <- ggplot(data = dens, na.rm=T) + plot_parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle(paste("Residue", field, sep=" ")) +
    scale_x_continuous(field, limit = c(0, 1.0))
  plot_field(p, paste("residue", field, "dens", sep="_"))
  
  group = c("sample_source", "interface")
  dens <- estimate_density_1d(data[data$dSASA > 0,],  group, field)
  p <- ggplot(data = dens, na.rm=T) + plot_parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle(paste("Residue", field, sep=" ")) +
    scale_x_continuous(field, limit = c(0, 1.0))
  plot_field(p, paste("residue", field, "dens_by_interface", sep="_"), grid=interface ~ .)
  
  #Energies
  fields = c("dG", "energy_int", "energy_sep")
  for(field in fields){
    group = c("sample_source")
    dens <- estimate_density_1d(data[-15 < data[field] & data[field] < 15,],  group, field)
    p <- ggplot(data = dens, na.rm=T) + plot_parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(paste("Residue", field, sep=" ")) +
      scale_x_continuous("REU", limit = c(-15, 15))
    plot_field(p, paste("residue", field, "dens", sep="_"))
  
    group = c("sample_source", "interface")
    dens <- estimate_density_1d(data[-15 <data[field] & data[field] < 15,],  group, field)
    p <- ggplot(data = dens, na.rm=T) + plot_parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(paste("Residue", field, sep=" ")) +
      scale_x_continuous("REU", limit = c(-15, 15))
    plot_field(p, paste("residue", field, "dens_by_interface", sep="_"), grid=interface ~ .)
  }
  #dG where dSASA is 0:
  field = "dG"
  group = c("sample_source")
  dens <- estimate_density_1d(data[-15 < data[field] & data[field] < 15,],  group, field)
  p <- ggplot(data = dens, na.rm=T) + plot_parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle(paste("Residue", field, sep=" ")) +
    scale_x_continuous(field, limit=c(-15, 15))
  plot_field(p, paste("residue", field, "@0dSASA_dens", sep="_"))
  
  group = c("sample_source", "interface")
  dens <- estimate_density_1d(data[-15 < data[field] & data[field] < 15,],  group, field)
  p <- ggplot(data = dens, na.rm=T) + plot_parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle(paste("Residue", field, sep=" ")) +
    scale_x_continuous(field, limit=c(-15, 15))
  plot_field(p, paste("residue", field, "@0dSASA_dens_by_interface", sep="_"), grid=interface ~ .)
  
  #Scatterplots
  
  #dSASA vs dSASA fraction
  parts = list(
    geom_point(size=.5, pch="o"),
    stat_smooth(method = lm),
    geom_density2d(),
    scale_x_continuous("dSASA"),
    scale_y_continuous("dSASA fraction", limit = c(0, 1.0)),
    theme_bw(),
    ggtitle("Residue dSASA vs dSASA Fraction"))

  p <- ggplot(data=data[data$dSASA > 0 & data$dSASA_fraction > 0,], aes(x=dSASA, y=dSASA_fraction)) +
    parts
  plot_field(p, "dSASA_vs_dSASA_fraction", grid = sample_source ~ .)
  
  p <- ggplot(data=data[data$dSASA > 0 & data$dSASA_fraction >0,], aes(x=dSASA, y=dSASA_fraction)) +
    parts
  plot_field(p, "dSASA_vs_dSASA_fraction_by_interface", grid = sample_source ~ interface)
  
  #->ss_overlay functions complain cannot coerce type 'symbol' to vector of type 'double' for some reason
  
  p <- ggplot(data=data[data$dSASA > 0 & data$dSASA_fraction > 0,], aes(x=dSASA, y=dSASA_fraction, color=factor(sample_source))) +
    parts
  plot_field(p, "dSASA_vs_dSASA_fraction_ss_overlay")

  p <- ggplot(data=data[data$dSASA > 0 & data$dSASA_fraction > 0,], aes(x=dSASA, y=dSASA_fraction, color=factor(sample_source))) +
    parts
  plot_field(p, "dSASA_vs_dSASA_fraction_by_interface_ss_overlay", grid=interface ~ .)
  
  #dSASA vs dG
  parts = list(
    xlab("dSASA"),
    scale_y_continuous("REU", limit = c(-15, 15)),
    geom_point(size=.5, pch="o"),
    stat_smooth(method = lm),
    stat_density2d(),
    theme_bw(),
    ggtitle("Residue dSASA vs dG"))
  
  p <- ggplot(data=data[data$dSASA > 0 & -15 < data[field] & data[field] < 15,], aes(x=dSASA, y=dG)) +
    parts
  plot_field(p, "dSASA_vs_dG", grid = sample_source ~ .)
  
  p <- ggplot(data=data[data$dSASA > 0 & -15 < data[field] & data[field] < 15,], aes(x=dSASA, y=dG)) +
    parts
  plot_field(p, "dSASA_vs_dG_by_interface", grid = sample_source ~ interface)

  p <- ggplot(data=data[data$dSASA > 0 & -15 < data[field] & data[field] < 15,], aes(x=dSASA, y=dG, color=sample_source)) +
    parts
  plot_field(p, "dSASA_vs_dG_ss_overlay")
  
  p <- ggplot(data=data[data$dSASA > 0 & -15 < data[field] & data[field] < 15,], aes(x=dSASA, y=dG, color=sample_source)) +
    parts
  plot_field(p, "dSASA_vs_dG_by_interface_ss_overlay", grid=interface ~ .)
  
  #dG vs dSASA_fraction
  parts = list(
    scale_x_continuous("dSASA fraction", limit = c(0, 1.0)),
    scale_y_continuous("REU", limit=c(-15, 15)),
    geom_point(size=.5, pch="o"),
    stat_smooth(method = lm),
    stat_density2d(),
    theme_bw(),
    ggtitle("Residue dSASA fraction vs dG"))
  
  p <- ggplot(data=data[data$dSASA_fraction > 0 & -10 < data[field] & data[field] < 15,], aes(y=dG, x=dSASA_fraction, color=sample_source)) +
    parts
  plot_field(p, "dSASA_fraction_vs_dG_ss_overlay")
  
  p <- ggplot(data=data[data$dSASA_fraction > 0 & -15 < data[field] & data[field] < 15,], aes(y=dG, x=dSASA_fraction, color=sample_source)) +
    parts
  plot_field(p, "dSASA_fraction_vs_dG_by_interface_ss_overlay", grid=interface ~ .)
  
  p <- ggplot(data=data[data$dSASA_fraction > 0 & -15 < data[field] & data[field] < 15,], aes(y=dG, x=dSASA_fraction)) +
    parts
  plot_field(p, "dSASA_fraction_vs_dG", grid = sample_source ~ .)
  
  p <- ggplot(data=data[data$dSASA_fraction > 0 & -15 < data[field] & data[field] < 15,], aes(y=dG, x=dSASA_fraction)) +
    parts
  plot_field(p, "dSASA_fraction_vs_dG_by_interface", grid = sample_source ~ interface)
  
  #Per residue data.  This may get crazy.
}))