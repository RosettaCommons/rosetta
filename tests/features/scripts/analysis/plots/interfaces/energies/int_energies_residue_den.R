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
id = "int_energies-by_residue_den",
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
    interface_residues.dG as dG,
    interface_residues.dSASA as dSASA,
    interface_residues.energy_int as energy_int,
    interface_residues.energy_sep as energy_sep
  FROM
    interface_residues"
  
  #Density plots
  
  data = query_sample_sources(sample_sources, sele, char_as_factor=F)
  ##Overall plots for all residues: Add Side data once we have this.
  
  #Densities
  
  
  #Energies
  fields = c("dG", "energy_int", "energy_sep")
  for(field in fields){
    group = c("sample_source")
    dens <- estimate_density_1d(data,  group, field)
    p <- ggplot(data = dens, na.rm=T) + plot_parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(paste("Residue", field, sep=" ")) +
      xlab("REU")
      #scale_x_continuous("REU", limit = c(-15, 15))
    plot_field(p, paste(field, "residue_dens_by_all", sep="_"))
  
    group = c("sample_source", "interface")
    dens <- estimate_density_1d(data,  group, field)
    p <- ggplot(data = dens, na.rm=T) + plot_parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(paste("Residue", field, sep=" ")) +
      xlab("REU")
      #scale_x_continuous("REU", limit = c(-15, 15))
    plot_field(p, paste(field, "residue_dens_by_interface", sep="_"), grid=interface ~ .)
  }
  
  
  #dG where dSASA is 0:
  #data[-15 < data[field] & data[field] < 15,]
  
  field = "dG"
  group = c("sample_source")
  dens <- estimate_density_1d(data[data$dSASA == 0,],  group, field)
  p <- ggplot(data = dens, na.rm=T) + plot_parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle(paste("Residue", field, sep=" ")) +
    xlab("REU")
    #scale_x_continuous(field, limit=c(-15, 15))
  plot_field(p, paste(field, "residue_@0dSASA_dens", sep="_"))
  
  group = c("sample_source", "interface")
  dens <- estimate_density_1d(data[data$dSASA == 0,],  group, field)
  p <- ggplot(data = dens, na.rm=T) + plot_parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle(paste("Residue", field, sep=" ")) +
    xlab("REU")
    #scale_x_continuous(field, limit=c(-15, 15))
  plot_field(p, paste(field, "residue_@0dSASA_dens_by_interface", sep="_"), grid=interface ~ .)
  
  #dG where dSASA > 0 :
  field = "dG"
  group = c("sample_source")
  dens <- estimate_density_1d(data[data$dSASA >0,],  group, field)
  p <- ggplot(data = dens, na.rm=T) + plot_parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle(paste("Residue", field, sep=" ")) +
    xlab("REU")
    #scale_x_continuous(field, limit=c(-15, 15))
  plot_field(p, paste(field, "residue_>0dSASA_dens", sep="_"))
  
  group = c("sample_source", "interface")
  dens <- estimate_density_1d(data[data$dSASA > 0,],  group, field)
  p <- ggplot(data = dens, na.rm=T) + plot_parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle(paste("Residue", field, sep=" ")) +
    xlab("REU")
    #scale_x_continuous(field, limit=c(-15, 15))
  plot_field(p, paste(field, "residue_>0dSASA_dens_by_interface", sep="_"), grid=interface ~ .)
  #Per residue data.  This may get crazy.
}))