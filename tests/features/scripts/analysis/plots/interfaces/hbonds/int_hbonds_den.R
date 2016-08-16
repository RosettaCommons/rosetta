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
id = "int_hbonds_den",
author = "Jared Adolf-Bryfogle",
brief_description = "Graphs basic hbond densities for interface - interface hbonds",
feature_reporter_dependencies = c("InterfaceFeatures/AntibodyFeatures", "HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
  
  #Thanks to Matt O'Meara's help for this query -  Very slow:
  

  sele = "
    SELECT
    hb.energy as energy,
    don_res.interface as interface,
    don_res.struct_id as struct_id,
    hb_geom.AHdist as dis
  FROM
    interface_residues AS don_res,
    interface_residues AS acc_res,
    hbond_sites AS don,
    hbond_sites AS acc,
    hbonds AS hb,
    hbond_geom_coords as hb_geom
  WHERE
    ((don_res.side== 'side1' AND
    acc_res.side == 'side2') OR
    (don_res.side=='side2' AND
    acc_res.side=='side1')) AND
    acc_res.interface == don_res.interface AND
    don.resNum == don_res.resNum AND 
    acc.resNum == acc_res.resNum AND
    hb.don_id == don.site_id AND
    hb.acc_id == acc.site_id AND
    hb.hbond_id == hb_geom.hbond_id AND
    don_res.struct_id == acc_res.struct_id AND
    acc_res.struct_id == acc.struct_id AND
    acc.struct_id == don.struct_id AND
    don.struct_id == hb.struct_id AND
    hb.struct_id == hb_geom.struct_id
  "
  
  data = query_sample_sources(sample_sources, sele, char_as_factor=F)
  #print(sum(data$struct_id==1))
  #print(sum(data$struct_id==2))
  #print(sum(data$struct_id==3))
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
  

  #Hbond Energy density
  field = "energy"
  group = c("sample_source")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("REU") +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Cross Interface Hydrogen Bond Energies")
  plot_field(p, paste("hbond", field, "den_by_all", sep="_"))
    
  group = c("sample_source", "interface")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("REU") +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Cross Interface Hydrogen Bond Energies")
  plot_field(p, paste("hbond", field, "den_by_interface", sep="_"),grid=~interface)
  
  #Hbond Distances
  field = "dis"
  group = c("sample_source")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("Angstroms") +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Cross Interface Hydrogen Bond Distances")
  plot_field(p, paste("hbond", field, "den_by_all", sep="_"))
    
  group = c("sample_source", "interface")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("Angstroms") +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Cross Interface Hydrogen Bond Distances")
  plot_field(p, paste("hbond", field, "den_by_interface", sep="_"),grid=~interface)
  
  #Hbonds/model or per interface
  #There is probably a better way to do this.
  
  hbond_counts <- ddply(data, .(interface, sample_source, struct_id), function(int_data){
    n = length(int_data$energy > 0)
    df = data.frame(n=n)
  })
  
  #print(head(hbond_counts))
  field = "n"
  group = c("sample_source")
  dens <- estimate_density_1d(hbond_counts,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("hbonds") +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Cross Interface Hydrogen Bonds")
  plot_field(p, "hbond_den_by_all")
    
  group = c("sample_source", "interface")
  dens <- estimate_density_1d(hbond_counts,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("hbonds") +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Cross Interface Hydrogen Bonds")
  plot_field(p, "hbond_den_by_interface",grid=~interface)
  
  
  #Histogram
  p <- ggplot(data=hbond_counts, na.rm=T) + 
    geom_bar(aes(x=n, y = ..density.. , fill=sample_source), position="dodge", binwidth=1) +
    scale_y_continuous(label=percent) +
    xlab("hbonds") +
    ggtitle("Average Cross Interface Hydrogen Bonds")
  plot_field(p, "hbond_hist_by_all")
  
  p <- ggplot(data=hbond_counts, na.rm=T) + 
    geom_bar(aes(x=n, y = ..density.. , fill=sample_source), position="dodge", binwidth=1) +
    scale_y_continuous(label=percent) +
    xlab("hbonds") +
    ggtitle("Average Cross Interface Hydrogen Bonds")
  plot_field(p, "hbond_hist_by_interface", grid=~interface)
})) # end FeaturesAnalysis