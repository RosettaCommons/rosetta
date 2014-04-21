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
id = "int_composition-hbond_stats_by_restype",
author = "Jared Adolf-Bryfogle",
brief_description = "Graphs basic hbond densities for interface - interface hbonds",
feature_reporter_dependencies = c("InterfaceFeatures/AntibodyFeatures", "HBondFeatures","ResidueFeatures", "ResidueTypesFeatures", "PdbDataFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
  
  #Thanks to Matt O'Meara's help for this query -  Very slow:
  
  #########BROKEN############## 
  
  sele = "
  SELECT
    hb.energy as energy,
    don_res.interface as interface,
    don_restype.name3 as don_name3,
    acc_restype.name3 as acc_name3
  FROM
    interface_residues AS don_res,
    interface_residues AS acc_res,
    residues AS don_restype,
    residues AS acc_restype,
    hbond_sites AS don,
    hbond_sites AS acc,
    hbonds AS hb
  WHERE
    ((don_res.side == 'side1' AND
    acc_res.side == 'side2') OR
    (don_res.side =='side2' AND
    acc_res.side =='side1')) AND
    acc_res.interface == don_res.interface AND
    don.resNum == don_res.resNum == don_restype.resNum AND 
    acc.resNum == acc_res.resNum == acc_restype.resNum AND
    hb.don_id == don.site_id AND
    hb.acc_id == acc.site_id AND
    don_res.struct_id ==acc_res.struct_id == acc.struct_id AND
    acc.struct_id == don.struct_id AND
    don.struct_id == hb.struct_id AND
    don.struct_id == don_restype.struct_id AND
    don_restype.struct_id == acc_restype.struct_id
    
  "

  
  data = query_sample_sources(sample_sources, sele, char_as_factor=F)
  data$pair_name = paste(data$don_name3, data$acc_name3, sep="-")
  
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
  
  plot_field_wrap = function(p, plot_id, grid, columns = 4) {
    p <- p + facet_wrap(grid, ncol=columns)
    if(nrow(sample_sources) <= 3){
      p <- p + theme(legend.position="bottom", legend.direction="horizontal")
    }
    save_plots(self, plot_id, sample_sources, output_dir, output_formats)
  }
  
  
  #Hbond Energy density
  field = "energy"
  group = c("sample_source", "don_name3")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("REU") +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Interface hbond energies by donor")
  plot_field_wrap(p, paste("hbond", field, "den_by_donor_by_all", sep="_"), ~don_name3)
    
  group = c("sample_source", "acc_name3")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("REU") +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Interface hbond energies by donor")
  plot_field_wrap(p, paste("hbond", field, "den_by_acceptor_by_all", sep="_"), ~acc_name3)
  
  group = c("sample_source", "pair_name")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("REU") +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Interface hbond energies by pair")
  plot_field_wrap(p, paste("hbond", field, "den_by_don_acc_pair_by_all", sep="_"), ~pair_name)
  
  
  
  #Hbonds/model or per interface
  #There is probably a better way to do this.
  
  donor_counts <- ddply(data, .(interface, sample_source, struct_id, acc_name3), function(int_data){
    n = length(int_data$energy)
    df = data.frame(n=n)
  })
  
  acceptor_counts <- ddply(data, .(interface, sample_source, struct_id, don_name3), function(int_data){
    n = length(int_data$energy)
    df = data.frame(n=n)
  })
  
  pair_counts <- ddply(data, .(interface, sample_source, struct_id, pair_name), function(int_data){
    n = length(int_data$energy)
    df = data.frame(n=n)
  })
  
  
  field = "n"
  group = c("sample_source", "acc_name3")
  dens <- estimate_density_1d(donor_counts,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("n")
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Interface hbonds by acceptor")
  plot_field_wrap(p, paste("hbond_den_by_acceptor_by_all", sep="_"), ~acc_name3)
  
  group = c("sample_source", "don_name3")
  dens <- estimate_density_1d(acceptor_counts,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("n")
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Interface hbonds by donor")
  plot_field_wrap(p, paste("hbond_den_by_donor_by_all", sep="_"), ~don_name3)

  group = c("sample_source", "pair_name")
  dens <- estimate_density_1d(pair_counts,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("n")
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Interface hbonds by pair")
  plot_field_wrap(p, paste("hbond_den_by_don_acc_pair_by_all", sep="_"), ~pair_name)


  

  
})) # end FeaturesAnalysis