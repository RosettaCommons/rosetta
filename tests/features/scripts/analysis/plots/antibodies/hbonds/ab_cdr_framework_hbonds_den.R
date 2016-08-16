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
id = "ab_hbonds-cdr_framework_den",
author = "Jared Adolf-Bryfogle",
brief_description = "CDR -> Framework Hbonds",
feature_reporter_dependencies = c("AntibodyFeatures", "HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
  
  
  sele_don = "
SELECT
    DISTINCT
    hb.energy as energy,
    don.struct_id as struct_id,
    don_pdb_info.residue_number as resnum,
    don_pdb_info.chain_id as chainid1,
    acc_pdb_info.residue_number as resnum2,
    acc_pdb_info.chain_id as chainid2,
    hb_geom.AHdist as distance,
    don_cdr_res.CDR
  FROM
    interface_residues as acc_res,
    hbond_sites AS don,
    hbond_sites AS acc,
    hbonds AS hb,
    hbond_geom_coords as hb_geom,
    cdr_residues as don_cdr_res,
    cdr_residues as acc_cdr_res,
    residue_pdb_identification as acc_pdb_info,
    residue_pdb_identification as don_pdb_info
  WHERE
    acc.struct_id == don.struct_id AND
    don.struct_id == hb.struct_id AND
    hb.struct_id == hb_geom.struct_id AND
    hb.struct_id == don_cdr_res.struct_id == acc_cdr_res.struct_id AND
    acc.struct_id == acc_pdb_info.struct_id AND
    acc.struct_id == don_pdb_info.struct_id AND
    don.resNum == don_cdr_res.resNum AND
    acc.resNum == acc_pdb_info.residue_number AND
    don.resNum == don_pdb_info.residue_number AND
    acc.resNum != acc_cdr_res.resNum AND
    (acc_pdb_info.chain_id == 'L' OR acc_pdb_info.chain_id == 'H') AND
    hb.don_id == don.site_id AND
    hb.acc_id == acc.site_id AND
    hb.hbond_id == hb_geom.hbond_id
  "
  
  sele_acc = "
SELECT
    hb.energy as energy,
    acc.struct_id as struct_id,
    acc.resNum as resnum,
    hb_geom.AHdist as distance,
    cdr_residues.CDR
  FROM
    interface_residues as don_res,
    hbond_sites AS don,
    hbond_sites AS acc,
    hbonds AS hb,
    hbond_geom_coords as hb_geom,
    cdr_residues
  WHERE
    hb.don_id == don.site_id AND
    hb.acc_id == acc.site_id AND
    hb.hbond_id == hb_geom.hbond_id AND
    acc.struct_id == don.struct_id AND
    don.struct_id == hb.struct_id AND
    hb.struct_id == hb_geom.struct_id AND
    hb.struct_id == cdr_residues.struct_id AND
    acc.resNum == cdr_residues.resNum AND
    don.resNum == don_res.resNum AND
    don.struct_id == don_res.struct_id AND
    don_res.side == 'side2'
  "
  
  don_data = query_sample_sources(sample_sources, sele_don, char_as_factor=F)
  acc_data = query_sample_sources(sample_sources, sele_acc, char_as_factor=F)
  
  
  #print(sum(data$struct_id==1))
  #print(sum(data$struct_id==2))
  #print(sum(data$struct_id==3))
  plot_parts <- list(
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
  
  don_data$type = "don"
  acc_data$type = "acc"
  data = rbind(don_data, acc_data)
  
  #Hbond Energy density
  field = "energy"
  group = c("sample_source", "CDR")
  dens <- estimate_density_1d(don_data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("REU") +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Cross Ab-Ag Hydrogen Bond Energies, CDR Donor")
  plot_field(p, paste("hbond", field, "don_den_by_cdr", sep="_"), ~ CDR)
  
  dens <- estimate_density_1d(acc_data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("REU") +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Cross Ab-Ag Hydrogen Bond Energies, CDR Acceptor")
  plot_field(p, paste("hbond", field, "acc_den_by_cdr", sep="_"), ~ CDR)
  
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("REU") +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Cross Ab-Ag Hydrogen Bond Energies")
  plot_field(p, paste("hbond", field, "tog_den_by_cdr", sep="_"), ~ CDR)
  
  #Hbond Distances
  field = "distance"
  group = c("sample_source", "CDR")
  dens <- estimate_density_1d(don_data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("Angstroms") +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Cross Ab-Ag Hydrogen Bond Distances, CDR Donor")
  plot_field(p, paste("hbond", field, "don_den_by_cdr", sep="_"), ~ CDR)
  
  dens <- estimate_density_1d(acc_data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("Angstroms") +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Cross Ab-Ag Hydrogen Bond Distances, CDR Acceptor")
  plot_field(p, paste("hbond", field, "acc_den_by_cdr", sep="_"), ~ CDR)
  
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("Angstroms") +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Cross Ab-Ag Hydrogen Bond Distances")
  plot_field(p, paste("hbond", field, "tog_den_by_cdr", sep="_"), ~ CDR)
  
  #Hbond Counts per CDR
  don_counts <- ddply(don_data, .(sample_source, struct_id, CDR), function(data){
    n = length(data$energy > 0)
    df = data.frame(n=n)
  })
  
  acc_counts <- ddply(acc_data, .(sample_source, struct_id, CDR), function(data){
    n = length(data$energy > 0)
    df = data.frame(n=n)
  })
  
  counts <- ddply(data, .(sample_source, struct_id, CDR), function(data){
    n = length(data$energy > 0)
    df = data.frame(n=n)
  })
  
  #print(head(hbond_counts))
  field = "n"
  group = c("sample_source", "CDR")
  
  dens <- estimate_density_1d(don_counts,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("hbonds") +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Cross Ab-Ag Hydrogen Bonds, CDR Donor")
  plot_field(p, "hbond_don_den_by_cdr", ~ CDR)

  dens <- estimate_density_1d(acc_counts,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("hbonds") +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Cross Ab-Ag Hydrogen Bonds, CDR Acceptor")
  plot_field(p, "hbond_acc_den_by_cdr", ~ CDR)
  
  dens <- estimate_density_1d(counts,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("hbonds") +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Cross Interface hbonds")
  plot_field(p, "hbond_tog_den_by_cdr", ~ CDR)
    
  don_avgs <- ddply(don_data, .(sample_source, CDR), function(data){
    data.frame(sample_source=data$sample_source, CDR=data$CDR, m=mean(data$n))
  })
  
  acc_avgs <- ddply(acc_data, .(sample_source, struct_id, CDR), function(data){
    data.frame(sample_source=data$sample_source, CDR=data$CDR, m=mean(data$n))
  })
  
  avgs <- ddply(data, .(sample_source, struct_id, CDR), function(data){
    data.frame(sample_source=data$sample_source, CDR=data$CDR, m=mean(data$n))
  })
  
  #Histograms
  p <- ggplot(data=don_counts, na.rm=T) + 
    geom_bar(aes(x=CDR, y = m, fill=sample_source), position="dodge", stat='identity') +
    scale_y_continuous(label=percent) +
    xlab("hbonds") +
    ggtitle("Cross Interface Hydrogen Bonds, CDR Donor")
  plot_field(p, "hbond_don_hist_by_cdr")

  p <- ggplot(data=acc_counts, na.rm=T) + 
    geom_bar(aes(x=CDR, y = m , fill=sample_source), position="dodge", stat='identity') +
    scale_y_continuous(label=percent) +
    xlab("hbonds") +
    ggtitle("Cross Interface Hydrogen Bonds, CDR Acceptor")
  plot_field(p, "hbond_acc_hist_by_cdr")
  
  p <- ggplot(data=counts, na.rm=T) + 
    geom_bar(aes(x=CDR, y=m, fill=sample_source), position="dodge", stat='identity') +
    scale_y_continuous(label=percent) +
    xlab("hbonds") +
    ggtitle("Cross Interface Hbonds")
  plot_field(p, "hbond_tog_hist_by_cdr")
  
})) # end FeaturesAnalysis