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
id = "ab_hbonds-cdr_ag_den",
author = "Jared Adolf-Bryfogle",
brief_description = "CDR - Antigen Hbonds.  Must have LH_A analyzed by features reporter for this to work",
feature_reporter_dependencies = c("AntibodyFeatures", "HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

    if ("FALSE" %in% opt$options$include_cdr4 & "FALSE" %in% opt$options$cdr4_only){
      sele_don = "
  SELECT 
      hb.energy as energy,
      hb.struct_id as struct_id,
      don.resNum as resnum,
      hb_geom.AHdist as distance,
      cdr_residues.CDR
  FROM
      interface_residues as acc_res,
      hbond_sites AS don,
      hbond_sites AS acc,
      hbonds AS hb,
      hbond_geom_coords as hb_geom,
      cdr_residues
  WHERE
      acc.struct_id == hb.struct_id AND
      don.struct_id == hb.struct_id AND
      hb_geom.struct_id = hb.struct_id AND
      cdr_residues.struct_id = hb.struct_id AND
      acc_res.struct_id = hb.struct_id AND
      acc.struct_id == hb.struct_id AND 
      don.site_id = hb.don_id AND
      acc.site_id = hb.acc_id  AND
      hb_geom.hbond_id = hb.hbond_id  AND
      cdr_residues.resNum = don.resNum AND
      acc_res.resNum = acc.resNum AND
      acc_res.side == 'side2' AND
      CDR NOT LIKE '%Proto%'
   "
      
      sele_acc = "
  SELECT
    hb.energy as energy,
    hb.struct_id as struct_id,
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
    acc.struct_id = hb.struct_id AND
    don.struct_id = hb.struct_id AND
    hb_geom.struct_id = hb.struct_id AND
    cdr_residues.struct_id = hb.struct_id AND
    don_res.struct_id = hb.struct_id AND
    don.site_id = hb.don_id AND
    acc.site_id = hb.acc_id AND
    hb_geom.hbond_id = hb.hbond_id AND
    cdr_residues.resNum = acc.resNum AND
    don_res.resNum = don.resNum  AND
    don_res.side = 'side2' AND
    CDR NOT LIKE '%Proto%'
  "
      
      sele_total_cdrs = "
  SELECT 
    struct_id,
    CDR
  FROM
    cdr_residues where CDR NOT LIKE '%Proto%'"
      
    }
    
    if ("TRUE" %in% opt$options$include_cdr4){
      sele_don = "
  SELECT 
      hb.energy as energy,
      hb.struct_id as struct_id,
      don.resNum as resnum,
      hb_geom.AHdist as distance,
      cdr_residues.CDR
      FROM
      interface_residues as acc_res,
      hbond_sites AS don,
      hbond_sites AS acc,
      hbonds AS hb,
      hbond_geom_coords as hb_geom,
      cdr_residues
      WHERE
      acc.struct_id == hb.struct_id AND
      don.struct_id == hb.struct_id AND
      hb_geom.struct_id = hb.struct_id AND
      cdr_residues.struct_id = hb.struct_id AND
      acc_res.struct_id = hb.struct_id AND
      acc.struct_id == hb.struct_id AND 
      don.site_id = hb.don_id AND
      acc.site_id = hb.acc_id  AND
      hb_geom.hbond_id = hb.hbond_id  AND
      cdr_residues.resNum = don.resNum AND
      acc_res.resNum = acc.resNum AND
      acc_res.side == 'side2'"
      
      sele_acc = "
      SELECT
      hb.energy as energy,
      hb.struct_id as struct_id,
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
      acc.struct_id = hb.struct_id AND
      don.struct_id = hb.struct_id AND
      hb_geom.struct_id = hb.struct_id AND
      cdr_residues.struct_id = hb.struct_id AND
      don_res.struct_id = hb.struct_id AND
      don.site_id = hb.don_id AND
      acc.site_id = hb.acc_id AND
      hb_geom.hbond_id = hb.hbond_id AND
      cdr_residues.resNum = acc.resNum AND
      don_res.resNum = don.resNum  AND
      don_res.side = 'side2'"
      
      sele_total_cdrs = "
      SELECT 
      struct_id,
      CDR
      FROM
      cdr_residues"
      
    }
    
    if ("TRUE" %in% opt$options$cdr4_only){
      sele_don = "
      SELECT 
      hb.energy as energy,
      hb.struct_id as struct_id,
      don.resNum as resnum,
      hb_geom.AHdist as distance,
      cdr_residues.CDR
      FROM
      interface_residues as acc_res,
      hbond_sites AS don,
      hbond_sites AS acc,
      hbonds AS hb,
      hbond_geom_coords as hb_geom,
      cdr_residues
      WHERE
      acc.struct_id == hb.struct_id AND
      don.struct_id == hb.struct_id AND
      hb_geom.struct_id = hb.struct_id AND
      cdr_residues.struct_id = hb.struct_id AND
      acc_res.struct_id = hb.struct_id AND
      acc.struct_id == hb.struct_id AND 
      don.site_id = hb.don_id AND
      acc.site_id = hb.acc_id  AND
      hb_geom.hbond_id = hb.hbond_id  AND
      cdr_residues.resNum = don.resNum AND
      acc_res.resNum = acc.resNum AND
      acc_res.side == 'side2' AND
      CDR LIKE '%Proto%'"
      
      sele_acc = "
      SELECT
      hb.energy as energy,
      hb.struct_id as struct_id,
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
      acc.struct_id = hb.struct_id AND
      don.struct_id = hb.struct_id AND
      hb_geom.struct_id = hb.struct_id AND
      cdr_residues.struct_id = hb.struct_id AND
      don_res.struct_id = hb.struct_id AND
      don.site_id = hb.don_id AND
      acc.site_id = hb.acc_id AND
      hb_geom.hbond_id = hb.hbond_id AND
      cdr_residues.resNum = acc.resNum AND
      don_res.resNum = don.resNum  AND
      don_res.side = 'side2' AND
      CDR LIKE '%Proto%'"
      
      sele_total_cdrs = "
      SELECT 
      struct_id,
      CDR
      FROM
      cdr_residues where CDR = 'Proto_H4' or CDR = 'Proto_L4'"
      
    }
  
  don_data = query_sample_sources(sample_sources, sele_don, char_as_factor=F)
  acc_data = query_sample_sources(sample_sources, sele_acc, char_as_factor=F)
  total_data = query_sample_sources(sample_sources, sele_total_cdrs, char_as_factor=F)
  
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
  
  #Hbond Counts per CDR
  print("Calculating hbond counts")
  counts <- ddply(total_data, .(sample_source, struct_id, CDR), function(totals){
    total_cdrs = nrow(totals)
    ndon = nrow(data[data$sample_source == totals$sample_source & data$struct_id == totals$struct_id & data$CDR == totals$CDR & data$type=="don",])
    nacc = nrow(data[data$sample_source == totals$sample_source & data$struct_id == totals$struct_id & data$CDR == totals$CDR & data$type=="acc",])
    
    if (is.null(ndon)){ndon = 0}
    if (is.null(nacc)){nacc = 0}
    
    don = data.frame(n=ndon, sample_source = as.character(totals$sample_source[1]), struct_id = totals$struct_id[1], CDR = totals$CDR[1], type="don")
    acc = data.frame(n=nacc, sample_source = as.character(totals$sample_source[1]), struct_id = totals$struct_id[1], CDR = totals$CDR[1], type="acc")
    
    counts = rbind(don, acc)
    counts
  })
  
  
  counts$nc = as.character(counts$n)
  
  types = c("don", "acc")
  
  #print(head(hbond_counts))
  field = "n"
  group = c("sample_source", "CDR")
  
  dens <- estimate_density_1d(counts,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("hbonds") +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Cross Interface hbonds")
  plot_field(p, "hbond_counts_den_by_cdr_tog", ~ CDR)
  
  for (type in types){
    dens <- estimate_density_1d(counts[counts$type == type,],  group, field)
    p <- ggplot(data=dens, na.rm=T) + plot_parts +
      xlab("hbonds") +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(paste("Cross Ab-Ag Hydrogen Bonds, CDR", type))
    plot_field(p, paste("hbond_counts_den_by_cdr", type, sep="_"), ~ CDR)

    dens <- estimate_density_1d(counts[counts$type == type,],  group, field)
    p <- ggplot(data=dens, na.rm=T) + plot_parts +
      xlab("hbonds") +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(paste("Cross Ab-Ag Hydrogen Bonds, CDR", type))
    plot_field(p, paste("hbond_counts_den_by_cdr", type, sep="_"), ~ CDR)
  }

    
  avgs <- ddply(counts, .(sample_source, CDR, type), function(d){
    data.frame( m=mean(d$n))
  })
  
  perc <- ddply(counts, .(sample_source, CDR, type), function(d){
    perc <- ddply(d, .(n), function(d2){
      perc = nrow(d2)/nrow(d) * 100
      data.frame(perc)
    })
  })
  
  perc_all <- ddply(counts, .(sample_source, CDR), function(d){
    perc <- ddply(d, .(n), function(d2){
      perc = nrow(d2)/nrow(d) * 100
      data.frame(perc)
    })
  })
  
  
  #Histograms by number of hbonds:
  for (type in types){

    p <- ggplot(data=perc[counts$perc == type,], na.rm=T) + 
      geom_bar(aes(x=as.character(n), y = perc, fill=sample_source), position="dodge", stat='identity') +
      xlab("hbonds") +
      ylab("% of Sample Source") +
      ggtitle(paste("Cross Interface H-bonds, CDR", type))
    plot_field(p, paste("hbond_counts_hist_by_cdr", type, sep="_"), ~ CDR)
  }

  p <- ggplot(data=perc_all, na.rm=T) + 
    geom_bar(aes(x=as.character(n), y=perc, fill=sample_source), position="dodge", stat='identity') +
    xlab("hbonds") +
    ylab("% of Sample Source") +
    ggtitle("Cross Interface H-bonds")
  plot_field(p, "hbond_counts_hist_by_cdr_tog", ~CDR)
  
  
  #Average Histograms
  for (type in types){

    p <- ggplot(data=avgs[avgs$type == type,], na.rm=T) + 
      geom_bar(aes(x=CDR, y = m, fill=sample_source), position="dodge", stat='identity') +
      xlab("CDR") +
      ylab("avg n") +
      ggtitle(paste("Average Cross Interface H-Bonds, CDR", type))
    plot_field(p, paste("avg_hbond_counts_hist_by_cdr", type, sep="_"))
  }

  p <- ggplot(data=avgs, na.rm=T) + 
    geom_bar(aes(x=CDR, y=m, fill=sample_source), position="dodge", stat='identity') +
    xlab("CDR") +
    ylab("avg n") +
    ggtitle("Average Cross Interface H-Bonds")
  plot_field(p, "avg_hbond_counts_hist_by_cdr_tog")
  
  p <- ggplot(data=avgs, na.rm=T) + 
    geom_bar(aes(x=type, y=m, fill=sample_source), position="dodge", stat='identity') +
    xlab("CDR") +
    ylab("avg n") +
    ggtitle("Average Cross Interface H-Bonds")
  plot_field(p, "avg_hbond_counts_hist_acc_vs_don_by_cdr", ~ CDR)
  
  
  #Hbond Energy density
  field = "energy"
  group = c("sample_source", "CDR")
  dens <- estimate_density_1d(don_data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("REU") +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Cross Ab-Ag Hydrogen Bond Energies, CDR Donor")
  plot_field(p, paste("hbond", "den_by_cdr", field, "don", sep="_"), ~ CDR)
  
  dens <- estimate_density_1d(acc_data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("REU") +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Cross Ab-Ag Hydrogen Bond Energies, CDR Acceptor")
  plot_field(p, paste("hbond", "den_by_cdr", field, "acc", sep="_"), ~ CDR)
  
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("REU") +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Cross Ab-Ag Hydrogen Bond Energies")
  plot_field(p, paste("hbond", "den_by_cdr", field, "tog", sep="_"), ~ CDR)
  
  #Hbond Distances
  field = "distance"
  group = c("sample_source", "CDR")
  dens <- estimate_density_1d(don_data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("Angstroms") +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Cross Ab-Ag Hydrogen Bond Distances, CDR Donor")
  plot_field(p, paste("hbond", "den_by_cdr", field, "don", sep="_"), ~ CDR)
  
  dens <- estimate_density_1d(acc_data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("Angstroms") +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Cross Ab-Ag Hydrogen Bond Distances, CDR Acceptor")
  plot_field(p, paste("hbond", "den_by_cdr", field, "acc", sep="_"), ~ CDR)
  
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("Angstroms") +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Cross Ab-Ag Hydrogen Bond Distances")
  plot_field(p, paste("hbond", "den_by_cdr",field, "tog",sep="_"), ~ CDR)
  
  
})) # end FeaturesAnalysis