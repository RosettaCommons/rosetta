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
id = "ab_hbonds-cdr_cdr_den",
author = "Jared Adolf-Bryfogle",
brief_description = "CDR - CDR Hbonds",
feature_reporter_dependencies = c("AntibodyFeatures", "HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
  
  sele_don = "
SELECT
    hb.energy as energy,
    don.struct_id as struct_id,
    don.resNum as resnum,
    hb_geom.AHdist as distance,
    don_c.CDR as CDR1,
    acc_c.CDR as CDR2,
    don.atmType as don_atm,
    acc.atmType as acc_atm,
    don.HBChemType as don_type,
    acc.HBChemType as acc_type,
    don_ss.dssp as don_dssp,
    acc_ss.dssp as acc_dssp
  FROM
    hbond_sites AS don,
    hbond_sites AS acc,
    hbonds AS hb,
    hbond_geom_coords as hb_geom,
    cdr_residues as  don_c,
    cdr_residues as  acc_c,
    residue_secondary_structure as don_ss,
    residue_secondary_structure as acc_ss
  WHERE
    acc.struct_id == don.struct_id AND
    don.struct_id == hb.struct_id AND
    hb_geom.struct_id = hb.struct_id AND
    don_c.struct_id = hb.struct_id AND
    acc_c.struct_id = hb.struct_id AND
    don_ss.struct_id = hb.struct_id AND
    acc_ss.struct_id = hb.struct_id AND
    don.site_id =  hb.don_id AND
    acc.site_id = hb.acc_id  AND
    hb_geom.hbond_id = hb.hbond_id AND
    don_c.resNum = don.resNum AND
    acc_c.resNum = acc.resNum AND
    acc_ss.resNum = acc.resNum AND
    don_ss.resNum = don.resNum AND
    don_c.CDR != acc_c.CDR AND
    NOT (don_ss.dssp = acc_ss.dssp AND don_ss.dssp = 'E' AND acc_ss.dssp='E' AND 
    don.HBChemType = 'hbdon_PBA' and acc.HBChemType='hbacc_PBA')
  "

  sele_acc = "
SELECT
    hb.energy as energy,
    don.struct_id as struct_id,
    don.resNum as resnum,
    hb_geom.AHdist as distance,
    acc_c.CDR as CDR1,
    don_c.CDR as CDR2,
    don.atmType as don_atm,
    acc.atmType as acc_atm,
    don.HBChemType as don_type,
    acc.HBChemType as acc_type,
    don_ss.dssp as don_dssp,
    acc_ss.dssp as acc_dssp
  FROM
    hbond_sites AS don,
    hbond_sites AS acc,
    hbonds AS hb,
    hbond_geom_coords as hb_geom,
    cdr_residues as  don_c,
    cdr_residues as  acc_c,
    residue_secondary_structure as don_ss,
    residue_secondary_structure as acc_ss
  WHERE
    acc.struct_id == don.struct_id AND
    don.struct_id == hb.struct_id AND
    hb.struct_id == hb_geom.struct_id AND
    don_c.struct_id = hb.struct_id AND
    acc_c.struct_id = hb.struct_id AND
    don_ss.struct_id = hb.struct_id AND
    acc_ss.struct_id = hb.struct_id AND
    hb.don_id == don.site_id AND
    hb.acc_id == acc.site_id AND
    hb.hbond_id == hb_geom.hbond_id AND
    don_c.resNum = don.resNum AND
    acc_c.resNum = acc.resNum AND
    don.resNum = don_ss.resNum AND
    acc.resNum = acc_ss.resNum AND
    acc_ss.resNum = acc.resNum AND
    don_ss.resNum = don.resNum AND
    don_c.CDR != acc_c.CDR AND
    NOT (don_ss.dssp = acc_ss.dssp AND don_ss.dssp = 'E' AND acc_ss.dssp='E' AND 
    don.HBChemType = 'hbdon_PBA' and acc.HBChemType='hbacc_PBA')
  "
  
  sele_total_cdrs = "
  SELECT 
    struct_id,
    CDR
  FROM
    cdr_residues
  "
  #NOT (don.HBChemType == 'hbdon_PBA' AND acc.HBChemType == 'hbacc_PBA')
  
  don_data = query_sample_sources(sample_sources, sele_don, char_as_factor=F)
  acc_data = query_sample_sources(sample_sources, sele_acc, char_as_factor=F)
  total_cdrs = query_sample_sources(sample_sources, sele_total_cdrs, char_as_factor=F)
  
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
  don_data$pair = paste(don_data$CDR1, don_data$CDR2, sep="_")
  acc_data$type = "acc"
  acc_data$pair = paste(acc_data$CDR1, acc_data$CDR2, sep="_")
  data = rbind(don_data, acc_data)
  
  cdrs = c("L1", "L2", "L3", "H1", "H2", "H3")
  types = c("don", "acc")
  
  #Hbond Counts per CDR

  
  don_counts = data.frame()
  acc_counts = data.frame()
  
  print("Calculating hbond counts")
  counts <- ddply(data, .(sample_source, struct_id, type), function(data){
    counts = data.frame()
    for (outer_cdr in cdrs){
      for (inner_cdr in cdrs){
        if(inner_cdr == outer_cdr){next}
        
        
        n = nrow(data[data$energy <0 & data$CDR1==outer_cdr & data$CDR2 == inner_cdr,])
        c = data.frame(sample_source = as.character(data$sample_source[1]), struct_id = data$struct_id[1], n=n, CDR2= inner_cdr, CDR1 = outer_cdr, type=data$type[1])
        counts = rbind(counts, c)
      }
    }
    counts
  })
  
  print(head(counts))
  counts$nc = as.character(counts$n)
  field = "n"
  group = c("sample_source", "CDR1", "CDR2")
  
  #In order to get spacing correct, we need to use identity:
  perc <- ddply(counts, .(sample_source, CDR1, CDR2, type), function(data){
    perc <- ddply(data, .(n), function(d2){
      perc = nrow(d2)/nrow(data)
      data.frame(perc = perc)
    })
  })
  perc$nc = as.character(perc$n)
  
  perc_all <- ddply(counts, .(sample_source, CDR1, CDR2), function(data){
    perc <- ddply(data, .(n), function(d2){
      perc = nrow(d2)/nrow(data)
      data.frame(perc = perc)
    })
  })
  perc_all$nc = as.character(perc_all$n)
  
  print(head(perc))
  for (outer_cdr in cdrs){
      
      p <- ggplot(data=perc_all[perc_all$CDR1==outer_cdr,], na.rm=T) + 
      geom_bar(aes(x=nc, y= perc , fill=sample_source), position="dodge", stat='identity') +
      xlab("hbonds") +
      ylab("% of Sample Source") +
      scale_y_continuous(label=percent) +
      ggtitle(paste("Cross CDR H-Bond Counts", outer_cdr)) 
    plot_field(p, paste("hbond_counts", "hist_by_cdr", outer_cdr,"tog", sep="_"), ~ CDR2)
    
      p <- ggplot(data=perc_all[perc_all$CDR1==outer_cdr,], na.rm=T) + 
      geom_bar(aes(x=nc, y = perc, fill=sample_source), position="dodge", stat = 'identity') +
      xlab("hbonds") +
      ylab("% of Sample Source") +
      scale_y_continuous(label=percent) +
      ggtitle(paste("Cross CDR H-Bond Counts", outer_cdr))
    plot_field(p, paste("hbond_counts", "hist_by_all", outer_cdr,"tog", sep="_"))
    
    for (type in types){
      p <- ggplot(data=perc[perc$CDR1==outer_cdr & counts$type==type,], na.rm=T) + 
        geom_bar(aes(x=nc, y= perc, fill=sample_source), position="dodge", stat = 'identity') +
        xlab("hbonds") +
        ylab("% of Sample Source") +
        scale_y_continuous(label=percent) +
        ggtitle(paste("Cross CDR H-Bond Counts,", outer_cdr, "as", type))
      plot_field(p, paste("hbond_counts", "hist_by_cdr", outer_cdr,type, sep="_"), ~ CDR2)
    
        p <- ggplot(data=perc[perc$CDR1==outer_cdr & counts$type == type,], na.rm=T) + 
        geom_bar(aes(x=nc, y = perc, fill=sample_source), position="dodge", stat = 'identity') +
        xlab("hbonds") +
        ylab("% of Sample Source") +
        scale_y_continuous(label=percent) +
        ggtitle(paste("Cross CDR H-Bond Counts,", outer_cdr, "as", type))
      plot_field(p, paste("hbond_counts", "hist_by_all", outer_cdr,type, sep="_"))
    }

  }
  
    #Total Averages for each cdr
  
  avgs <- ddply(counts, .(sample_source, CDR1, CDR2, type), function(data){
    data.frame(sample_source=data$sample_source, CDR1=data$CDR1, CDR2=data$CDR2, m=mean(data$n), type=data$type)
  })
  
  #Histograms
  for (outer_cdr in cdrs){
    
    p <- ggplot(data=avgs[avgs$CDR1==outer_cdr & avgs$CDR2!=outer_cdr,], na.rm=T) + 
      geom_bar(aes(x=CDR2, y = m, fill=sample_source),position="dodge", stat='identity') +
      xlab("CDR2") +
      ylab("Avg hbonds") +
      ggtitle(paste("Average Cross CDR H-Bond Counts", outer_cdr))
    plot_field(p, paste("avg_hbond_counts", "hist_by_cdr_as_tog", outer_cdr, sep="_"))

    for (type in types){
      p <- ggplot(data=avgs[avgs$CDR1==outer_cdr & avgs$CDR2!=outer_cdr & avgs$type==type,], na.rm=T) + 
        geom_bar(aes(x=CDR2, y = m , fill=sample_source), position="dodge", stat='identity') +
        xlab("CDR2") +
        ylab("Avg hbonds") +
        ggtitle(paste("Average Cross CDR H-Bond Counts,",  outer_cdr,  "as", type))
      plot_field(p, paste("avg_hbond_counts", "hist_by_cdr", outer_cdr, type, sep="_"))
    }  
  }

  
  #Combined everages - CDR and everything else

  avgs <- ddply(counts, .(sample_source, CDR1, type), function(data){
    data.frame(sample_source=data$sample_source, CDR1=data$CDR1, CDR2=data$CDR2, m=mean(data$n), type = data$type)
  })
  
    p <- ggplot(data=avgs, na.rm=T) + 
      geom_bar(aes(x=CDR1, y=m, fill=sample_source), position="dodge", stat='identity') +
      xlab("CDR2") +
      ylab("Avg hbonds") +
      ggtitle(paste("Average Cross CDR H-Bond Counts"))
    plot_field(p, paste("avg_hbond_counts", "hist_all_tog", sep="_"))

  for (type in types){
    p <- ggplot(data=avgs[avgs$type==type,], na.rm=T) + 
      geom_bar(aes(x=CDR1, y = m , fill=sample_source), position="dodge", stat='identity') +
      xlab("CDR2") +
      ylab("Avg hbonds") +
      ggtitle(paste("Average Cross CDR H-Bond Counts as", type))
    plot_field(p, paste("avg_hbond_counts","hist_all",type, sep="_"))
  }
  
  
  #Hbond Energies
  field = "energy"
  group = c("sample_source", "CDR1", "CDR2")
  group2 = c("sample_source", "CDR1")

  for (outer_cdr in cdrs){

    
    dens <- estimate_density_1d(data[data$CDR1==outer_cdr,],  group, field)
    if(nrow(dens)>=1){
    p <- ggplot(data=dens, na.rm=T) + plot_parts +
      xlab("REU") +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(paste(outer_cdr, "Cross CDR H-Bond Energies", outer_cdr))
    plot_field(p, paste("hbond", field, "den_by_cdr", outer_cdr, "tog", sep="_"), ~ CDR2)
    }
    
    if(nrow(dens)>=1){
    dens <- estimate_density_1d(data[data$CDR1==outer_cdr,],  group2, field)
    p <- ggplot(data=dens, na.rm=T) + plot_parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(paste(outer_cdr, "Cross CDR H-Bond Energies", outer_cdr))
    plot_field(p, paste("hbond", field, "den_by_all", outer_cdr, "tog", sep="_"))
    }
    
    for (type in types){
    dens <- estimate_density_1d(data[data$CDR1==outer_cdr & data$type == type,],  group, field)
    if(nrow(dens) >=1){
      p <- ggplot(data=dens, na.rm=T) + plot_parts +
        xlab("REU") +
        geom_line(aes(x, y, colour=sample_source), size=1.2) +
        ggtitle(paste(outer_cdr,"Cross CDR H-Bond Energies,", outer_cdr, "as", type))
      plot_field(p, paste("hbond", field, "den_by_cdr", outer_cdr, type, sep="_"), ~ CDR2)
      }
    
      dens <- estimate_density_1d(data[data$CDR1==outer_cdr & data$type==type,],  group2, field)
      if(nrow(dens) >=1){
      p <- ggplot(data=dens, na.rm=T) + plot_parts +
        xlab("REU") +
        geom_line(aes(x, y, colour=sample_source), size=1.2) +
        ggtitle(paste(outer_cdr,"Cross CDR H-Bond Energies,", outer_cdr, "as", type))
      plot_field(p, paste("hbond", field,  "den_by_all", outer_cdr, type ,sep="_") )
      }
    }
    
    

  }
  
  #Hbond Distances
  field = "distance"
  group = c("sample_source", "CDR1", "CDR2")
  group1 = c("sample_source", "CDR1")
  for (outer_cdr in cdrs){
    
    dens <- estimate_density_1d(data[data$CDR1==outer_cdr,],  group, field)
    if(nrow(dens)>=1){
    p <- ggplot(data=dens, na.rm=T) + plot_parts +
      xlab("Angstroms") +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(paste(outer_cdr, "Cross CDR H-Bond Distances", outer_cdr))
    plot_field(p, paste("hbond", field, "den_by_cdr", outer_cdr, "tog",  sep="_"), ~ CDR2)
    }
    
    dens <- estimate_density_1d(data[data$CDR1==outer_cdr,],  group2, field)
    if(nrow(dens)>=1){
    p <- ggplot(data=dens, na.rm=T) + plot_parts +
      xlab("Angstroms") +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(paste(outer_cdr, "Cross CDR H-Bond Distances", outer_cdr))
    plot_field(p, paste("hbond", field, "den_by_all", outer_cdr, "tog",  sep="_"))
    }
    
    for (type in types){
      dens <- estimate_density_1d(data[data$CDR1==outer_cdr & data$type == type,],  group, field)
      if(nrow(dens)>=1){
      p <- ggplot(data=dens, na.rm=T) + plot_parts +
        xlab("Angstroms") +
        geom_line(aes(x, y, colour=sample_source), size=1.2) +
        ggtitle(paste(outer_cdr,"Cross CDR H-Bond Distances,", outer_cdr,  "as", type))
      plot_field(p, paste("hbond", field, "den_by_cdr", outer_cdr, type, sep="_"), ~ CDR2)
      }
      if(nrow(dens)>=1){
      dens <- estimate_density_1d(data[data$CDR1==outer_cdr & data$type == type,],  group2, field)
    
      p <- ggplot(data=dens, na.rm=T) + plot_parts +
        xlab("Angstroms") +
        geom_line(aes(x, y, colour=sample_source), size=1.2) +
        ggtitle(paste(outer_cdr,"Cross CDR H-Bond Distances,", outer_cdr, "as", type))
      plot_field(p, paste("hbond", field, "den_by_all", outer_cdr, type, sep="_"))
      }
    }
    
    

    
  }
  
})) # end FeaturesAnalysis