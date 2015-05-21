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
id = "ag_ab_contact_den",
author = "Jared Adolf-Bryfogle",
brief_description = "VL VH packing angle metrics",
feature_reporter_dependencies = c("AntibodyFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
  
  #First we run on all the interfaces in the database

  

  sele = "
  SELECT
    CDR,
    struct_id,
    length,
    ag_ab_contacts_total,
    ag_ab_contacts_nres
  FROM
    cdr_metrics
  "
  
  data = query_sample_sources(sample_sources, sele, char_as_factor=F)
 
  parts <- list(
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
  
  #Contacts Antigen
  binary_data <- ddply(data, .(sample_source, CDR), function(data){
    percent = (nrow(data[data$ag_ab_contacts_total >= 1,]))/length(data$ag_ab_contacts_total)
    data.frame(percent=percent)
  })
  p <- ggplot(data=binary_data, na.rm=T, aes(x=CDR)) +
    geom_bar(position="dodge", aes(y=percent, fill=sample_source), stat='identity') +
    ggtitle("CDR Makes Antigen Contact") +
    scale_y_continuous(label=percent) +
    ylab("% of sample source")
  plot_field(p, "cdr_makes_contact_hist")
  save_tables(self, binary_data, "cdr_makes_contact_table", sample_sources, output_dir, output_formats,
    caption="CDR makes contact", caption.placement="top", quote_strings=F)
      
  #Avg Contacts per Residue per CDR
  #Testing to make sure this is working:
  data$avg = data$ag_ab_contacts_total/data$length
  means = ddply(data, .(sample_source, CDR), function(data){
    data.frame(m= mean(data$avg))
  })
  p <- ggplot(data=means, na.rm=T, aes(x=CDR)) +
    geom_bar(position="dodge", aes(y=m, fill=sample_source), stat='identity') +
    ggtitle("Average Contacts per Residue") +
    xlab("CDR") +
    ylab("n")
  plot_field(p, "avg_contacts_per_residue_per_cdr_hist")
  save_tables(self, means, "avg_contacts_per_residue_per_cdr_table", sample_sources, output_dir, output_formats,
    caption="Avg contacts per residue per cdr", caption.placement="top", quote_strings=F)

  #Residues in contact per CDR
  means = ddply(data, .(sample_source, CDR), function(data){
    data.frame(m= mean(data$ag_ab_contacts_nres))
  })
  p <- ggplot(data=means, na.rm = T, aes(x=CDR)) +
    geom_bar(position="dodge", aes(y=m, fill=sample_source), stat='identity') +
    ggtitle("Average Residues in contact") +
    xlab("CDR") +
    ylab("n")
  plot_field(p, "avg_residues_in_contact_per_cdr_hist")
  save_tables(self, means, "avg_residues_in_contact_per_cdr_table", sample_sources, output_dir, output_formats,
    caption="Avg residues in contact per cdr", caption.placement="top", quote_strings=F)
  
  #This is for Brain -  average percentage of contacts that come from a CDR.  If the antibody is not in contact with antigen, we skip it.
  data2 = ddply(data, .(sample_source, struct_id), function(d2){
    total_contacts = sum(d2$ag_ab_contacts_total)
    data.frame(total_contacts = total_contacts, CDR = d2$CDR, ag_ab_contacts_total = d2$ag_ab_contacts_total)
  })
  
  avgs = ddply(data2[data2$total_contacts > 0,], .(sample_source, struct_id, CDR), function(d2){
    contacts = d2$ag_ab_contacts_total[1]/d2$total_contacts
    print(paste(contacts, d2$total_contacts))
    perc = contacts
    data.frame(perc = perc)
  })
  print(head(avgs))
  
  avg_perc = ddply(avgs, .(sample_source, CDR), function(d2){
    
    data.frame(m_perc = mean(d2$perc))
  })
  print(head(avg_perc))
  
  p <- ggplot(data=avg_perc, na.rm = T, aes(x=CDR)) +
    geom_bar(position="dodge", aes(y=m_perc, fill=sample_source), stat='identity') +
    ggtitle("Average Percent of total contacts") +
    xlab("CDR") +
    ylab("Avg %") +
    scale_y_continuous(label=percent)
  plot_field(p, "avg_perc_total_contacts_hist")
  save_tables(self, avg_perc, "avg_perc_total_contacts_table", sample_sources, output_dir, output_formats,
    caption="Avg Percent of total contacts", caption.placement="top", quote_strings=F)
  
})) # end FeaturesAnalysis