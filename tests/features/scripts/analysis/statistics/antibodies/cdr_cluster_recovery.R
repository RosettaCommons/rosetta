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
id = "cdr_cluster_recovery",
author = "Jared Adolf-Bryfogle",
brief_description = "Used for length and cluster recovery of CDRs.  Mainly for AbDesign program based on North Clusters,",
long_description = "First sample source should be Natives.  This is the reference.  Other sample sources are collections of decoys from different experiments.
Decoys should have the native's name in input_tag.",
feature_reporter_dependencies = c("CDRClusterFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


  len_sele <- "
  SELECT
    new.cdr_clusters.CDR as CDR,
    new.cdr_clusters.length as length,
    new.cdr_clusters.normDis_deg as normDis_deg
  FROM
    ref.cdr_clusters, new.cdr_clusters,
    ref.structures, new.structures
  WHERE
    new.cdr_clusters.struct_id=new.structures.struct_id AND new.structures.input_tag LIKE :like_tag AND
    ref.cdr_clusters.struct_id=ref.structures.struct_id AND ref.structures.input_tag = :tag AND
    new.cdr_clusters.CDR = ref.cdr_clusters.CDR AND
    new.cdr_clusters.length = ref.cdr_clusters.length;
  "

  clu_sele <- "
  SELECT
    new.cdr_clusters.CDR as CDR,
    new.cdr_clusters.fullcluster as cluster,
    new.cdr_clusters.normDis_deg as normDis_deg
  FROM
    ref.cdr_clusters, new.cdr_clusters,
    ref.structures, new.structures
  WHERE
    new.cdr_clusters.struct_id=new.structures.struct_id AND new.structures.input_tag LIKE :like_tag AND
    ref.cdr_clusters.struct_id=ref.structures.struct_id AND ref.structures.input_tag = :tag AND
    new.cdr_clusters.CDR = ref.cdr_clusters.CDR AND
    new.cdr_clusters.fullcluster = ref.cdr_clusters.fullcluster;
  "
  
  total_sele <- "
  SELECT
    cdr_clusters.fullcluster 
  FROM
    structures, cdr_clusters 
  WHERE
    structures.input_tag LIKE :like_tag AND structures.struct_id = cdr_clusters.struct_id AND
    cdr_clusters.CDR= :c;
  "
  
  #Need individual and combined data of CDR/PDB/sample_source. Next time, we do this the slower usual way.
  
  cdrs = c("L1", "L2", "L3", "H1", "H2", "H3")
  
  #This is a fixer function as I really should not have used ddply for everything. Inserts zeros into DF where needed.
  create_zero_data = function(current_data){
    sele = "SELECT DISTINCT cdr_clusters.CDR FROM cdr_clusters"
    #print(sample_sources)
    for (i in 2:nrow(sample_sources)){
      #ss_id = as.character(ss["sample_source"])
      ss = sample_sources[i,]
      summary(ss)
      ss_id = as.character(ss$sample_source)
      result = query_sample_source(ss, sele, char_as_factor=F)
      for (cdr in result$CDR){
        if (! any(current_data$CDR==cdr & current_data$sample_source == ss_id)){
          new_df = data.frame(CDR=as.character(cdr), sample_source = ss_id, normDis_deg=0, fullcluster="NA", stringsAsFactors=F)
          current_data = merge(current_data, new_df, all=T)
        }
      }
    }
    return(current_data)
  }
  #Get input_tags to match natives
  result = query_sample_source(sample_sources[1,], "SELECT input_tag from structures", char_as_factor=F)
  native_tags = result$input_tag
  get_and_write_recovery = function(sele, type){
    all_data = adply(native_tags, 1, function(native_tag) {
      
      native_tag_sp = unlist(strsplit(native_tag, "/"))
      pdb_sp = strsplit(native_tag_sp[length(native_tag_sp)], "\\.")
      tag = unlist(pdb_sp)[1] #2J88.pdb -> take 2J88
      
      cat("\nWorking on", tag, type, "recovery", "\n", sep=" ")
      
      match = paste("%", tag, "%", sep="")
      tag_frame = data.frame(like_tag=match, tag=native_tag)
    
      data = query_sample_sources_against_ref(sample_sources, sele, sele_args_frame=tag_frame, char_as_factor=F) 
      data = create_zero_data(data)
      #Type is length or cluster here:
      combine_data <- function(data, type){
        res_by_ss = ddply(data, "sample_source", function(data_by_ss, type){
          
          res_by_cdr = ddply(data_by_ss, "CDR", function(data_by_cdr, type) {
            
            cdr = data_by_cdr$CDR[1]
            total_data = query_sample_sources(sample_sources, total_sele, sele_args_frame=data.frame(like_tag=match, c=cdr), char_as_factor=F)
            total_decoys = length(total_data$fullcluster[total_data$sample_source == as.character(data_by_ss$sample_source[1])])
            
            if (length(data_by_cdr$normDis_deg[data_by_cdr$normDis_deg > 0]) == 0){
              recovery = 0
              recovery_total = 0
              angle_mean = 0
              angle_sd = 0
            }
            else{
              recovery = length(data_by_cdr$CDR)/total_decoys
              recovery_total = length(data_by_cdr$CDR)
              angle_mean = mean(data_by_cdr$normDis_deg[data_by_cdr$normDis_deg > 0])
              angle_sd = sd(data_by_cdr$normDis_deg[data_by_cdr$normDis_deg >0])
            }
            #cat("rec: ", recovery, "rec_total: ", recovery_total, "\n")
          
            result = data.frame(native=tag, sample_source=as.character(data_by_cdr$sample_source[1]), recovery = recovery, 
              total_rec=recovery_total, total=total_decoys, angle_mean=angle_mean, angle_sd = angle_sd)
            return(result)
          })
        })
      } #End combine_data
    
    rec_by_ss = combine_data(data, type)
    return(rec_by_ss)
    })
  
    
    #Recovery for CDR by individual PDB
    #print(all_data)
    all_data$X1 = NULL #Remove extra column
    grouped_data = sort(all_data, by = ~ native + CDR + sample_source)
    
    save_tables(self,
              grouped_data, paste("cdr_recovery_by_pdb", type, sep="_"), sample_sources, output_dir, output_formats,
              caption=paste("CDR", type, "recovery", sep=" "), caption.placement="top", quote_strings=F)
  
    grouped_data = sort(all_data, by = ~recovery + total_rec + CDR + sample_source)
  
    save_tables(self,
              grouped_data, paste("cdr_recovery_by_pdb_best",type, sep="_"), sample_sources, output_dir, output_formats,
              caption=paste("CDR", type, "recovery", sep=" "), caption.placement="top", quote_strings=F)
  
    calc_df = function(data){
      recovery = sum(data$total_rec)/sum(data$total)
      recovery_total = sum(data$total_rec)
      total = sum(data$total)
    
      #Note not totally correct mean and SD here. (mean of mean)
      angle_mean = mean(data$angle_mean[data$angle_mean>0])
      angle_sd = sd(data$angle_mean[data$angle_mean>0])
      df = data.frame(recovery = recovery, 
          total_rec=recovery_total, total=total, angle_mean=angle_mean, angle_sd = angle_sd)
    }
  
    
    #Recovery by CDR by source
    grouped_data = ddply(all_data, "sample_source", function(ss_data){
        f = ddply(ss_data, "CDR", function(cdr_data){
          df = calc_df(cdr_data)
        })
    })
    grouped_data = sort(grouped_data, by = ~ sample_source + CDR)
  
    save_tables(self,
              grouped_data, paste("cdr_recovery_by_sample_source", type, sep="_"), sample_sources, output_dir, output_formats,
              caption=paste("CDR", type, "recovery", sep=" "), caption.placement="top", quot_strings=F)
    
    grouped_data = sort(grouped_data, by= ~ recovery +total_rec)
  
    save_tables(self,
              grouped_data, paste("cdr_recovery_by_sample_source_best",type, sep="_"), sample_sources, output_dir, output_formats,
              caption=paste("CDR", type, "recovery", sep=" "), caption.placement="top", quote_strings=F)
  
  
    
    #Recovery by CDR:
    grouped_data = ddply(all_data, "CDR", function(cdr_data){
      df = calc_df(cdr_data)
    })
 
    grouped_data = sort(grouped_data, by= ~ CDR + recovery)
  
    save_tables(self,
              grouped_data, paste("overall_cdr_recovery",type, sep="_"), sample_sources, output_dir, output_formats,
              caption=paste("CDR", type, "recovery", sep=" "), caption.placement="top", quote_strings=F)
  
    
    #Reovery by Native:
    grouped_data = ddply(all_data, "native", function(cdr_data) {
      df = calc_df(cdr_data)
    })
  
    grouped_data = sort(grouped_data, by=~ native + recovery)
    save_tables(self,
              grouped_data, paste("overall_native_recovery",type, sep="_"), sample_sources, output_dir, output_formats,
              caption=paste("CDR", type, "recovery", sep=" "), caption.placement="top", quote_strings=F)
  
    
    #Recovery by Source:
    grouped_data = ddply(all_data, "sample_source", function(cdr_data) {
      df = calc_df(cdr_data)
    })
    grouped_data = sort(grouped_data, by= ~ sample_source + recovery)
    save_tables(self,
             grouped_data, paste("overall_recovery_by_sample_source",type, sep="_"), sample_sources, output_dir, output_formats,
             caption=paste("CDR", type, "recovery", sep=" "), caption.placement="top", quote_strings=F)
  }

  get_and_write_recovery(len_sele, "length")
  get_and_write_recovery(clu_sele, "cluster")

 #DIR: cdr_cluster_recovery
    #5 Tables?:
     #cdr_by_pdb, cdr_by_source, by_cdr, by_pdb, by_source  
  
    
})) # end FeaturesAnalysis
  
  