#! /usr/bin/Rscript --vanilla
# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


library(plyr)
library(ggplot2)
library(rjson)

config <- list(
  fromJSON(file="_arguments.py")
  fromJSON(file="config.py"))

# Rosetta runs generate a .sc score file. This extracts the score data
# from it.
read_score_file_data <- function(score_filename){
  sc_file <- grep("^SCORE:", readLines(score_filename), value=T)
  header_row <- sc_file[1]
  data_rows <- grep(header_row, sc_file, value=T, invert=T)
  data <- read.table(
    con <- textConnection(data_rows),
    strip.white=T,
    col.names=strsplit(header_row, "\\s+", perl=T)[[1]])
  close(con)
  data
}



score_filenames <- list.files(
	config$output_score_path,
	paste("*", config$output_scorefile_extension, sep=""))
df <- ldply(score_files, function(score_filename){
  read_score_file_data(score_filename)
})
df$target <- substr(df$description, 0, 9)


plot_id <-"score_vs_rms.png" 
p <- ggplot(df) + theme_bw() +
 geom_point(aes(x=rms, y=score)) +
 facet_wrap( ~ target, scales="free_y") +
 opts(title="Score vs RMS for Docking Local Refine Targets")
ggsave(file.path(config$output_analysis_path, plot_id), p)

plot_id <-"interface_score_vs_interface_rms.png" 
p <- ggplot(df) + theme_bw() +
 geom_point(aes(x=Irms, y=I_sc)) +
 facet_wrap( ~ target, scales="free_y") +
 opts(title="Interface Score vs Interface RMS for Docking Local Refine Targets")
ggsave(file.path(config$output_analysis_path, plot_id), p)


# summary statistic:
# the minimum interface rms over the top n scoring structures
opt_df <- ddply(df, .(target), function(target_df){
	top_df <- head(
		target_df[order(target_df$I_sc),],
		n=config$num_top_structures)
	data.frame(
		target=target_df$target[1],
		min_Irms=min(top_df$Irms))
})

avg_opt_df <- list(
  AverageMinimumInterfaceRMSOverTopScoringInterfaces=mean(opt_df$min_Irms))

cat(toJSON(avg_opt_df), config$test_results_yaml)
