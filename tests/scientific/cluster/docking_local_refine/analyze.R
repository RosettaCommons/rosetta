#! /usr/bin/Rscript --vanilla
# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


library(plyr)
library(ggplot2)
library(rjson)

# Extract the score data from it a rosetta generated .sc or silent
# file
read_score_file_data <- function(score_filename){
	cat("extracting score data from ", score_filename, " ...\n", sep = "")
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

get_all_score_files_in_directory <- function(path, extension){
	z <- list.files(path=path,	pattern=paste("*", extension, sep=""))
	if(length(z) == 0){
		cat("WARNING: No score files were found in the '", path, "', with extension '", extension, "'.", sep="")
		cat("That path does contain these files, however,\n",
			paste(dir(path), collapse="\n\t"), sep="")
	}
	z
}

extract_data_from_score_files <- function(filenames){
	df <- ldply(score_filenames, function(score_filename){
	  z <- read_score_file_data(score_filename)
		if(nrow(z) < 2L){
			cat("WARNING: No score rows were extracted from '", score_filename, "'.", sep="")
		}
		z
	})
	df$target <- factor(
		substr(
			df$description, 0L,
			nchar(as.character(df$description[1L])) - 5L))
	df
}


generate_score_vs_rmsd_plots <- function(
	data,
	output_dir) {
	plot_id <-"score_vs_rms.png"
	if( !("rms" %in% names(data))){
		cat("WARNING: Attepting to generate score vs rmsd plot but the data does not contain the column 'rms'", sep="")
		print(summary(data))
		stop()
	}
	if( !("score" %in% names(data))){
		cat("WARNING: Attepting to generate score vs rmsd plot but the data does not contain the column 'score'", sep="")
		print(summary(data))
		stop()
	}

	# this makes a lattice of plots by target (assumes ~50 targets)
	p <- ggplot(data) + theme_bw() +
		geom_point(aes(x=rms, y=score), size=.7) +
		facet_wrap( ~ target, scales="free_y") +
		opts(title="Score vs RMS for Docking Local Refine Targets") +
		scale_y_continuous("Score") +
		scale_x_continuous("RMSD", limits=c(0,5.5))
	ggsave(file.path(output_dir, plot_id), p, width=19, height=15)
}


generate_interface_score_vs_interface_rmsd_plots <- function(
	data,
	output_dir
){
	plot_id <-"interface_score_vs_interface_rms.png"

	if( !("Irms" %in% names(data))){
		cat("WARNING: Attepting to generate interface score vs interface rmsd plot but the data does not contain the column 'Irms'", sep="")
		print(summary(data))
		stop()
	}
	if( !("score" %in% names(data))){
		cat("WARNING: Attepting to generate interface score vs interface rmsd plot but the data does not contain the column 'I_sc'", sep="")
		print(summary(data))
		stop()
	}

	# this makes a lattice of plots by target (assumes ~50 targets)
	p <- ggplot(data) + theme_bw() +
		geom_point(aes(x=Irms, y=I_sc), size=.7) +
		facet_wrap( ~ target, scales="free_y") +
		opts(title="Interface Score vs Interface RMS for Docking Local Refine Targets") +
		scale_y_continuous("Interface Score") +
		scale_x_continuous("Interface RMSD", limits=c(0,.8))
	ggsave(file.path(output_dir, plot_id), p, width=15, height=15)
}

generate_sample_depth_plots <- function(
	data,
	output_dir
){

	sample_sizes <-
		floor(
			seq(10, max( ddply(data, .(target), nrow)$V1), length.out=100))

	df <- ldply(sample_sizes, function(sample_size){
		ddply(data, .(target), function(per_target_df){
			z <- per_target_df[
				sample(nrow(per_target_df), min(nrow(df), sample_size), replace=T),]
			print(summary(z))
			data.frame(
				target=z$target[1],
				min_rms = min(z$rms),
				min_Irms = min(z$Irms),
				sample_size = sample_size)

		})
	})

	plot_id <- "sample_depth_rmsd.png"
	p <- ggplot(df) + theme_bw() +
		stat_smooth(aes(x=sample_size, y=min_rms, colour=target, group=target), se=F) +
		opts(title="Sample Depth: Min RMSD by Sample Size") +
		scale_colour_discrete(legend=F) +
		scale_x_continuous("Sample Size") +
		scale_y_continuous("Minimum RMSD")
	ggsave(file.path(output_dir, plot_id), p)

	plot_id <- "sample_depth_interface_rmsd.png"
	p <- ggplot(df) + theme_bw() +
		stat_smooth(aes(x=sample_size, y=min_Irms, colour=target, group=target), se=F) +
		opts(title="Sample Depth: Min Interface RMSD by Sample Size") +
		scale_colour_discrete( legend=F) +
		scale_x_continuous("Sample Size") +
		scale_y_continuous("Minimum Interface RMSD")
	ggsave(file.path(output_dir, plot_id), p)

}

# Generate a yaml file with summary statistics
#
# summary statistic:
# the minimum interface rms over the top 10 scoring structures
generate_summary_statistics <- function(
	data,
	output_dir,
	output_filename
){
	opt_df <- ddply(data, .(target), function(target_df){
		top_df <- head(
			target_df[order(target_df$I_sc),],
			n=10)
		data.frame(
			target=target_df$target[1],
			min_Irms=min(top_df$Irms))
	})
	mean_opt_df <- list(
		MeanMinimumInterfaceRMSDOverTop10ScoringInterfaces=
			mean(opt_df$min_Irms))
	cat(toJSON(mean_opt_df), "\n",
		sep="", file=file.path(output_dir, output_filename), append=F)

}



###########################################

config <- fromJSON(file="config.py")

score_filenames <- get_all_score_files_in_directory(
	config$output_score_path,
	config$output_scorefile_extension)

sc_data <- extract_data_from_score_files(score_filenames)

generate_score_vs_rmsd_plots(
	data=sc_data,
	output_dir=config$output_analysis_path)

generate_interface_score_vs_interface_rmsd_plots(
	data=sc_data,
	output_dir=config$output_analysis_path)

generate_sample_depth_plots(
	data=sc_data,
	output_dir=config$output_analysis_path)

generate_summary_statistics(
	data=sc_data,
	output_dir=".",
	output_filename=config$test_results_yaml)


