#!/usr/bin/env Rscript
# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

# The base_dir is usually rosetta/rosetta_tests/features
initialize_base_dir <- function(){
	command_args <- commandArgs(trailingOnly = FALSE)
	base_dir <- dirname(substring(command_args[grep("--file=", command_args)], 8))
	if(length(base_dir) == 0){
		base_dir = getwd()
	}
	base_dir
}

# The workdir is where the command was run from in particular the
# output dir can be specified relative to the work_dir
initialize_work_dir <- function(){
	getwd()
}

# Running ./compare_sample_sources.R will auto-generate
# compare_sample_sources_iscript.R which can be run in an interactive
# R session to re-run the last execution, like this:
#
#   source("compare_sample_sources_iscript.R")
#
initialize_iscript <- function(base_dir){
	source(paste(base_dir, "scripts/methods/generate_compare_sample_sources_iscript.R", sep="/"))
	iscript_header()
	iscript_base_dir(base_dir)
}

# Some packages are needed to setup and read the command line options
# in order to support interatively loading packages load most of the
# packages after the options have been initialized.
initialize_packages_before_options <- function(base_dir){
	# load_packages() will help the user to install the packages if they are missing
	source(paste(base_dir, "scripts/methods/load_packages.R", sep="/"))
	libraries <- c(
		"plyr",
		"optparse"
	)

	load_packages(
		libraries,
		fail_on_missing_package=T)
	iscript_libraries(libraries)
}

# Some method scripts are needed to setup and read the command line options
# in order to support interatively loading packages load most of the
# packages after the options have been initialized.
initialize_method_scripts_before_options <- function(base_dir) {
	includes <- c(
		"scripts/methods/output_formats.R"
	)
	for(inc in includes){
		tryCatch(source(paste(base_dir, inc, sep="/")), error=function(e){
			cat(paste(
				"ERROR: Failed to parse the methods script: '", inc, "'",
				"ERROR: It failed with the following error message:\n", e,
				sep=""))
			stop()
		})
	}
	iscript_includes(base_dir, includes)
}



initialize_command_line_options <- function() {
	option_list <- list()

	# Run level options
	option_list <- c(option_list, list(
		make_option(c("-v", "--verbose"), action="store_true", default=FALSE, dest="verbose",
			help="Print extra output [Default \"%default\"]"),
		make_option(c("--fail_on_missing_packages"), action="store_true", default=FALSE, dest="fail_on_missing_packages",
			help="Do not try to interactively install missing packages; simply fail. [Default \"%default\"]"),
		make_option(c("--dry_run"), action="store_true", type="logical", default=FALSE, dest="dry_run",
			help="Debug the analysis scripts but do not run them.  [Default \"%default\"]"),
		make_option(c("--db_cache_size"), action="store_true", type="integer", default=10000, dest="db_cache_size",
			help="Number of 1k pages of cache to use for database queries.  [Default \"%default\"]"),
		make_option(c("--add_footer"), action="store_true", type="logical", default=TRUE, dest="add_footer",
			help="Add footer to plots saying the analysis script and run date."),
		make_option(c("--generate_website"), action="store_true", type="logical", default=FALSE, dest="generate_website",
			help="Add footer to plots saying the analysis script and run date. [Default \"%default\"]"),
		make_option(c("--config"), action="store", type="character", default=NULL, dest="config_filename",
			help=". [Default \"%default\"]"),
		make_option(c("--ncores"), action="store_true", type="integer", default=1, dest="ncores",
			help="Run in parallel with specified number of cores.")))

	# Density estimation options
	option_list <- c(option_list, list(
		make_option(c("--adjust_kernel"), action="store_true", type="integer", default=1, dest="general_kernel_adjust",
			help="Multiplicative factor in the kernel bandwith for density estimation. [Default \"%default\"]")))

	# Which analysis scripts should be run
	option_list <- c(option_list, list(
		make_option(c("-i", "--script"), action="store", type="character", default=NULL, dest="script",
			help="Path to a single analysis script.  [Default \"%default\"]"),
		make_option(c("-a", "--analysis_dir"), type="character", default=NULL, dest="analysis_dir",
			help="Directory where the analysis scripts are located. The supplied directory is searched recursively for files of the form \"*.R\".  [Default \"%default\"]")))

	# Where should the results be stored?
	option_list <- c(option_list, list(
		make_option(c("-o", "--output_dir"), type="character", default="build", dest="output_dir",
			help="Directory where the output plots and statistics will be generated.  [Default \"%default\"]")))
  
  option_list <- c(option_list, list(
    make_option(c("--out_disable_ss_names"), action="store_true", default=FALSE, dest="out_disable_ss_names",
      help="Disable concatonating the ss names for output.  Be sure to combine this with out_dir or all plots will go to build. [Default \"%default\"]")))
  
  # Antibody Analysis Options
  option_list <- c(option_list, list(
  	make_option(c("--include_cdr4"), action="store_true", default=FALSE, dest="include_cdr4",
		  help="Include Proto CDR4 data in the analysis framework for specific antibody analysis scripts"),
		make_option(c("--cdr4_only"), action="store_true", default=FALSE, dest="cdr4_only",
		  help="make plots based only on Proto CDR4 data for specific antibody analysis scripts")))
  
	option_list <- c(option_list, make_output_formats_options_list(all_output_formats))

	opt <- parse_args(OptionParser(option_list=option_list), positional_arguments=TRUE)
}

check_basic_input <- function(opt){
	data_sources <- opt$args
	if(!opt$options$dry_run & length(data_sources) == 0 & is.null(opt$options$config_filename)){
		cat(
			"ERROR: No configuration file or sample_source databases was supplied",
			"",
			"##################################################",
			"NAME",
			"    compare_sample_source.R",
			"",
			"SYNOPSIS",
			"    ./compare_sample_sources.R [OPTIONS] --script <analysis_script> features_<ss_id1>.db3 [features_<ss_id2>.db3 ...]",
			"    ./compare_sample_sources.R [OPTIONS] --analysis_dir <analysis_dir> features_<ss_id1>.db3 [features_<ss_id2>.db3 ...]",
			"",
			"    From within an R session:",
			"       source(\"compare_sample_sources_iscript.R\") # This will re-run the last './compare_sample_sources.R' interactively",
			"",
			"DESCRIPTION",
			"    Compare structural features coming from different sample sources.",
			"",
			"    See ./compare_sample_sources.R --help for available options.",
			"    See https://wiki.rosettacommons.org/index.php/FeaturesScientificBenchmark for application documentation.",
			"    To cite: contact mattjomeara@gmail.com as the work is in progress.",
			"",
			"EXAMPLE",
			"    To compare how Rosetta distorts native structures, extract feature databases for a set of structures from the pdb",
			"    and the same set of structures optimized with your favorite prediction protocol. Assume the resulting feature databases",
			"    are stored in 'features_natives.db3' and 'features_rosetta.db3'. To compare the length of hydrogen bonds",
			"    conditional on the donor and acceptor chemical types between the two sample sources, run:",
			"",
			"    ./compare_sample_sources.R --script scripts/analysis/plots/hbonds/AHdist_chem_type.R features_rosetta.db3 features_rosetta.db3",
			"",
			"AUTHOR",
			"    Matthew O'Meara (mattjomeara@gmail.com)\n", sep="\n")
	  quit()
	}
}



initialize_packages <- function(opt, base_dir){
	# load_packages() will help the user to install the packages if they are missing
	source(paste(base_dir, "scripts/methods/load_packages.R", sep="/"))
	libraries <- c(
		"methods",
		"reshape",
		"plyr",
		"proto",
		"ggplot2",
		"RSQLite",
		"plotrix",
		"polynom",
		"rjson",
		"xtable",
    "scales")

#		"logspline", # APL TEMP!


	if(!is.null(opt$ncores)){
		libraries <- c(libraries, "doMC")
	}

	load_packages(
		libraries,
		fail_on_missing_package=opt$options$fail_on_missing_packages)
	iscript_libraries(libraries)
}

initialize_method_scripts <- function(base_dir) {
	includes <- c(
		"scripts/methods/sample_sources.R",
		"scripts/methods/methods.R",
		"scripts/methods/density_estimation.R",
		"scripts/methods/primary_modes.R",
		"scripts/methods/features_analysis.R",
		"scripts/methods/ggplot2_geom_indicator.R",
		"scripts/methods/ggplot2_scales.R",
		"scripts/methods/instancer.R",
		"scripts/methods/output_formats.R",
		"scripts/methods/save_plots.R",
		"scripts/methods/save_tables.R",
		"scripts/methods/coordinate_normalizations.R",
		"scripts/methods/comparison_statistics.R",
		"scripts/methods/generate_plot_webpage.R",
		"scripts/methods/vector_math.R",
		"scripts/methods/color_palettes.R",
    "scripts/methods/util.R")
	for(inc in includes){
		tryCatch(source(paste(base_dir, inc, sep="/")), error=function(e){
			cat(paste(
				"ERROR: Failed to parse the methods script: '", inc, "'",
				"ERROR: It failed with the following error message:\n", e,
				sep=""))
			stop()
		})
	}
	iscript_includes(base_dir, includes)
}

initialize_parallel_backend <- function(opt, base_dir){
	if(!is.null(opt$ncores)){
		registerDoMC(opt$ncores)
		use_parallel <<- TRUE
	} else {
		use_parallel <<- FALSE
	}
	iscript_parallel_backend(base_dir, opt$ncores)
}


initialize_config_file <- function(opt){
	config_filename <- opt$options$config_filename
	if(is.null(config_filename)){
		return(NULL)
	}
	if(!file.exists(config_filename)){
		cat("ERROR: Config file '", config_filename, "' does not exist.\n")
		stop(1)
	}

	tryCatch({
		configuration <- fromJSON(file=config_filename)
	}, error=function(e){
		cat(paste(
			"ERROR: Unable to parse configuration file '", config_filename, "'\n",
			"failed with the following error:\n",
			e, sep=""))
	})
	return(configuration)
}

initialize_analysis_scripts_from_options <- function(opt, base_dir){
	if("script" %in% names(opt$options)){
		if(file.exists(opt$options$script)){
			analysis_scripts <- c(opt$options$script)
		} else if(file.exists(file.path(base_dir, opt$options$script))){
			analysis_scripts <- c(file.path(base_dir, opt$options$script))
		} else {
			stop(paste("Analysis script '", opt$options$script, "', does not exist.", sep=""))
		}
	} else if("analysis_dir" %in% names(opt$options)) {
		if(substr(opt$options$analysis_dir,1,1) == "/"){
			analysis_dir <- opt$options$analysis_dir
		} else {
			analysis_dir <- paste(base_dir, opt$options$analysis_dir, sep="/")
		}
		if(!file.exists(analysis_dir)){
			cat("ERROR: Analysis script directory, '", analysis_script_dir, "', does not exist.\n", sep="")
			stop(1)
		}
		analysis_scripts <- dir(analysis_dir, "*.R$", full.names=TRUE, recursive=TRUE)
	} else {
		analysis_scripts <- c()
	}
	analysis_scripts
}

add_command_line_options_to_configuration <- function(
	configuration,
	opt,
	base_dir
) {
	configuration$build_dir <- opt$options$build_dir

  
	sample_sources <- alply(get_sample_sources(opt$args), 1, function(ss){
		list(
			database_path=as.character(ss$fname),
			id=as.character(ss$sample_source))
	})
	names(sample_sources) <- NULL
	attributes(sample_sources) <- NULL

	analysis_scripts <- initialize_analysis_scripts_from_options(opt, base_dir)

	if(length(analysis_scripts)){
		configuration$sample_source_comparisons <- c(
			configuration$sample_source_comparisons,
			list(list(
				sample_sources=sample_sources,
				analysis_scripts=analysis_scripts)))
	}

	configuration
}

add_configuration_options_to_options <- function (
  configuration,
  opt
) {
  
  if (! is.null(configuration$output_dir)){
    opt$output_dir <- configuration$output_dir
  }
  opt
}


initialize_output_dir <- function(work_dir, opt, ss_cmp){
	# if the output dir is a relative path, make it relative to the
	# work_dir, which is the path from which the script was executed
	if(substr(opt$options$output_dir,1,1) == "/"){
		output_dir <- opt$options$output_dir
	} else {
		output_dir <- file.path(work_dir, opt$options$output_dir)
	}

  if(! opt$options$out_disable_ss_names){
	  file.path(
		  output_dir,
		  paste(
			  llply(ss_cmp$sample_sources, function(ss) ss$id),
			  sep="/", collapse="_"))
  }
  output_dir
}

summarize_configuration <- function(
	output_dir,
	output_formats,
	sample_source_comparison){

	cat(
		"Sample Source Comparison:\n",
		"  Output Directory <- '", output_dir, "'\n", sep="")
	cat(
		"    Output Formats <- ", paste(output_formats$id, collapse=", "), "\n\n",
		sep="")

	cat("  Sample Sources:\n")
	l_ply(sample_source_comparison$sample_sources, function(ss) {
		cat("  ", ss$id, " <- ", ss$database_path, "\n", sep="")
	})
	cat("\n  Analysis_scripts:\n")
	cat(" ", paste(sample_source_comparison$analysis_scripts, sep="", colapse="\n "))
	cat("\n")
	iscript_summarize_configuration(sample_source_comparison)
}

initialize_sample_sources <- function(ss_cmp){
	sample_sources <- ldply(ss_cmp$sample_sources, function(ss){

		if(is.na(ss$id) || is.null(ss$id) || length(ss$id) == 0){
			stop(paste(
				"ERROR: The 'id' for the sample source is empty.", sep=""))
		}

		if(is.na(ss$id) || is.null(ss$id) || length(ss$database_path) == 0){
			stop(paste(
				"ERROR: The 'database_path' for sample source '", ss$id, "' is empty.", sep=""))
		}

		if(substr(ss$database_path,1,1) == "/"){
			database_path <- ss$database_path
		} else {
			database_path <- file.path(getwd(), ss$database_path)
		}

		if(!file.exists(database_path)){
			stop(paste(
				"ERROR: The database path '", database_path, "' ",
				"for sample source '", ss$id, "' does not exist.", sep=""))
		}

		#TODO filter null values because they lead to errors
		data.frame(
			fname=database_path,
			sample_source=ss$id,
			ss[names(ss) != "id"])

	})

	if(nrow(sample_sources) != length(unique(sample_sources$sample_source))){
		cat("ERROR: The sample sources must have unique identifiers:\n")
		print(sample_sources)
		stop()
	}

	iscript_sample_sources(sample_sources)
	sample_sources
}


initialize_output_formats <- function(opt, ss_cmp){
	if("output_formats" %in% names(ss_cmp)){
		output_formats <- get_output_formats_from_comparison(
			ss_cmp, all_output_formats)
	} else {

		#only using the output formats from the command line if they are
		#not specified in the analysis configuration. Note this is a
		#little goofy, feel free to refactor...
		output_formats <- get_output_formats(opt$options, all_output_formats)
	}


	# If we're generating the features website for analysis we better
	# have the plots in the formats that we expect.
	if(opt$options$generate_website){
		opt$options$output_web_raster <- T
		opt$options$output_web_icon <- T
	}

	if(nrow(output_formats) == 0){
		stop(paste(
			"ERROR: No output formats were specified.  To specify an output format, either add the appropriate tags the sample source configuration file or specify one or more output formats on the command line on the command line."))
	}

	iscript_output_formats(output_formats)

	#Setup footer specification in output_formats
	output_formats$add_footer <- output_formats$accepts_footer & opt$options$add_footer
	iscript_setup_add_footer_to_output_formats(opt$options$add_footer)
	output_formats
}


initialize_database_configuration <- function(opt){
	#Setup db_cache_size
	database_configuration <- list(
		db_cache_size = opt$options$db_cache_size
	)
	iscript_db_cache_size(database_configuration$db_cache_size)
	database_configuration
}

initialize_density_estimation_kernel <- function(opt){
	#Setup density estimation options
	general_kernel_adjust <- opt$options$general_kernel_adjust
	iscript_general_kernel_adjust(general_kernel_adjust)
	general_kernel_adjust
}

parse_analysis_scripts <- function(base_dir, analysis_scripts){
	#Read in all the features analysis scripts
	iscript_source_scripts(base_dir, analysis_scripts)

	# This is a vector of FeaturesAnalysis objects that are defined each features analysis script
	feature_analyses <- c()
	num_feature_analyses_before <- 0
	for(analysis_script in analysis_scripts){

		# parse all the analysis scripts
		if(file.exists(analysis_script)){
			# analysis script is ok
		} else if (file.exists(file.path(base_dir, analysis_script))){
			analysis_script <- file.path(base_dir, analysis_script)
		} else {
			cat(paste(
				"ERROR: The features analysis script '",
				analysis_script,"' does not exist\n", sep=""))
		}

		tryCatch({
			source(analysis_script, local=T)
		}, error=function(e){
			cat(paste(
				"ERROR: Failed to parse the Features Analysis '",
				analysis_script,"' with the following error message:\n", e, sep=""))
		})

		# assign the filename to each feature analysis
		num_new_feature_analyses = length(feature_analyses) - num_feature_analyses_before
		for(feature_analysis in
			feature_analyses[
				seq(num_feature_analyses_before+1, length.out=num_new_feature_analyses)]){
			feature_analysis@filename <- analysis_script
		}
		num_feature_analyses_before <- length(feature_analyses)
	}
	feature_analyses
}

run_feature_analyses <- function(
	feature_analyses,
	sample_sources,
	opt,
	base_dir,
	output_dir,
	output_formats
) {
	#Run all the features analysis scripts
	iscript_run_feature_analyses(output_dir)
	if(!opt$options$dry_run){
		for(features_analysis in feature_analyses){
			cat(paste("Features Analysis: ", features_analysis@id, "\n", sep=""))

			# set current working directory the base_dir so that way scripts
			# can reference other scripts in a canonical way

			setwd(base_dir)
			tryCatch({
				features_analysis@run(
					features_analysis,
					sample_sources,
					output_dir,
					output_formats)
			}, error=function(e){
				cat(paste("ERROR: Failed to run the Features Analysis '", features_analysis@id, "' with the follwing error message:\n", e, sep=""))
			})
			cat("\n")
		}
	}
}

generate_webpages <- function(opt, output_dir){
	if(opt$options$generate_website){
		cat("Generating Features website: '", output_dir, "'\n", sep="")
		generate_sample_source_comparison_webpages(output_dir)
	}
}


########################################################

base_dir <- initialize_base_dir()
work_dir <- initialize_work_dir()
initialize_iscript(base_dir)
initialize_packages_before_options(base_dir)
initialize_method_scripts_before_options(base_dir)

opt <- initialize_command_line_options()
check_basic_input(opt)

initialize_packages(opt, base_dir)
initialize_method_scripts(base_dir)
initialize_parallel_backend(opt, base_dir)

database_configuration <- initialize_database_configuration(opt)
general_kernel_adjust <- initialize_density_estimation_kernel(opt)

configuration <- initialize_config_file(opt)
configuration <- add_command_line_options_to_configuration(
	configuration,
	opt,
	base_dir)

opt <- add_configuration_options_to_options(configuration, opt)

l_ply(
	.data=configuration$sample_source_comparisons,
	.parallel=use_parallel,
	.fun=function(ss_cmp){

	sample_source_output_dir <- initialize_output_dir(work_dir, opt, ss_cmp)
	output_formats <- initialize_output_formats(opt, ss_cmp)

	summarize_configuration(sample_source_output_dir, output_formats, ss_cmp)

	sample_sources <- initialize_sample_sources(ss_cmp)

	feature_analyses <- parse_analysis_scripts(base_dir, ss_cmp$analysis_scripts)

	run_feature_analyses(
		feature_analyses,
		sample_sources,
		opt,
		base_dir,
		sample_source_output_dir,
		output_formats)

	#generate_webpages(opt, sample_source_output_dir)

})




