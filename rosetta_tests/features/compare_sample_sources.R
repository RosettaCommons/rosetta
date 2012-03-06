#! /usr/bin/Rscript --vanilla
# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

# Get path to 'rosetta_tests/features', where this script lives
command_args <- commandArgs(trailingOnly = FALSE)
base_dir <- dirname(substring(command_args[grep("--file=", command_args)], 8))
if(length(base_dir) == 0){
	base_dir = getwd()
}
source(paste(base_dir, "scripts/methods/generate_compare_sample_sources_iscript.R", sep="/"))
iscript_header()

# load_packages() will help the user to install the packages if they are missing
source(paste(base_dir, "scripts/methods/load_packages.R", sep="/"))
libraries <- c("reshape", "plyr", "optparse", "ggplot2", "RSQLite", "logspline", "plotrix", "polynom", "xtable")
load_packages(libraries)
iscript_libraries(libraries)

iscript_base_dir(base_dir)

includes <- c(
	"scripts/methods/sample_sources.R",
	"scripts/methods/methods.R",
	"scripts/methods/density_estimation.R",
	"scripts/methods/features_analysis.R",
	"scripts/methods/analysis_manager.R",
	"scripts/methods/ggplot2_geom_indicator.R",
	"scripts/methods/ggplot2_scales.R",
	"scripts/methods/instancer.R",
	"scripts/methods/output_formats.R",
	"scripts/methods/save_plots.R",
	"scripts/methods/coordinate_normalizations.R",
	"scripts/methods/comparison_statistics.R",
	"scripts/methods/vector_math.R",
	"scripts/methods/color_palettes.R")
for(inc in includes) source(paste(base_dir, inc, sep="/"))
iscript_includes(includes)


#Setup options
option_list <- list()

# Run level options
option_list <- c(option_list, list(
	make_option(c("-v", "--verbose"), action="store_true", default=FALSE, dest="verbose",
							help="Print extra output [default]"),
	make_option(c("--dry_run"), action="store_true", type="logical", default=FALSE, dest="dry_run",
							help="Debug the analysis scripts but do not run them.  [Default \"%default\"]"),
	make_option(c("--db_cache_size"), action="store_true", type="integer", default=10000, dest="db_cache_size",
							help="Number of 1k pages of cache to use for database queries.  [Default \"%default\"]"),
	make_option(c("--add_footer"), action="store_true", type="logical", default=TRUE, dest="add_footer",
							help="Add footer to plots saying the analysis script and run date.")))


# Which analysis scripts should be run
option_list <- c(option_list, list(
	make_option(c("-i", "--script"), action="store", type="character", default=NULL, dest="script",
							help="Path to a single analysis script.  [Default \"%default\"]"),
	make_option(c("-a", "--analysis_dir"), type="character", default="scripts/analysis", dest="analysis_dir",
							help="Directory where the analysis scripts are located. The supplied directory is searched recursively for files of the form \"*.R\".  [Default \"%default\"]")))

# Where should the results be stored?
option_list <- c(option_list,
								 list(
	make_option(c("-o", "--output_dir"), type="character", default="build", dest="output_dir",
							help="Directory where the output plots and statistics will be generated.  [Default \"%default\"]")))

option_list <- c(option_list, make_output_formats_options_list(all_output_formats))


# Analysis manager options
option_list <- c(option_list, list(
	make_option(c("--analysis_manager_db"), type="character", default="analysis_manager.db3", dest="analysis_manager_db",
							help="Store information about the results of the features analysis in the analysis manager database.  [Default \"<output_dir>/%default\"]"),
	make_option(c("--store_plots_in_analysis_manager_db"), action="store_true", type="logical", default=FALSE, dest="store_plots_analysis_manager_db",
							help="Should the plots themselves be stored in the the analysis manager database?  [Default \"%default\"]")))

opt <- parse_args(OptionParser(option_list=option_list), positional_arguments=TRUE)

#Setup sample sources
data_sources <- opt$args
if(length(data_sources) == 0){
	cat(
		"ERROR: No sample_source databases were supplied",
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
sample_sources <- get_sample_sources(data_sources)

#setup ouput_dir
output_dir <- opt$options$output_dir
sample_source_output_dir <- file.path(
	output_dir,
	paste(sample_sources$sample_source, collpase="/", sep=""))

if(!file.exists(sample_source_output_dir)){
	cat("Creating output directory: '", sample_source_output_dir, "' ...\n", sep="")

	dir.create(sample_source_output_dir, recursive=TRUE)

	if(!file.exists(sample_source_output_dir)){
		print("ERROR: Unable to create output directory.")
		stop(1)
	}

	file.copy(
		from=file.path(base_dir, "scripts/methods/features_web"),
		to=file.path(sample_source_output_dir),
		recursive=TRUE)
}
iscript_output_dir(output_dir)



#Setup analysis manager
analysis_manager_db_path <- paste(
	opt$options$output_dir,
	paste(sample_sources$sample_source, collapse="/"),
	opt$options$analysis_manager_db, sep="/")
iscript_setup_analysis_manager(analysis_manager_db_path)
analysis_manager_con <- initialize_analysis_manager_db(
	analysis_manager_db_path)

add_sample_sources_to_analysis_manager(analysis_manager_con, sample_sources)
iscript_sample_sources(sample_sources)
cat("Comparing the following sample sources:\n")
a_ply(sample_sources,1, function(ss){
	cat("   ", as.character(ss$sample_source), " <- ", as.character(ss$fname), "\n", sep="")
})

#Setup analysis_scripts
if("script" %in% names(opt$options)){
	if(file.exists(opt$options$script)){
		analysis_scripts <- c(opt$options$script)
	} else {
		stop(paste("Analysis script '", opt$options$script, "', does not exist.", sep=""))
	}
} else {
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
}
cat("\n")
cat("Running the following analysis scripts:\n   ")
cat(paste(analysis_scripts, sep="", colapse="\n  "))
cat("\n")


#Setup output formats
output_formats <- get_output_formats(opt$options, all_output_formats)
iscript_output_formats(output_formats)

#Setup footer specification in output_formats
output_formats$add_footer <- output_formats$accepts_footer & opt$options$add_footer
iscript_setup_add_footer_to_output_formats(opt$options$add_footer)


#Setup db_cache_size
db_cache_size <- opt$options$db_cache_size
iscript_db_cache_size(db_cache_size)

#Read in all the features analysis scripts
iscript_source_scripts(analysis_scripts)

# This is a vector of FeaturesAnalysis objects that are defined each features analysis script
feature_analyses <- c()
num_feature_analyses_before <- 0
for(analysis_script in analysis_scripts){

	# parse all the analysis scripts
	tryCatch(source(analysis_script), error=function(e){
		cat(paste("ERROR: Failed to parse the Features Analysis '",analysis_script,"' with the following error message:\n",e,sep=""))
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


#Run all the features analysis scripts
iscript_run_feature_analyses()
if(!opt$options$dry_run){
	for(features_analysis in feature_analyses){
		cat(paste("Features Analysis: ", features_analysis@id, "\n", sep=""))

		tryCatch({
			add_features_analysis_to_analysis_manager(
				analysis_manager_con, features_analysis)
		}, error=function(e){
			cat(paste("ERROR: Failed to store the Features Analysis '", features_analysis@id, "' in the analysis manager and it failed wih the following error message:\n", e, sep=""))
		})

		tryCatch({
			features_analysis@run(features_analysis)
		}, error=function(e){
			cat(paste("ERROR: Failed to run the Features Analysis '", features_analysis@id, "' with the follwing error message:\n", e, sep=""))
		})
		cat("\n")
	}
}



# close connection to the analysis manager
result <- dbDisconnect(analysis_manager_con)
if(!result){
	cat("ERROR: Failed to close the analysis manager database connection.\n")
}
