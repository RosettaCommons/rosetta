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
libraries <- c("optparse", "ggplot2", "RSQLite", "logspline", "plotrix", "plyr", "polynom")
load_packages(libraries)
iscript_libraries(libraries)

iscript_base_dir(base_dir)

includes <- c(
	"scripts/methods/sample_sources.R",
	"scripts/methods/methods.R",
	"scripts/methods/density_estimation.R",
	"scripts/methods/ggplot2_geom_indicator.R",
	"scripts/methods/ggplot2_scales.R",
	"scripts/methods/instancer.R",
	"scripts/methods/output_formats.R",
	"scripts/methods/coordinate_normalizations.R",
	"scripts/methods/vector_math.R",
	"scripts/methods/color_palettes.R")
for(inc in includes) source(paste(base_dir, inc, sep="/"))
iscript_includes(includes)

option_list <- list(
	make_option(c("-v", "--verbose"), action="store_true", default=FALSE, dest="verbose",
							help="Print extra output [default]"),
	make_option(c("-i", "--script"), action="store", type="character", default=NULL, dest="script",
							help="Path to a single analysis script.  [Default \"%default\"]"),
	make_option(c("-a", "--analysis_dir"), type="character", default="scripts/analysis", dest="analysis_dir",
							help="Directory where the analysis scripts are located. The supplied directory is searched recursively for files of the form \"*.R\".  [Default \"%default\"]"),
	make_option(c("-s", "--sample_source_dir"), type="character", default="sample_sources", dest="sample_source_dir",
							help="Directory where the feature databases are located. The supplied directory is searched recursively for SQLite database files of the form \"features_<sample_source_id>_<date_code>.db3\".  [Default \"%default\"]"),
	make_option(c("-o", "--output_dir"), type="character", default="build", dest="output_dir",
							help="Directory where the output plots and statistics will be generated.  [Default \"%default\"]"),
	make_option(c("--output_web_raster"), action="store_true", type="logical", default=TRUE, dest="output_web_raster",
							help="Generate output plots suitable for the web in .png format.  [Default \"%default\"]"),
	make_option(c("--output_web_vector"), action="store_true", type="logical", default=FALSE, dest="output_web_vector",
							help="Generate output plots suitable for the web in .svg format.  [Default \"%default\"]"),
	make_option(c("--output_slide_raster"), action="store_true", type="logical", default=FALSE, dest="output_slide_raster",
							help="Generate output plots suitable for slides in .png format.  [Default \"%default\"]"),
	make_option(c("--output_slide_vector"), action="store_true", type="logical", default=FALSE, dest="output_slide_vector",
							help="Generate output plots suitable for slides in .svg format.  [Default \"%default\"]"),
	make_option(c("--output_slide_pdf"), action="store_true", type="logical", default=FALSE, dest="output_slide_pdf",
							help="Generate output plots suitable for slides in .pdf format.  [Default \"%default\"]"),
	make_option(c("--output_print_raster"), action="store_true", type="logical", default=FALSE, dest="output_print_raster",
							help="Generate output plots suitable for printing in .png format.  [Default \"%default\"]"),
	make_option(c("--output_print_vector"), action="store_true", type="logical", default=FALSE, dest="output_print_vector",
							help="Generate output plots suitable for printing in .svg format.  [Default \"%default\"]"),
	make_option(c("--output_print_pdf"), action="store_true", type="logical", default=FALSE, dest="output_print_pdf",
							help="Generate output plots suitable for printing in .pdf format.  [Default \"%default\"]"),
	make_option(c("--output_huge_raster"), action="store_true", type="logical", default=FALSE, dest="output_huge_raster",
							help="Generate output plots suitable for hugeing in .png format.  [Default \"%default\"]"),
	make_option(c("--output_huge_vector"), action="store_true", type="logical", default=FALSE, dest="output_huge_vector",
							help="Generate output plots suitable for hugeing in .svg format.  [Default \"%default\"]"),
	make_option(c("--output_huge_pdf"), action="store_true", type="logical", default=FALSE, dest="output_huge_pdf",
							help="Generate output plots suitable for hugeing in .pdf format.  [Default \"%default\"]"),
	make_option(c("--db_cache_size"), action="store_true", type="integer", default=10000, dest="db_cache_size",
							help="Number of 1k pages of cache to use for database queries.  [Default \"%default\"]"))

opt <- parse_args(OptionParser(option_list=option_list), positional_arguments=TRUE)


#Setup sample sources
if(length(opt$args) > 0){
  #use feature databases specified as positional arguments
	cat("Using supplied feature databases", opt$args, "\n")
	data_sources <- opt$args
} else {
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
		"       source(\"compare_sample_sources_iscript.R\", # This will re-run the last './compare_sample_sources.R' interactively",
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
iscript_sample_sources(sample_sources)

#Setup analysis_scripts
if("script" %in% names(opt$options)){
	if(file.exists(opt$options$script)){
		analysis_scripts <- c(opt$options$script)
	} else {
		stop(paste("Analysis script '", opt$options$script, "', does not exist."))
	}
} else {
	if(substr(opt$options$analysis_dir,1,1) == "/"){
		analysis_dir <- opt$options$analysis_dir
	} else {
		analysis_dir <- paste(base_dir, opt$options$analysis_dir, sep="/")
	}
	if(!file.exists(analysis_dir)){
		cat("ERROR: Analysis script directory, '", analysis_script_dir, "', does not exist.\n")
		stop(1)
	}
	analysis_scripts <- dir(analysis_dir, "*.R$", full.names=TRUE, recursive=TRUE)
}
cat("\n")
cat("Running the following analysis scripts:\n   ")
cat(paste(analysis_scripts, sep="", colapse="\n  "))
cat("\n")



#Validate ouput_dir
if(!file.exists(opt$options$output_dir)){
	print(paste("Creating output directory: '",opt$options$output_dir,"'...",sep=""))
	dir.create(opt$options$output_dir)
	if(!file.exists(opt$options$output_dir)){
		print("ERROR: Unable to create output directory.")
		stop(1)
	}
}
output_dir <- opt$options$output_dir
iscript_output_dir(output_dir)

#Setup output formats
output_formats <- get_output_formats(opt$options)
iscript_output_formats(output_formats)

#Setup db_cache_size
db_cache_size <- opt$options$db_cache_size
iscript_db_cache_size(db_cache_size)

iscript_scripts(analysis_scripts)
l_ply(analysis_scripts, function(analysis_script){
	cat(paste("Begin running '", analysis_script,"'.\n", sep=""))
	tryCatch(source(analysis_script), error=function(e){
		cat(paste("ERROR: The analysis script '",analysis_script,"' failed with error:\n",e,sep=""))
	})
	cat("\n")
})

