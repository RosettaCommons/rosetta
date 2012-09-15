# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


date_code <- function(d=NA){
	# reference http://www.r-cookbook.com/node/17
	if(is.na(d)) d <- Sys.Date()
	pattern <- '20([[:digit:]]{2})-([[:digit:]]{2})-([[:digit:]]{2})'
	paste(
		sub(pattern, '\\1', d), sub(pattern, '\\2', d), sub(pattern, '\\3', d),
		sep="")
}


# Save a data.frame as a table. For each output format,
# generate a table and put in the output directory
save_tables <- function(
	features_analysis,
	table,
	table_id,
	sample_sources,
	output_dir,
	output_formats,
	...
) {
	tryCatch(table_id, error=function(e){
		stop(paste(
			"ERROR: Unable to save the table because ",
			"the 'table_id' is not specified.\n", e, sep=""))
	})

	tryCatch(features_analysis, error=function(e){
		stop(paste(
			"ERROR: Unable to save the table '", table_id,"' ",
			"because the specified 'features_analysis' is not valid.\n",
			e, sep=""))
	})

	tryCatch(sample_sources, error=function(e){
		stop(paste(
			"ERROR: Unable to save the table '", table_id, "' ",
			"because the specified 'sample_sources' is not valid.\n",
			e, sep=""))
	})

	if(nrow(sample_sources)==0){
		stop(paste(
			"ERROR: Unable to save the table '", table_id, "' ",
			"because no sample_sources were specified.\n", e, sep=""))
	}

	tryCatch(output_dir, error=function(e){
		stop(paste(
			"ERROR: Unable to save the table '", table_id, "' ",
			"because the specified 'output_dir' ",
			"is not a valid variable.\n",
			e, sep=""))
	})

	tryCatch(output_formats, error=function(e){
		stop(paste(
			"ERROR: Unable to save the table '", table_id, "' ",
			"because the 'output_formats' parameter is not valid.\n",
			e, sep=""))
	})
	table_formats <- output_formats[output_formats$type == "table",]

	if(nrow(table_formats)==0){
		stop(paste(
			"ERROR: Unable to save the table '", table_id, "' ",
			"because no output formats were specified.", sep=""))
	}

	a_ply(table_formats, 1, function(fmt){
		full_output_dir <- file.path(output_dir, features_analysis@id, fmt$id)
		if(!file.exists(full_output_dir)){
			dir.create(full_output_dir, recursive=TRUE)
		}
		date <- date_code()
		fname <- paste(table_id, date, sep="_")
		full_path <- file.path(full_output_dir, paste(fname, fmt$extension, sep=""))
		type <- substring(fmt$extension, 2)
		cat("Saving Table of type ", type, ": ", full_path, sep="")
		timing <- system.time({
			tryCatch({
				print(xtable(
					table),
					file=full_path,
					type=type,
					...)
			}, error=function(e){
				cat("\n")
				cat(paste(
					"ERROR: Generating and saving the table:\n",
					e, sep=""))
			})
		})
		cat(" ... ", as.character(round(timing[3],2)), "s\n", sep="")
	})
}


