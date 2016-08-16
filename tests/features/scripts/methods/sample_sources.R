# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


#TODO THIS FUNCTION IS BROKEN
#Each feature database contains a table sample_source with containing
#meta data for the sample source.
get_sample_source_meta_data <- function(fname){
	if(!file.exists(fname)){
		stop(paste("ERROR: The sample source SQLite database, '", fname, "', does not exist.", sep=""))
	}

	tryCatch({
		con <- dbConnect("SQLite", fname)
	}, error=function(e){
		stop(paste("ERROR: Unable to connect with the sample source SQLite database, '", fname, "'.\nERROR: ", e, sep=""))
	})

	tryCatch({
		df <- dbGetQuery(con, "SELECT * FROM sample_source LIMIT 1;")
	}, error=function(e){
		stop(paste("ERROR: Querying the 'sample_source' table in the sample source SQLite database '", fname, "' failed.\nERROR: ", e, sep=""))
	})

	if(nrow(df)==0){
		stop(paste("ERROR: The sample source SQLite database, '", fname, "', has no rows in the sample_source table."))
	}

	dbDisconnect(con)
	# get the meta-data from the first row of the sample sources table
	return( df[1,] )
}

# this searchers for feature databases given a base directory
get_data_sources <- function(sample_source_dir, extension=".db3"){
	if(!file.exists(sample_source_dir)){
		stop(paste("ERROR: Sample source directory, '", sample_source_dir, "', does not exist.", sep=""))
	}

	full_sample_source_dir <- paste(getwd(), sample_source_dir, sep="/")
	data_sources <- dir(full_sample_source_dir, paste(extension, "$", sep=""), full.names=TRUE, recursive=TRUE)

	if(length(data_sources) < 1){
		warning(paste("No sample sources were found in the '", full_sample_source_dir, "' directory.", sep=""))
	}
	data_sources
}

# given a list of feature database filenames, return a data.frame with
# information each one that can be used in query_sample_sources(...)
get_sample_sources <- function(data_sources){
	ss <- ldply(data_sources, function(fname){
		if(fname == ""){
			stop("Unable to get meta data from sample source because no database file name was supplied.")
		}
		sample_source <- factor( strsplit(
						strsplit( basename(fname), "^features_")[[1]][2],".db3$")[[1]])
		if(is.na(sample_source) || is.null(sample_source) || sample_source==""){
			stop(paste("Unable to get 'sample_source_id' from sample source database file name '",fname,"', verify that it is of the form 'features_<sample_source_id>.db3'", sep=""))
		}
		#cbind(data.frame(fname, sample_source), get_sample_source_meta_data(fname))
		data.frame(fname, sample_source)
	})
	rownames(ss) <- ss$sample_source
	ss
}

sample_source_titles <- function(sample_sources){
	paste(laply(sample_sources, function(ss) ss$sample_source), collapse=" ")
}

