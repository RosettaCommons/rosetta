# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


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



add_sample_sources_to_analysis_manager <- function(
		con, sample_sources){

	a_ply(sample_sources, 1, function(ss) {
		#TODO sanitize ss$fname and ss$sample_source and ss$description
		sql <- paste("INSERT OR REPLACE INTO sample_sources VALUES ('",
			paste(as.character(ss$sample_source), as.character(ss$fname),
				as.character(ss$description), sep="', '"), "');", sep="")
		dbGetQuery(con, sql)

		ss_con <- dbConnect(engine, as.character(ss$fname))
		df <- dbGetQuery(ss_con, "SELECT count(*) AS has_features_reporters_table FROM sqlite_master WHERE name='features_reporters';")
		dbDisconnect(ss_con)
		if(df$has_features_reporters_table[1] == 0){
			#print(WARNING: The feature database '", ss$fname, "' does not have a features_reporter table. This may be because the code you used to extract the features database is outdated."))
		} else {
			dbGetQuery(con, paste("ATTACH DATABASE '", ss$fname, "' AS ss;", sep=""))
			dbGetQuery(con, "INSERT OR IGNORE INTO features_reporters
SELECT * FROM ss.features_reporters;")
			dbGetQuery(con, paste("INSERT OR IGNORE INTO sample_source_features_reporters
		SELECT '", ss$sample_source, "', features_reporter.features_reporter_type_name
		FROM features_reporters AS features_reporter;", sep=""))
			dbGetQuery(con, "DETACH DATABASE ss;")
		}
#		dbDisconnect(ss_con)
	})
}

