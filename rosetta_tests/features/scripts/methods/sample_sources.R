# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#Each feature database contains a table sample_source with containing
#meta data for the sample source.
get_sample_source_meta_data <- function(fname){
	cat(paste("getting meta data from fname '", fname,"'\n",sep=""))
	con <- dbConnect("SQLite", fname)
	df <- dbGetQuery(con, "SELECT * FROM sample_source;")
	dbDisconnect(con)

	# get the meta data from the first row of the sample sources table
	return( df[1,] )
}

# this searchers for feature databases given a base directory
get_data_sources <- function(sample_source_dir, extension=".db3"){
	if(!file.exists(sample_source_dir)){
		stop(paste("ERROR: Sample source directory, '", sample_source_dir, "', does not exist."))
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
			stop(paste("Unable to get 'sample_source_id' from sample source database file name '",fname,"', verify that it is of the form 'features_<sample_source_id>.db3'"))
		}
		data.frame(fname, sample_source, get_sample_source_meta_data(fname))
	})
	rownames(ss) <- ss$sample_source
	ss
}

sample_source_titles <- function(sample_sources){
	paste(laply(sample_sources, function(ss) ss$sample_source), collapse=" ")
}


add_sample_sources_to_analysis_manager <- function(
	con, sample_sources){

	sql <- "INSERT INTO sample_sources (?,?,?,?);"
	dbBeginTransaction(con)
	l_ply(sample_sources, function(ss) dbGetPreparedQuery(sql, ss))
	dbCommit(con)

}

add_feature_reporters_to_analysis_manager <- function(
		con, sample_sources){

	dbBeginTransaction(con)
	l_ply(sample_sources, function(ss) {
		dbSendPreparedQuery(con,
			"INSERT OR IGNORE INTO sample_sources VALUES (?,?,?,?);", ss);
		#TODO sanitize ss$fname and ss$sample_source
		sql <- paste("
ATTACH DATABASE '", ss$fname, "' AS ss;

INSERT OR IGNORE INTO feature_reporters
	SELECT * FROM ss.feature_reporters;

INSERT OR IGNORE INTO feature_analysis_tables
	SELECT * FROM ss.feature_analysis_tables;

INSERT OR IGNORE INTO sample_source_feature_reporters
	SELECT
		'", ss$sample_source, "',
		feature_reporter.feature_reporter_id
	FROM feature_reporters AS feature_reporter;

DETACH DATABASE ss;", sep="")

		dbSendQuery(con, sql)
	})
	dbCommit(con)
}

