#-*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

engine <- SQLite();

sample_rows <- function(df, n=100){
	if( n >= nrow(df) ) return( df )
	df[sample(nrow(df),n),]
}

setup_output_directory <- function(output_dir){
	if(!file.exists(file.path(output_dir))){
		dir.create(file.path(output_dir), recursive=TRUE)
	}
}

check_setup <- function(){
#	tryCatch(sample_sources, error=function(e){
#		stop("ERROR: The variable 'sample_sources' is not defined. See compare_sample_sources.R")
#	})
#	if(!is.data.frame(sample_sources)){
#		stop("ERROR: The variable 'sample_sources' is not a data frame.")
#	}
#	if(nrow(sample_sources) < 1){
#		stop("ERROR: The variable 'sample_sources' contains no sample sources.")
#	}
#
#	tryCatch(sample_source_output_dir, error=function(e){
#		stop("ERROR: The variable 'sample_source_output_dir' is not defined. See compare_sample_sources.R")
#	})
#	tryCatch(output_formats, error=function(e){
#		stop("ERROR: The variable 'output_formats' is not defined. See compare_sample_sources.R")
#	})
#	tryCatch(database_configuration$db_cache_size, error=function(e){
#		stop("ERROR: The variable 'database_configuration$db_cache_size' is not defined. See compare_sample_sources.R")
#	})
}


set_db_cache_size <- function(con, cache_size){
	if(is.null(cache_size)){
		stop("ERROR: unable to set database cache size because the cache_size is null.")
	}
	dbGetQuery(con,
		paste("PRAGMA cache_size=",as.integer(cache_size),";",sep=""))
}

query_sample_source <- function(
  sample_source,
  sele,
  sele_args_frame = NULL,
  cache_size=database_configuration$db_cache_size,
  char_as_factor=T
  ){
  tryCatch(sele,error=function(e){
  	cat("ERROR: The select statement is not defined.\n")
	})
  
  #Get the features from the sample_source
	tryCatch(c(sample_source),error=function(e){
		cat("ERROR: The specified sample source is not defined.\n")
	})

	cat("loading:", as.character(sample_source$sample_source), "... ")
	if( is.na(sample_source$sample_source[1]) ){
		stop("Specified sample source is not defined")
	}
	con <- dbConnect(engine, as.character(sample_source$fname))
	set_db_cache_size(con, cache_size);

	#Allow select statements to be prefaced with arbitrary statements.
	#This allows the creation of temporary tables, indices, etc.
	sele_split <- paste(strsplit(sele, ";\\W*", perl=TRUE)[[1]], ";", sep="")
	l_ply(sele_split[-length(sele_split)], function(sele){
      dbGetQuery(con, sele)
	})

	last_stmt <- sele_split[length(sele_split)] #Real statement for data
  if (! is.null(sele_args_frame)){
    features <- dbGetPreparedQuery(con, sele, bind.data=sele_args_frame)
	} 
  else {
		features <- dbGetQuery(con, last_stmt)
  }

	dbDisconnect(con)
	#cat(as.character(round(timing[3], 2)),"s\n")

	if(nrow(features)==0){
		cat("WARNING: The following query returned no rows:\n")
		cat(sele)
	}
  if(char_as_factor){
	  for(col in names(features)){
		  if(is.character(features[,col])){
			  features[,col] <- factor(features[,col])
		  }
	  }
	}
	features
}

query_sample_sources <- function(
	sample_sources,
	sele,
  sele_args_frame = NULL,
	cache_size=database_configuration$db_cache_size,
  char_as_factor=T
  ){
	tryCatch(sele,error=function(e){
		cat("ERROR: The select statement is not defined.\n")
	})
	features <- ddply(sample_sources, c("sample_source"), function(ss){
		tryCatch(c(ss),error=function(e){
			cat("ERROR: The specified sample source is not defined.\n")
		})

		cat("loading:", as.character(ss$sample_source), "... ")
		if( is.na(ss$sample_source[1]) ){
			stop("Specified sample source is not defined")
		}
		timing <- system.time({
			con <- dbConnect(engine, as.character(ss$fname))
			set_db_cache_size(con, cache_size);

			#Allow select statements to be prefaced with arbitrary statements.
			#This allows the creation of temporary tables, indices, etc.
			sele_split <- paste(strsplit(sele, ";\\W*", perl=TRUE)[[1]], ";", sep="")
			l_ply(sele_split[-length(sele_split)], function(sele){
          dbGetQuery(con, sele)
			})

			last_stmt <- sele_split[length(sele_split)] #Real statement for data
      if (! is.null(sele_args_frame)){
        df <- dbGetPreparedQuery(con, sele, bind.data=sele_args_frame)
			} 
      else {
		    df <- dbGetQuery(con, last_stmt)
      }
      
			dbDisconnect(con)
		})
		cat(as.character(round(timing[3], 2)),"s\n")
		df
	})
	if(nrow(features)==0){
		cat("WARNING: The following query returned no rows:\n")
		cat(sele)
	}
  if(char_as_factor){
	  for(col in names(features)){
		  if(is.character(features[,col])){
			  features[,col] <- factor(features[,col])
		  }
	  }
	}
	features
}

query_sample_sources_against_ref <- function(
	sample_sources,
	sele,
  sele_args_frame = NULL,
	cache_size=database_configuration$db_cache_size,
  char_as_factor=T
  ){
	tryCatch(sele,error=function(e){
		cat("ERROR: The select statement ", sele, " is not defined.\n")
	})

	if(nrow(sample_sources) < 2) {
		stop(paste("Please provide 2 or more sample sources.

sample_sources provided:
	'", paste(sample_sources$sample_source_ids, collapse= "',\n\t'"), "'

This function will execute the query for all but the first sample source
	The first sample source will be available as a database named 'ref'
	The second sample source will be available as a database named 'new'.

In the returned data.frame the there will be the following columns:
	'ref_sample_source' -> id of ref sample_source
	'new_sample_source' -> id of new sample_source
", sep=""))
	}

	ref_ss <- sample_sources[1,]
	con <- dbConnect(engine)
	set_db_cache_size(con, cache_size);
	dbGetQuery(con, paste("ATTACH DATABASE '", ref_ss$fname, "' AS ref;", sep=""))

	features <- ddply(sample_sources[seq(2,nrow(sample_sources)),], c("sample_source"), function(ss){
		tryCatch(c(ss),error=function(e){
			cat("ERROR: The specified sample source is not defined.\n")
		})
		cat("loading: ref:", as.character(ref_ss$sample_source),
				" new:", as.character(ss$sample_source), " ... ", sep="")
		if( is.na(ss$sample_source[1]) ){
			stop("Specified sample source is not defined")
		}
		timing <- system.time({
			dbGetQuery(con, paste("ATTACH DATABASE '", ss$fname, "' AS new;", sep=""))

			#Allow select statements to be prefaced with arbitrary statements.
			#This allows the creation of temporary tables, indices, etc.
			sele_split <- paste(strsplit(sele, ";\\W*", perl=TRUE)[[1]], ";", sep="")
			l_ply(sele_split[-length(sele_split)], function(sele){
				dbGetQuery(con, sele)
			})
			last_stmt <- sele_split[length(sele_split)] #Real statement for data
      if(! is.null(sele_args_frame)){
        df <- dbGetPreparedQuery(con, last_stmt, bind.data=sele_args_frame)
      }
      else{
			  df <- dbGetQuery(con, last_stmt)
      }
		})
		dbGetQuery(con, "DETACH DATABASE new;")
		cat(as.character(round(timing[3],2)),"s\n")
		df
	})
	dbDisconnect(con)

	if(nrow(features)==0){
		cat("WARNING: The following query returned no rows:\n")
		cat(sele)
    return(features)
	}
  if(char_as_factor){
	  for(col in names(features)){
		  if(is.character(features[,col])){
			  features[,col] <- factor(features[,col])
		  }
	  }
	}
	data.frame(
		ref_sample_source = factor(ref_ss$sample_source[1]),
		new_sample_source = factor(features$sample_source),
		sample_source = factor(features$sample_source),
		subset(features, select= -sample_source))
}

# Add a column to the data.frame called "counts" that for the total
# number of rows in each group, where the groups are determined by
# having the same values in the id.vars columns
group_counts <- function(f, id.vars) {
	ddply(f, id.vars, function(df) {
		df$counts = nrow(df)
		df
	})
}

# Add a column to the data.frame called "mean" that has the mean value
# of the measure.var, where the groups are determined by having the
# same values in the id.vars columns
group_means <- function(f, id.vars, measure.var, precision=4) {
	ddply(f, id.vars, transform, mean = round(mean(names(f)[1])))
}


locate_rosetta_application <- function(
	app_name,
	rosetta_base_path = NULL,
	platform=NULL,
	extras="default",
	compiler="clang",
	mode="release") {

	# find application is located at
	# ${rosetta_base_dir}/rosetta_source/bin/${app_name}.${extras}.${platform}${compiler}${mode}
	# detect what is not specified

	if(is.null(rosetta_base_path)){
		rosetta_base_path <- file.path(base_dir, "..", "..")
	}

	if(is.null(platform)){
		sysname <- Sys.info()[1]
		if (sysname == "Linux") {
			platform = "linux"
		} else if (sysname == "Darwin") {
			platform = "mac"
		} else {
			stop(paste("Unable to determine platform in trying to locate the application '", app_name, ". Try specifying the platform explicitly. (e.g. platform='linux' or platform='mac')", sep=""))
		}
	}

	full_app_path = file.path(rosetta_base_path, "rosetta_source", "bin", paste(app_name, ".", extras, ".", platform, compiler, mode, sep=""))
	if(!file.exists(full_app_path)){
		stop(paste("Looking for application '", app_name, "' at the path '", full_app_path, "', but it does not exist.", sep=""))
	}
	full_app_path
}
