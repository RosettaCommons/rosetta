# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
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

date_code <- function(d=NA){
  # reference http://www.r-cookbook.com/node/17
  if(is.na(d)) d <- Sys.Date()
  pattern <- '20([[:digit:]]{2})-([[:digit:]]{2})-([[:digit:]]{2})'
  paste(sub(pattern, '\\1', d),
        sub(pattern, '\\2', d),
        sub(pattern, '\\3', d),
        sep="")
}

check_setup <- function(){
  tryCatch(sample_sources, error=function(e){
    stop("ERROR: The variable 'sample_sources' is not defined. See compare_sample_sources.R")
  })
	if(!is.data.frame(sample_sources)){
		stop("ERROR: The variable 'sample_sources' is not a data frame.")
	}
	if(nrow(sample_sources) < 1){
		stop("ERROR: The variable 'sample_sources' contains no sample sources.")
	}

	tryCatch(output_dir, error=function(e){
    stop("ERROR: The variable 'output_dir' is not defined. See compare_sample_sources.R")
  })
	tryCatch(output_formats, error=function(e){
		stop("ERROR: The variable 'output_formats' is not defined. See compare_sample_sources.R")
	})
	tryCatch(db_cache_size, error=function(e){
		stop("ERROR: The variable 'db_cache_size' is not defined. See compare_sample_sources.R")
	})
}

save_plots <- function(
  plot_id,
  sample_sources,
  output_dir,
  output_formats,
  ...
) {
  a_ply(output_formats, 1, function(fmt){
    if(!file.exists(file.path(output_dir, fmt$id))){
      dir.create(file.path(output_dir, fmt$id), recursive=TRUE)
    }
    ss_ids <- paste(sample_sources$sample_source,collapse="_")
    fname <- paste(plot_id, date_code(), "with", ss_ids, sep="_")
    full_path <- file.path(output_dir, fmt$id, paste(fname, fmt$extension, sep=""))
    cat("Saving Plot: ", full_path, "\n")
    ggsave(
      filename=full_path,
      width=fmt$width,
      height=fmt$height,
      dpi=fmt$dpi,
      scale=fmt$scale,
      ...)
  })
}

query_sample_sources <- function(
  sample_sources,
  sele,
  cache_size=db_cache_size){
  tryCatch(sele,error=function(e){
    cat("ERROR: The select statement ", sele, " is not defined.\n")
  })
  features <- ddply(sample_sources, c("sample_source"), function(ss){

    tryCatch(c(ss),error=function(e){
      cat("ERROR: The specified sample source is not defined.\n")
    })
    cat("loading:", as.character(ss$sample_source), "... ")
    if( is.na(ss$sample_source[1]) ){
      stop("Specified sample source is not defined")
    }
    con <- dbConnect(engine, as.character(ss$fname))

    set_db_cache_size(con, cache_size);

		timing <- system.time({
	    #Allow select statements to be prefaced with arbitrary statements.
	    #This allows the creation of temporary tables, indices, etc.
	    sele_split <- paste(strsplit(sele, ";\\W*", perl=TRUE)[[1]], ";", sep="")
	    l_ply(sele_split[-length(sele_split)], function(sele){
	#      report_query_plan(con,sele);
	      dbSendQuery(con, sele)
	    })
	    last_stmt <- sele_split[length(sele_split)]
	#    report_query_plan(con, last_stmt)
	    df <- dbGetQuery(con, last_stmt)
		})
		cat(as.character(timing[3]),"s\n")
    df
  })
	for(col in names(features)){
	  if(is.character(features[,col])){
		  features[,col] <- factor(features[,col])
    }
  }
	features
}

set_db_cache_size <- function(con, cache_size){
  res <- dbSendQuery(con,
    paste("PRAGMA cache_size=",as.integer(cache_size),";",sep=""))
  dbClearResult(res)
}

n_pts <- 200
compute_density <- function(pts, wts){
  density(x=pts, from=xlim[1], to=xlim[2], n=n_pts, weights=wts)
}

estimate_density_1d <-function(
  data,
  ids,
  variable,
  weight_fun=uniform_normalization,
  min_count=20,
  n_pts=200,
  histogram=FALSE,
  ...){
	if(!(class(data) == "data.frame")){
		stop(paste("The data argument must be a data.frame, instead it is of class '", class(data), "'"))
	}
	for(id in ids){
		if(!(id %in% names(data))){
			stop(paste("The id variable '", id, "' is not a column name of the data. The ids are used to group the data instances for computing the density estimation.", sep=""))
		}
	}
	if(!(variable %in% names(data))){
		stop(paste("The value variable '", variable, "' is not a column name of the data. The value variable is used to compute the density estimation.", sep=""))
	}

  xlim <- range(data[,variable])
  compute_density <- function(factor_df){
    if (nrow(factor_df) < min_count){
      return( data.frame(x=seq(xlim[1], xlim[2], n_pts), y=0))
    } else {
      weights <- weight_fun(factor_df[,variable])
      if(histogram){
        breaks = seq(from=xlim[1], to=xlim[2], length=n_pts)
        d <- weighted.hist(x=factor_df[,variable], w=weights, plot=FALSE)
        return(data.frame(x=d$mids, y=d$density, counts=nrow(factor_df)))
      } else {
        d <- density(x=factor_df[,variable], from=xlim[1], to=xlim[2], n=n_pts,
          weights=weights,
          ...)
          return(data.frame(x=d$x, y=d$y, counts=nrow(factor_df)))
      }

    }
  }
  ddply(data, ids, compute_density)
}

estimate_density_1d_wrap <-function(
  data,
  ids,
  variable,
  weight_fun=uniform_normalization,
  min_count=20,
  n_pts=200,
	xlim=c(0,360),
  ...){
	if(!(class(data) == "data.frame")){
		stop(paste("The data argument must be a data.frame, instead it is of class '", class(data), "'"))
	}
	for(id in ids){
		if(!(id %in% names(data))){
			stop(paste("The id variable '", id, "' is not a column name of the data. The ids are used to group the data instances for computing the density estimation.", sep=""))
		}
	}
	if(!(variable %in% names(data))){
		stop(paste("The value variable '", variable, "' is not a column name of the data. The value variable is used to compute the density estimation.", sep=""))
	}

	extended_data <- mdply(
		c(xlim[1]-xlim[2], 0, xlim[2]-xlim[1]),function(s){
		y <- data;
		y[,variable] <- y[,variable] + s;
		y
	})
	extended_n_pts=n_pts*3
	extended_min_count=min_count*3
  compute_density <- function(factor_df){
    if (nrow(factor_df) < extended_min_count){
      return( data.frame(x=seq(xlim[1], xlim[2], length.out=n_pts), y=0))
    } else {
      weights <- weight_fun(factor_df[,variable])
			d <- density(x=factor_df[,variable], from=xlim[1], to=xlim[2], n=extended_n_pts,
				weights=weights,
        ...)
      return(data.frame(
				x=d$x[xlim[1] <= d$x & d$x <= xlim[2]],
				y=d$y[xlim[1] <= d$x & d$x <= xlim[2]]*3,
				counts=nrow(factor_df)/3))
    }
  }
  ddply(extended_data, ids, compute_density)
}


estimate_density_1d_reflect_boundary <-function(
  data,
  ids,
  variable,
  weight_fun=uniform_normalization,
  min_count=20,
  n_pts=200,
	reflect_left=FALSE,
	reflect_right=FALSE,
  left_boundary=NULL,
	right_boundary=NULL,
  ...){
	if(!(class(data) == "data.frame")){
		stop(paste("The data argument must be a data.frame, instead it is of class '", class(data), "'"))
	}
	for(id in ids){
		if(!(id %in% names(data))){
			stop(paste("The id variable '", id, "' is not a column name of the data. The ids are used to group the data instances for computing the density estimation.", sep=""))
		}
	}
	if(!(variable %in% names(data))){
		stop(paste("The value variable '", variable, "' is not a column name of the data. The value variable is used to compute the density estimation.", sep=""))
	}

	if(is.null(left_boundary)){
		left_boundary=min(data[,variable])
	}
	if(is.null(right_boundary)){
		right_boundary=max(data[,variable])
	}

	extended_factor = 1

	if(reflect_left==TRUE){
		data_lower <- data
		data_lower[,variable] <- 2*left_boundary - data[,variable]
		data <- rbind(data, data_lower)
		extended_factor = extended_factor + 1
	}

	if(reflect_right==TRUE){
		data_upper <- data
		data_upper[,variable] <- 2*right_boundary - data[,variable]
		data <- rbind(data, data_upper)
		extended_factor = extended_factor + 1
	}
  compute_density <- function(factor_df){
			if (nrow(factor_df) < min_count*extended_factor){
				return(data.frame(x=seq(left_boundary, right_boundary, length.out=n_pts), y=0, counts=nrow(factor_df)))
    } else {
      weights <- weight_fun(factor_df[,variable])
			d <- density(x=factor_df[,variable], from=left_boundary, to=right_boundary, n=n_pts*extended_factor,
				weights=weights,
        ...)

			return(data.frame(
				x=d$x[left_boundary <= d$x & d$x <= right_boundary],
				y=d$y[left_boundary <= d$x & d$x <= right_boundary]*extended_factor,
				counts=nrow(factor_df)/extended_factor))
    }
  }
	ddply(data, ids, compute_density)
}


estimate_density_1d_logspline <-function(
  data,
  ids,
  variable,
  min_count=20,
  n_pts=200,
  weight_fun=NULL,
  ...){

	if(!(class(data) == "data.frame")){
		stop(paste("The data argument must be a data.frame, instead it is of class '", class(data), "'", sep=""))
	}
	for(id in ids){
		if(!(id %in% names(data))){
			stop(paste("The id variable '", id, "' is not a column name of the data. The ids are used to group the data instances for computing the density estimation.", sep=""))
		}
	}
	if(!(variable %in% names(data))){
		stop(paste("The value variable '", variable, "' is not a column name of the data. The value variable is used to compute the density estimation.", sep=""))
	}

  xlim <- range(data[,variable])
  if(!is.null(weight_fun)){
    xlim_transformed <- weight_fun(xlim)
  }

  compute_density <- function(factor_df){
    if (nrow(factor_df) < min_count){
      return( data.frame(x=seq(xlim[1], xlim[2], n_pts), y=0))
    } else {
      if(!is.null(weight_fun)){
        values_transformed <- weight_fun(factor_df[,variable])
        lgs <- logspline(
          values_transformed,
          lbound=xlim_transformed[1],
          ubound=xlim_transformed[2])
        x_transformed <- seq(
          from=xlim_transformed[1],
          to=xlim_transformed[2],
          length.out=n_pts)
        y <- dlogspline(x_transformed, lgs)
        x <- seq(from=xlim[1], to=xlim[2], length.out=n_pts)

      } else {
        lgs <- logspline(factor_df[,variable], lbound=xlim[1], ubound=xlim[2])
        x <- seq(from=xlim[1], to=xlim[2], length.out=n_pts)
        y <- dlogspline(x, lgs)
      }
      return(data.frame(x=x, y=y, counts=nrow(factor_df)))
    }
  }
  ddply(data, ids, compute_density)
}


estimate_density_2d <-function(
  data,
  ids,
  xvariable,
	yvariable,
  min_count=20,
  n_pts=200,
  histogram=FALSE,
  ...){
	if(!(class(data) == "data.frame")){
		stop(paste("The data argument must be a data.frame, instead it is of class '", class(data), "'", sep=""))
	}
	for(id in ids){
		if(!(id %in% names(data))){
			stop(paste("The id variable '", id, "' is not a column name of the data. The ids are used to group the data instances for computing the density estimation.", sep=""))
		}
	}
	if(!(xvariable %in% names(data))){
		stop(paste("The value variable '", xvariable, "' is not a column name of the data. The xvariable and yvariable are used to compute the density estimation.", sep=""))
	}
	if(!(yvariable %in% names(data))){
		stop(paste("The value variable '", yvariable, "' is not a column name of the data. The xvariable and yvariable are used to compute the density estimation.", sep=""))
	}

	xlim <- range(data[,xvariable])
  xlim <- range(data[,yvariable])
	ddply(data, ids, function(df){
	  if (nrow(df) < min_count){
      d <- data.frame(x=NULL, y=NULL, z=NULL)
    } else {
			if(histogram){
			  h <- gplots::hist2d(
          x=as.matrix(df[,c(xvariable, yvariable)]), nbins=n_pts, show=FALSE)
        d <- with(h, data.frame(expand.grid(x=x, y=y), z=as.vector(counts), density=as.vector(counts)/nrow(df)))
        d$z <- melt(h$counts)$value
			} else {
				dm <- MASS::kde2d(
				  x=df[,xvariable], y=df[,yvariable], n_pts)
				d <- with(dm, data.frame(expand.grid(x=x, y=y), z=as.vector(z)))
      }
    }
    d
  })
}
