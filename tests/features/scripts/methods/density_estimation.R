# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.



estimate_density_1d <-function(
  data,
  ids,
  variable,
  weight_fun=uniform_normalization,
  min_count=10,
  n_pts=512,
  histogram=FALSE,
	sample_domain=NULL,
	adjust=1,
  ...){
	density.args <- list(...)
	data <- as.data.frame(data)
	if(!(class(data) == "data.frame")){
		stop(paste("The data argument must be a data.frame, instead it is of class '", class(data), "'"))
	}
	if(nrow(data) == 0){
		stop(paste("Unable to compute density estimation because the data argument has no rows."))
	}
	for(id in ids){
		if(!(id %in% names(data))){
			stop(paste("The id variable '", id, "' is not a column name of the data. The ids are used to group the data instances for computing the density estimation.", sep=""))
		}
	}
	if(!(variable %in% names(data))){
		stop(paste("The value variable '", variable, "' is not a column name of the data. The value variable is used to compute the density estimation.", sep=""))
	}

	if(is.null(sample_domain)){
		sample_domain <- range(data[,variable])
	}
  compute_density <- function(factor_df){
    if (nrow(factor_df) < min_count){
			#mjo I was having issues with whole columns/rows
			#disappearing. This worked for ggplot2 version 0.8.9 but now
			#gives warnings: asking "Do you need to adjust the group
			#aesthetic?" with version 0.9.0. I'm trying just returning an
			#empty data.frame to see if that works alright.
			#return( data.frame(x=seq(sample_domain[1], sample_domain[2], n_pts), y=0))

			return(data.frame())

    } else {
      weights <- weight_fun(factor_df[,variable])
      if(histogram){
        breaks = seq(from=sample_domain[1], to=sample_domain[2], length=n_pts)
        d <- weighted.hist(x=factor_df[,variable], w=weights, breaks=breaks, plot=FALSE)
        return(data.frame(x=d$mids, y=d$density, counts=nrow(factor_df)))
      } else {
        adjust <- adjust * general_kernel_adjust
        d <- do.call(density,
					c(list(x=factor_df[,variable], from=sample_domain[1], to=sample_domain[2], n=n_pts,
          weights=weights, adjust=adjust), density.args))
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
  n_pts=512,
	xlim=c(0,360),
	adjust=1,
  ...){
	density.args <- list(...)
	if(!(class(data) == "data.frame")){
		stop(paste("The data argument must be a data.frame, instead it is of class '", class(data), "'"))
	}
	if(nrow(data) == 0){
		stop(paste("Unable to compute density estimation because the data argument has no rows."))
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

	extended_n_pts <- n_pts*3
	extended_min_count <- min_count*3

	# there is three times the data and assume the the bin width scales linearly with the sample size
	adjust <- adjust/3
  compute_density <- function(factor_df){
    if (nrow(factor_df) < extended_min_count){
      return( data.frame(x=seq(xlim[1], xlim[2], length.out=n_pts), y=0))
    } else {
      weights <- weight_fun(factor_df[,variable])
			d <- do.call(density,
				c(list(x=factor_df[,variable], from=xlim[1], to=xlim[2], n=extended_n_pts,
					weights=weights, adjust=adjust), density.args))
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
  n_pts=512,
	reflect_left=FALSE,
	reflect_right=FALSE,
  left_boundary=NULL,
	right_boundary=NULL,
	adjust=1,
  ...){
	if(!(class(data) == "data.frame")){
		stop(paste("The data argument must be a data.frame, instead it is of class '", class(data), "'"))
	}
	if(nrow(data) == 0){
		stop(paste("Unable to compute density estimation because the data argument has no rows."))
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
	extended_factor <- 1

	# Do the density estimation on the covering variable but apply the
	# weight funtion to the original variable
	data$covering_variable <- data[,variable]

	if(reflect_left==TRUE){
		data_lower <- data
		data_lower$covering_variable <- 2*left_boundary - data[,variable]
		data <- rbind(data, data_lower)
		extended_factor = extended_factor + 1
	}

	if(reflect_right==TRUE){
		data_upper <- data
		data_upper$covering_variable <- 2*right_boundary - data[,variable]
		data <- rbind(data, data_upper)
		extended_factor = extended_factor + 1
	}

	adjust <- adjust * general_kernel_adjust
	adjust <- adjust/(1 + reflect_left + reflect_right)

  compute_density <- function(factor_df){
			if (nrow(factor_df) < min_count*extended_factor){
				return(data.frame(x=seq(left_boundary, right_boundary, length.out=n_pts), y=0, counts=nrow(factor_df)))
    } else {
      weights <- weight_fun(factor_df[,variable])
			d <- density(
				x=factor_df$covering_variable,
				adjust=adjust,
				weights=weights,
				n=n_pts*extended_factor,
				from=left_boundary,
				to=right_boundary,
        ...)

			return(data.frame(
				x=d$x[left_boundary <= d$x & d$x <= right_boundary],
				y=d$y[left_boundary <= d$x & d$x <= right_boundary]*extended_factor,
				counts=round(nrow(factor_df)/extended_factor),0))
    }
  }
	z <- ddply(data, ids, compute_density)
	z[,!(names(z) %in% "X0")]
}


estimate_density_1d_logspline <-function(
  data,
  ids,
  variable,
  min_count=20,
  n_pts=512,
  weight_fun=NULL,
  ...){

	if(!(class(data) == "data.frame")){
		stop(paste("The data argument must be a data.frame, instead it is of class '", class(data), "'", sep=""))
	}
	if(nrow(data) == 0){
		stop(paste("Unable to compute density estimation because the data argument has no rows."))
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
  n_pts=512,
  histogram=FALSE,
  ...){
	if(!(class(data) == "data.frame")){
		stop(paste("The data argument must be a data.frame, instead it is of class '", class(data), "'", sep=""))
	}
	if(nrow(data) == 0){
		stop(paste("Unable to compute density estimation because the data argument has no rows."))
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

# Compute quantiles for specified 'variable' grouping by 'ids' columns
# in 'variable'. If lengths(probs)==1, then use "Ordinates for
# Probability Plotting" with 'probs' probability points.
#
# Returns data.frame with the following columns
#   ids, probs, quantiles, counts
#
# # for example:
# q <- compute_quantiles(f, c("id_col1", "id_col2"), "var")
# ggplot(q) + geom_line(aes(x=quantiles*180/pi, y=probs))
compute_quantiles <- function(
	data,
	ids,
	variable,
	probs=ppoints(1000)
) {
	ddply(data, ids, function(df){
		data.frame(
			probs=probs, quantiles=quantile(df[,variable], probs=probs), counts=nrow(df))
	})
}
