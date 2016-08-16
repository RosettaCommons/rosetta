# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


# provides estimate_primary_mode_1d, which estimates the location of
# 	the primary mode by binary searching the bandwidth of kernel
# 	density estimation till there is just one mode.

find_modes_in_density <- function(dens_y){
	which(diff(sign(diff(dens_y)))==-2)
}

one_mode_test_fun <- function(x, ...){
	function(bandwidth_adjust){
		d <- density(x=x, adjust=bandwidth_adjust, ...)
		modes <- find_modes_in_density(d$y)
		if( length(modes) == 1){
			#kernel adjustment is too big
			-1
		} else {
			#kernel adjustment is too small
			1
		}
	}
}

binary_search <- function(
	test_fun,
	min=0,
	max=100,
	tol=.001,
	debug=FALSE){

	step = (max - min)/2
	cur = min + step
	while(step*2 > tol){
		if(debug) cat("binary_search: cur: ", cur, "\tstep: ", step, "\n")
		step <- step/2
		if(test_fun(cur) > 0){
			cur <- cur + step
		} else {
			cur <- cur - step
		}
	}
	cur
}


locate_primary_mode <- function(
	x,
	weight_fun=uniform_normalization,
	sample_domain=range(x),
	n_pts=1024,
	bandwidth_tolerance = .001,
	debug=FALSE,
	...
) {
	density.args <- list(...)
	weights <- weight_fun(x)
	t <- do.call(
		one_mode_test_fun,
		c(list(
			x=x,
			from=sample_domain[1],
			to=sample_domain[2],
			n=n_pts,
			weights=weights),
			density.args))
	min_kernel_adjust <-
		binary_search(t, tol=bandwidth_tolerance, debug=debug) + bandwidth_tolerance
	d <- do.call(
		density,
		c(list(
			x=x,
			from=sample_domain[1],
			to=sample_domain[2],
			n=n_pts,
			weights=weights,
			min_kernel_adjust),
			density.args))
	modes <- find_modes_in_density(d$y)
	if( length(modes) != 1){
		stop(paste("Attempted to find minimum kernel that has just one mode, but somehow it eneded up with '", length(modes), "' modes.", sep=""))
	}
	if(debug){
		return(data.frame(
			min_kernel_adjust=min_kernel_adjust,
			statistic=d$x[modes[1]]))
	} else {
		return(c(primary_mode=d$x[modes[1]]))
	}
}


# For each set of rows grouped by the columns in ids, binary search
# for the minimum kernel bandwidth that gives a single mode with
# kernel density estimation.
#
# returns a data.frame with columns c(ids, "statistic") where
#   "statistic" is the value of "variable" where the mode occurs in
#   each group defined by ids.
#
# data: A data.frame with columns ids and variable
#
# ids: A vector of strings that are the column names of data that is
# 	used to group the data to make multiple comparisons at once. See
# 	the Plyr package for details.
#
# variable: A string for a numeric column of data where the primary
# 	mode is to be computed.
#
# min_count: The minimum number of rows in group before the primary
# 	mode is computed. This helps over-interpreting groups with very
# 	few counts
#
# n_pts: The resolution of "variable" in measuring the location of the mode.
#
# sample_domain: The limits of "variable", when computing the primary
# 	mode. If it isn't specified then take the whole range of
# 	"variable" over all the data.
#
# banddwidth_tolerance: Stop the binary search of kernel density estimation
#
# debug: Print out convergence information in computing the primary
# 	mode, and add the column "min_kernel_adjust" to the returned
# 	data.frame of the value of the kernel bandwidth for each group.
#
# Extra parameters are passed to the density function.
estimate_primary_modes_1d <- function(
	data,
	ids,
	variable,
	sample_domain=NULL,
	min_count=20,
	...
){
	extra.args <- list(...)
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

	ddply(data, ids, function(df){
		if (nrow(df) < min_count){
			return(data.frame())
		}
		args <- c(list(x=df[,variable], sample_domain=sample_domain), extra.args)
		data.frame(primary_mode = do.call(locate_primary_mode, args))
	})
}

primary_modes_diff <- function(a, b, ...){
	extra.args <- list(...)
	data.frame(
		statistic_name=factor("Primary Mode Diff"),
		statistic=abs(
			do.call(locate_primary_mode, c(list(a), extra.args)) -
			do.call(locate_primary_mode, c(list(b), extra.args))),
		p.value=NA)
}
