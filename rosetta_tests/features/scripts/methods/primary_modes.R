# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


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

estimate_primary_modes_1d <- function(
	data,
	ids,
	variable,
	weight_fun=uniform_normalization,
	min_count=20,
	n_pts=1024,
	sample_domain=NULL,
	bandwidth_tolerance = .001,
	debug=FALSE,
	...
){
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
	locate_primary_mode <-function(factor_df){
		if (nrow(factor_df) < min_count){
			return(data.frame())
		}
		if(debug){
			print(factor_df[1, ids])
		}
		weights <- weight_fun(factor_df[,variable])
		t <- do.call(
			one_mode_test_fun,
			c(list(
				x=factor_df[,variable],
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
				x=factor_df[,variable],
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
		return(data.frame(
			min_kernel_adjust=min_kernel_adjust,
			mode_location=d$x[modes[1]],
			counts=nrow(factor_df)))
	}
	ddply(data, ids, locate_primary_mode)
}
