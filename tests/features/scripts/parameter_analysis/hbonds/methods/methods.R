# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


require(polynom)

valid_rosetta_database_path <- function(rosetta_database_path) {
	if(is.null(rosetta_database_path)){
		rosetta_database_path <- Sys.getenv("ROSETTA3_DB")
	}
	if(is.null(rosetta_database_path) || !file.exists(rosetta_database_path)){
		stop("Please specify the rosetta_database path via the --database command line argument or set the $ROSETTA3_DB environment variable.")
	}
	rosetta_database_path
}



#It's necessary to actually go into the folder because schema.sql
#tries to import the .csv files with in the current directory.
convert_parameter_set_into_database <- function(
	rosetta_database_path, parameter_set, debug=FALSE){

	parameter_set_path <- paste(rosetta_database_path, "scoring", "score_functions", "hbonds", parameter_set, sep="/")
	if(!file.exists(parameter_set_path)){
		stop(paste("Attempting to convert the hbond parameter set '",
			parameter_set_path, "' ",
			"set into a database, but it does not exist.", sep=""))
	}
	cwd <- getwd()
	setwd(parameter_set_path)
	system("sqlite3 params.db3 < schema.sql")
	setwd(cwd)
	parameter_set_path
}


#################  Polynomial Functions ###################

# return the HBPoly1D table
get_polynomials <- function(rosetta_database_path, parameter_set, debug=FALSE){
	con <- dbConnect("SQLite",
		paste(rosetta_database_path,
			"scoring/score_functions/hbonds", parameter_set, "params.db3", sep="/"))
	poly_params <- dbGetQuery(con, "SELECT * FROM HBPoly1D;")
	dbDisconnect(con)

	if(debug){
		cat("poly_params:\n")
		print(poly_params)
	}
	poly_params
}

#return a polynom::polynomial object that represents the polynomial
#with the specified name in the specified polynomials data.frame.
get_1d_polynomial <- function(polynomials, name){
	df <- polynomials[polynomials$name==name,]
	if (nrow(df) !=1){
		stop(paste("Polynomial, ", name,", could not be found.", sep=""))
	}
	#we don't know the order in which the rows columns selected will be returned!
	switch(df$degree,
		p <- df[1,c("c_a")],
		p <- df[1,c("c_b","c_a")],
		p <- df[1,c("c_c","c_b","c_a")],
		p <- df[1,c("c_d","c_c","c_b","c_a",)],
		p <- df[1,c("c_e","c_d","c_c","c_b","c_a")],
		p <- df[1,c("c_f","c_e","c_d","c_c","c_b","c_a")],
		p <- df[1,c("c_g","c_f","c_e","c_d","c_c","c_b","c_a")],
		p <- df[1,c("c_h","c_g","c_f","c_e","c_d","c_c","c_b","c_a")],
		p <- df[1,c("c_i","c_h","c_g","c_f","c_e","c_d","c_c","c_b","c_a")],
		p <- df[1,c("c_j","c_i","c_h","c_g","c_f","c_e","c_d","c_c","c_b","c_a")],
		p <- df[1,c("c_k","c_j","c_i","c_h","c_g","c_f","c_e","c_d","c_c","c_b","c_a")])

	polynomial(p)


}

get_1d_polynomial_from_types <- function(polynomials, sample_source, acc_chem_type, don_chem_type, separation, dimension){
	df <- polynomials[
		polynomials$sample_source == as.character(sample_source) &
		polynomials$don_chem_type == as.character(don_chem_type) &
		polynomials$acc_chem_type == as.character(acc_chem_type) &
		polynomials$separation == as.character(separation) &
		polynomials$dimension == as.character(dimension),]
	if(nrow(df) != 1){
		print("ERROR: in get_1d_polynomial_from_types")
		print(paste(sample_source, don_chem_type, acc_chem_type, separation, dimension))
		print(polynomials)
		print(df)
	}
	switch(df$degree,
		p <- df[1,c("c_a")],
		p <- df[1,c("c_b","c_a")],
		p <- df[1,c("c_c","c_b","c_a")],
		p <- df[1,c("c_d","c_c","c_b","c_a",)],
		p <- df[1,c("c_e","c_d","c_c","c_b","c_a")],
		p <- df[1,c("c_f","c_e","c_d","c_c","c_b","c_a")],
		p <- df[1,c("c_g","c_f","c_e","c_d","c_c","c_b","c_a")],
		p <- df[1,c("c_h","c_g","c_f","c_e","c_d","c_c","c_b","c_a")],
		p <- df[1,c("c_i","c_h","c_g","c_f","c_e","c_d","c_c","c_b","c_a")],
		p <- df[1,c("c_j","c_i","c_h","c_g","c_f","c_e","c_d","c_c","c_b","c_a")],
		p <- df[1,c("c_k","c_j","c_i","c_h","c_g","c_f","c_e","c_d","c_c","c_b","c_a")])
	polynomial(p)
}



polynomial_range <- function(polynomials, name){
	if(!(name %in% polynomials$name)){
		stop(paste("Polynomial name, ", name,", is not recognized.", sep=""))
	}
	c(polynomials[polynomials$name==name, "xmin"], polynomials[polynomials$name==name, "xmax"])
}

display_polynomials <- function(parameter_set_path, polynomials){
	polynomials <- adply(polynomials, 1, function(polynomial){
		name <- polynomial[1, "name"]
		polynomial$equation <- as.character(signif(get_1d_polynomial(polynomials,name), digits=4), decreasing=TRUE)
		polynomial
	})
	cat("Plynomials defined in the parameter set '", parameter_set_path,"'\n", sep="")
	print(polynomials[,c("name", "equation")], right=FALSE)
}

evaluate_polynomial <- function(params, name){
	print(params)
	print(class(params))

	poly <- params[params$name == name,]
	x <- seq(poly$xmin, poly$xmax, length.out=100)
	p <- get_1d_polynomial(params, name)
	data.frame(name=name, x=x, y=aaply(x, 1, as.function(p)))
}

#################  Fade Interval Functions ###################

# return the HBFadeIntervals table
get_fade_intervals <- function(parameter_set_path, debug=FALSE){
	con <- dbConnect("SQLite", paste(parameter_set_path, "/params.db3", sep=""))
	fade_params <- dbGetQuery(con, "SELECT * FROM HBFadeIntervals;")
	dbDisconnect(con)
	
	if(debug){
		cat("fade_params:\n")
		print(fade_params)
	}
	fade_params
}

# the defining points of the piecewise linear for the fade interval
fade_interval_knots <- function(fade_intervals, name, lim){
	df <- fade_intervals[fade_intervals$name==name,]
	if (nrow(df) != 1){
		print(fade_intervals)
		stop(paste("Fade interval, ", name, ", could not be found.", sep=""))
	}

	if(name=="fade_zero"){
		return(data.frame(
			name=name,
			x=c(lim[1],(3*lim[1]+lim[2])/4,(lim[1]+3*lim[2])/4,lim[2]),
			y=c(0,0,0,0)))
	}

	data.frame(
		name=name,
		x=c(df[1,"min0"], df[1,"fmin"], df[1,"fmax"], df[1,"max0"]),
		y=c(0,1,1,0))
}


# return function that evaluates the fade value
# see core::scoring::hbonds::FadeInterval::value()
fade_function <- function(fade_interval_knots, x, smooth=FALSE){
	knots <- fade_interval_knots[,"x"]
	y = NA
	if(fade_interval_knots[1,"name"]=="fade_zero") return(0)
	if (x <= knots[3]){
		if(x <= knots[1]){
			y = 0
		} else if(x >= knots[2]) {
			y = 1
		} else {
			z= (x - knots[1])/(knots[2]-knots[1])
			if(smooth){
				y = -2*z^3 + 3*z^2
			} else {
				y = z
			}
		}
	} else {
		if(x >= knots[4]){
			y = 0
		} else {
			z = (x-knots[3])/(knots[4]-knots[3])
			if(smooth){
				y = 2*z^3 - 3*z^2 + 1
			} else {
				y = 1-z
			}
		}
	}
#	cat("fade_function: (", knots[1], " ", knots[2], " ", knots[3], " ", knots[4], ") x=", x, " y=", y, "\n")
	y
}

evaluate_fade_interval <- function(params, xmin, xmax, name){
  x <- seq(xmin, xmax, length.out=100)
	cat("evaluating fade interval for ", name,"\n")
	print(params)
	p <- data.matrix(params[,c("smooth", "min0", "fmin", "fmax", "max0")])
  data.frame(name=name, x=x, y=aaply(x, 1, eval_fade, p))
}

#################  HBond Parameter Functions  ###################
get_polynomial_parameters <- function(
	parameter_set_path, don_chem_type, acc_chem_type, separation, debug=FALSE){

	con <- dbConnect("SQLite", paste(parameter_set_path, "/params.db3", sep=""))
	poly_params <- function(which_poly){
		p <- dbGetQuery(con, paste("
SELECT
	p.degree, p.xmin, p.xmax,
  p.c_a, p.c_b, p.c_c, p.c_d, p.c_e, p.c_f, p.c_g, p.c_h, p.c_j, p.c_k
FROM
  HBPoly1D as p,
  HBEval as e
WHERE
  e.don_chem_type='", don_chem_type, "' AND
  e.acc_chem_type='", acc_chem_type, "' AND
	e.separation='", separation, "' AND
  e.", which_poly, "= p.name;", sep=""))
		p$name <- which_poly
		p
	}

	params <- rbind(
		poly_params("AHdist"),
		poly_params("cosBAH_short"),
		poly_params("cosBAH_long"),
		poly_params("cosAHD_short"),
		poly_params("cosAHD_long"))
	dbDisconnect(con)
	params
}

get_fade_parameters <- function(
	parameter_set_path, don_chem_type, acc_chem_type, separation, debug=FALSE){

	con <- dbConnect("SQLite", paste(parameter_set_path, "/params.db3", sep=""))
	fade_params <- function(which_fade){
		f <- dbGetQuery(con, paste("
SELECT
  CASE f.junction_type WHEN 'smooth' THEN 1 ELSE 0 END AS smooth,
  f.min0, f.fmin, f.fmax, f.max0
FROM
  HBFadeIntervals as f,
  HBEval as e
WHERE
  e.don_chem_type='", don_chem_type, "' AND
  e.acc_chem_type='", acc_chem_type, "' AND
	e.separation='", separation, "' AND
  e.", which_fade, "= f.name;", sep=""))
		f$which_fade = which_fade
		f
	}
	params <- rbind(
		fade_params("AHdist_short_fade"),
		fade_params("AHdist_long_fade"),
		fade_params("cosBAH_fade"),
		fade_params("cosAHD_fade"))
	dbDisconnect(con)
	params
}



# the donor chemical type, acceptor chemical type and sequence
# separation uniquely identify hydrogen bond evaluation classes.
# return a row from HBEval for the specified class.
get_evaluation_type <- function(
	parameter_set_path, don_chem_type, acc_chem_type, separation){

	con <- dbConnect("SQLite", paste(parameter_set_path, "/params.db3", sep=""))
	eval_type <- dbGetQuery(con, paste("
SELECT * FROM HBEval
WHERE
	don_chem_type='", don_chem_type, "' AND
	acc_chem_type='", acc_chem_type, "' AND
	separation='", separation, "';", sep=""))

	# now do error checking...
	if(nrow(eval_type) != 1){
		# check which part of the key is causing problems
		don_row <- dbGetQuery(con, paste("
SELECT * FROM HBDonChemType
WHERE don_chem_type='",don_chem_type,"';", sep=""))
		if(nrow(don_row) != 1){
			stop(paste("The donor chemical type,",
				don_chem_type, ", is not recognized."))
		}

		acc_row <- dbGetQuery(con, paste("
SELECT * FROM HBAccChemType
WHERE acc_chem_type='",acc_chem_type,"';", sep=""))
		if(nrow(acc_row) != 1){
			stop(paste("The acceptor chemical type,",
				acc_chem_type, ", is not recognized."))
		}

		sep_row <- dbGetQuery(con, paste("
SELECT * FROM HBSeqSep
WHERE name='",separation,"';", sep=""))
		if(nrow(sep_row) != 1){
			stop(paste("The sequence separation type,",
				separation, ", is not recognized."))
		}

		stop(paste("Unable to find evaluation type for:
			don_chem_type: '", don_chem_type, "'
			acc_chem_type: '", acc_chem_type, "'
			separation:    '", separation, "'", sep=""))
	}
	dbDisconnect(con)
	eval_type
}

get_evaluation_function <- function(
	polynomials, fade_intervals, eval_type){

	AHdist_poly <- get_1d_polynomial(polynomials, eval_type[1,"AHdist"])
	cosBAH_short_poly <- get_1d_polynomial(polynomials, eval_type[1,"cosBAH_short"])
	cosBAH_long_poly <- get_1d_polynomial(polynomials, eval_type[1,"cosBAH_long"])
	cosAHD_short_poly <- get_1d_polynomial(polynomials, eval_type[1,"cosAHD_short"])
	cosAHD_long_poly <- get_1d_polynomial(polynomials, eval_type[1,"cosAHD_long"])

	AHdist_range = polynomial_range(polynomials, eval_type[1, "AHdist"])
	cosBAH_range = polynomial_range(polynomials, eval_type[1, "cosBAH_short"])
	cosAHD_range = polynomial_range(polynomials, eval_type[1, "cosAHD_short"])

	AHdist_short_fade <- fade_interval_knots(
		fade_intervals, eval_type[1,"AHdist_short_fade"], AHdist_range)
	AHdist_long_fade <- fade_interval_knots(
		fade_intervals, eval_type[1,"AHdist_long_fade"], AHdist_range)
	cosBAH_fade <- fade_interval_knots(
		fade_intervals, eval_type[1, "cosBAH_fade"], cosBAH_range)
	cosAHD_fade <- fade_interval_knots(
		fade_intervals, eval_type[1, "cosAHD_fade"], cosAHD_range)

	function(AHdist, cosBAH, cosAHD){
		result = 0
		if(AHdist > AHdist_range[1] && AHdist < AHdist_range[2]){
			result = result +
				as.function(AHdist_poly)(AHdist) *
					fade_function(cosBAH_fade, cosBAH, opt$smooth) *
						fade_function(cosAHD_fade, cosAHD, opt$smooth)
		}
		if(cosBAH > cosBAH_range[1] && cosBAH < cosBAH_range[2]){
			result = result +
				fade_function(AHdist_short_fade, AHdist, opt$smooth) *
					as.function(cosBAH_short_poly)(cosBAH) *
						fade_function(cosAHD_fade, cosAHD, opt$smooth)
			result = result +
				fade_function(AHdist_long_fade, AHdist, opt$smooth) *
					as.function(cosBAH_long_poly)(cosBAH) *
						fade_function(cosAHD_fade, cosAHD, opt$smooth)
		}
		if(cosAHD > cosAHD_range[1] && cosAHD < cosAHD_range[2]){
			result = result +
				fade_function(AHdist_short_fade, AHdist, opt$smooth) *
					fade_function(cosBAH_fade, cosBAH, opt$smooth) *
						as.function(cosAHD_short_poly)(cosAHD)
			result = result +
				fade_function(AHdist_long_fade, AHdist, opt$smooth) *
					fade_function(cosBAH_fade, cosBAH, opt$smooth) *
						as.function(cosAHD_long_poly)(cosAHD)
		}
		result
	}
}


evaluate_function <- function(evaluation_function, AHdist, cosBAH, cosAHD){
	ddply(expand.grid(list(AHdist=AHdist, cosBAH=cosBAH, cosAHD=cosAHD)),
		.variables=c("AHdist", "cosBAH", "cosAHD"), function(df){
		data.frame(energy=evaluation_function(
			df[1,"AHdist"], df[1,"cosBAH"], df[1, "cosAHD"]))
	})
}


