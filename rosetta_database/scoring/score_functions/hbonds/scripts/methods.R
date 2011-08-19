# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


require(polynom)

#It's necessary to actually go into the folder because schema.sql
#tries to import the .csv files with in the current directory.
convert_parameter_set_into_database <- function(parameter_set_path){
	cwd <- getwd()
	setwd(parameter_set_path)
	system("sqlite3 params.db3 < schema.sql")
	setwd(cwd)
}

#################  Polynomial Functions ###################

# return the HBPoly1D table
get_polynomials <- function(parameter_set_path){
	con <- dbConnect("SQLite", paste(parameter_set_path, "/params.db3", sep=""))
	dbGetQuery(con, "SELECT * FROM HBPoly1D;")
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

polynomial_range <- function(polynomials, name){
	if(!(name %in% polynomials$name)){
		stop(paste("Polynomial name, ", name,", is not recognized.", sep=""))
	}
	c(polynomials[polynomials$name==name, "xmin"], polynomials[polynomials$name==name, "xmax"])
}

display_polynomials <- function(parameter_set_path, polynomials){
	polynomials <- adply(polynomials, 1, function(polynomial){
		name <- polynomial[1, "name"]
		polynomial$equation <- as.character(signif(get_1d_polynomial_by_name(polynomials,name), digits=4), decreasing=TRUE)
		polynomial
	})
	cat("Plynomials defined in the parameter set ", parameter_set_path,"\n")
	print(polynomials[,c("name", "equation")], right=FALSE)
}

evaluate_polynomial <-function(polynomials, name, lenght.out=50){
	p <- get_1d_polynomial(polynomials, name)
	xmin <- polynomials[polynomials$name==name, "xmin"]
	xmax <- polynomials[polynomials$name==name, "xmax"]
	data.frame(
		name=name,
		x=seq(xmin, xmax, length.out=50),
		y=predict(p, seq(xmin, xmax, length.out=50)))
}

#################  Fade Interval Functions ###################

# return the HBFadeIntervals table
get_fade_intervals <- function(parameter_set_path){
	con <- dbConnect("SQLite", paste(parameter_set_path, "/params.db3", sep=""))
	dbGetQuery(con, "SELECT * FROM HBFadeIntervals;")
}

# the defining points of the piecewise linear for the fade interval
evaluate_fade_interval <- function(fade_intervals, name, lim){
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
fade_function <- function(fade_interval, x){
	if(fade_interval[1,"name"]=="fade_zero") return(0)
	if (x <= fade_interval[3,"x"]){
		if(x <= fade_interval[1,"x"]) return(0)
		if(x >= fade_interval[2,"x"]) return(1)
		return((x - fade_interval[1,"x"])/(fade_interval[2,"x"]-fade_interval[1,"x"]))
	} else {
		if(x >= fade_interval[4,"x"]) return(0)
		return((fade_interval[4, "x"] - x)/(fade_interval[4,"x"]-fade_interval[3,"x"]))
	}
}



#################  HBond Parameter Functions  ###################

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
	eval_type
}

get_evaluation_function <- function(
	polynomials, fade_intervals, eval_type){

	AHdist_poly <- get_1d_polynomial(polynomials, eval_type[1,"AHdist"])
	cosBAH_short_poly <- get_1d_polynomial(polynomials, eval_type[1,"cosBAH_short"])
	cosBAH_long_poly <- get_1d_polynomial(polynomials, eval_type[1,"cosBAH_long"])
	cosAHD_short_poly <- get_1d_polynomial(polynomials, eval_type[1,"cosAHD_short"])
	cosAHD_long_poly <- get_1d_polynomial(polynomials, eval_type[1,"cosAHD_long"])
	AHdist_short_fade <- evaluate_fade_interval(
		fade_intervals, eval_type[1,"AHdist_short_fade"], polynomial_range(polynomials, eval_type[1, "AHdist"]))
	AHdist_long_fade <- evaluate_fade_interval(
		fade_intervals, eval_type[1,"AHdist_long_fade"], polynomial_range(polynomials, eval_type[1, "AHdist"]))
	cosBAH_fade <- evaluate_fade_interval(
		fade_intervals, eval_type[1, "cosBAH_fade"], polynomial_range(polynomials, eval_type[1, "cosBAH"]))
	cosAHD_fade <- evaluate_fade_interval(
		fade_intervals, eval_type[1, "cosAHD_fade"], polynomial_range(polynomials, eval_type[1, "cosAHD"]))

	function(AHdist, cosBAH, cosAHD){
		(as.function(AHdist_poly)(      AHdist)) * fade_function(cosBAH_fade,     cosBAH) * fade_function(cosAHD_fade,     cosAHD) +
		fade_function(AHdist_short_fade, AHdist) * as.function(cosBAH_short_poly)(cosBAH) * fade_function(cosAHD_fade,     cosAHD) +
		fade_function(AHdist_long_fade,  AHdist) * as.function(cosBAH_long_poly)( cosBAH) * fade_function(cosAHD_fade,     cosAHD) +
		fade_function(AHdist_short_fade, AHdist) * fade_function(cosBAH_fade,     cosBAH) * as.function(cosAHD_short_poly)(cosAHD) +
		fade_function(AHdist_long_fade,  AHdist) * fade_function(cosBAH_fade,     cosBAH) * as.function(cosAHD_long_poly)( cosAHD)
	}
}


evaluate_function <- function(evaluation_function, AHdist, cosBAH, cosAHD){
	ddply(expand.grid(list(AHdist=AHdist, cosBAH=cosBAH, cosAHD=cosAHD)),
		.variables=c("AHdist", "cosBAH", "cosAHD"), function(df){
		data.frame(energy=evaluation_function(
			df[1,"AHdist"], df[1,"cosBAH"], df[1, "cosAHD"]))
	})
}
