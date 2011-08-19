#! /usr/bin/Rscript --vanilla
# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


####################################################################
#
#plot the hbond evaluation function for given class of hydrogen bonds
#
#  cd minirosetta_database/scoring/score_functions/hbonds
#  ./scripts/plot_hbond_evaluation.R --help
#
#  see core::scoring::hbonds::hbond_compute_energy(...)
#
####################################################################


# install packages: with install.packages("<PACKAGE_NAME>")

# to convert parametrizations into database
library("RSQLite")

# represent polynomials
library("polynom")

# plot functionality
library("ggplot2")

# parse command line options
library("optparse")

# extract information from parameter database
source("scripts/methods.R")

# custom plot layer for printing indicator label in the plot cell
source("scripts/ggplot2_geom_indicator.R")

##### Setup command line options #########
option_list <- list(
	make_option(c("--parameter_set"), default="standard_params",
		help="Parameter set [default: %default]"),
	make_option(c("--don_chem_type"), default="hbdon_PBA",
		help="Donor chemical type [default: %default]"),
	make_option(c("--acc_chem_type"), default="hbacc_PBA",
		help="Acceptor chemical type [default: %default]"),
	make_option(c("--separation"), default="seq_sep_other",
		help="Acceptor chemical type [default: %default]"),
	make_option(c("--output_prefix"), default="polynomial_potential_",
		help="Prefix of filename where the plot is saved [default: %default]"),
	make_option(c("--output_suffix"), default=".png",
		help="Suffix of filename where the plot is saved [default: %default]"),
	make_option(c("--list"), action="store_true", default=FALSE,
		help="List all polynomials for specified parameter set and exit"))
opt <- parse_args(OptionParser(option_list=option_list))


##### Validate input ##########
cwd <- getwd();
if(substr(cwd, nchar(cwd)-29,nchar(cwd)) != "scoring/score_functions/hbonds"){
	stop("ERROR: Please run this script from $ROSETTA3_DB/scoring/score_functions/hbonds.")
}

if(!file.exists(opt$parameter_set)){
	stop(paste("Parameter set, ", opt$parameter_set, ", is not recognized.",sep=""))
}

########### Extract evaluation functions ###############
convert_parameter_set_into_database(opt$parameter_set)

polynomials <- get_polynomials(opt$parameter_set)
fade_intervals <- get_fade_intervals(opt$parameter_set)

eval_type <- get_evaluation_type(
	opt$parameter_set, opt$don_chem_type, opt$acc_chem_type, opt$separation)

AHdist_poly <- evaluate_polynomial(polynomials, eval_type[1,"AHdist"])
cosBAH_short_poly <- evaluate_polynomial(polynomials, eval_type[1,"cosBAH_short"])
cosBAH_long_poly <- evaluate_polynomial(polynomials, eval_type[1,"cosBAH_long"])
cosAHD_short_poly <- evaluate_polynomial(polynomials, eval_type[1,"cosAHD_short"])
cosAHD_long_poly <- evaluate_polynomial(polynomials, eval_type[1,"cosAHD_long"])

AHdist_short_fade <- evaluate_fade_interval(
	fade_intervals, eval_type[1,"AHdist_short_fade"], polynomial_range(polynomials, eval_type[1, "AHdist"]))
AHdist_long_fade <- evaluate_fade_interval(
	fade_intervals, eval_type[1,"AHdist_long_fade"], polynomial_range(polynomials, eval_type[1, "AHdist"]))
cosBAH_fade <- evaluate_fade_interval(
	fade_intervals, eval_type[1, "cosBAH_fade"], polynomial_range(polynomials, eval_type[1, "cosBAH"]))
cosAHD_fade <- evaluate_fade_interval(
	fade_intervals, eval_type[1, "cosAHD_fade"], polynomial_range(polynomials, eval_type[1, "cosAHD"]))

#########  Describe the plot that is going to be made ############
output_fname <- paste(
	opt$output_prefix, opt$parameter_set, "_",
	opt$don_chem_type, "_", opt$acc_chem_type, "_", opt$separation,
	opt$output_suffix, sep="")

cat("Making plot of polynomial parameter:\n")
cat("		parameter_set: ", opt$parameter_set, "\n")
cat("		don_chem_type: ", opt$don_chem_type, "\n")
cat("		acc_chem_type: ", opt$acc_chem_type, "\n")
cat("		separation:    ", opt$separation, "\n")
cat("  	output plot:   ", output_fname, "\n")

########### Assemble plot data ####################

evaluation_function <-rbind(
	data.frame(term=factor("AHdist"),       geo_dim=factor("AHdist"), type=factor("poly"), AHdist_poly),
	data.frame(term=factor("AHdist"),       geo_dim=factor("cosBAH"), type=factor("fade"), cosBAH_fade),
	data.frame(term=factor("AHdist"),       geo_dim=factor("cosAHD"), type=factor("fade"), cosAHD_fade),
	data.frame(term=factor("cosBAH_short"), geo_dim=factor("AHdist"), type=factor("fade"), AHdist_short_fade),
	data.frame(term=factor("cosBAH_short"), geo_dim=factor("cosBAH"), type=factor("poly"), cosBAH_short_poly),
	data.frame(term=factor("cosBAH_short"), geo_dim=factor("cosAHD"), type=factor("fade"), cosAHD_fade),
	data.frame(term=factor("cosBAH_long"),  geo_dim=factor("AHdist"), type=factor("fade"), AHdist_long_fade),
	data.frame(term=factor("cosBAH_long"),  geo_dim=factor("cosBAH"), type=factor("poly"), cosBAH_long_poly),
	data.frame(term=factor("cosBAH_long"),  geo_dim=factor("cosAHD"), type=factor("fade"), cosAHD_fade),
	data.frame(term=factor("cosAHD_short"), geo_dim=factor("AHdist"), type=factor("fade"), AHdist_short_fade),
	data.frame(term=factor("cosAHD_short"), geo_dim=factor("cosBAH"), type=factor("fade"), cosBAH_fade),
	data.frame(term=factor("cosAHD_short"), geo_dim=factor("cosAHD"), type=factor("poly"), cosAHD_short_poly),
	data.frame(term=factor("cosAHD_long"),  geo_dim=factor("AHdist"), type=factor("fade"), AHdist_long_fade),
	data.frame(term=factor("cosAHD_long"),  geo_dim=factor("cosBAH"), type=factor("fade"), cosBAH_fade),
	data.frame(term=factor("cosAHD_long"),  geo_dim=factor("cosAHD"), type=factor("poly"), cosAHD_long_poly))

######### Make the plot ################

p <- ggplot(evaluation_function) + theme_bw() +
	geom_hline(yintercept=0, size=.5) +
	geom_line(aes(x, y, size=type)) +
	geom_indicator(aes(indicator=name)) +
	facet_grid(term ~ geo_dim, scales="free_x") +
	opts(title=paste("Hydrogen Bond Evaluation Function: ", opt$don_chem_type, ", ", opt$acc_chem_type, ", ", opt$separation, "
1) evaluate each cell  2) multiply across  3) add results", sep="")) +
	scale_x_continuous("Geometric Degree Of Freedom") +
	scale_y_continuous("Energy Term", limit=c(-.7, 1.2)) +
	scale_size_manual(values = c("poly"=1.5, "fade"=.7))
ggsave(output_fname, height=8.3, width=10.8)


################# make plots of energy along each dimension ################3

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

evaluation_function <-function(AHdist, cosBAH, cosAHD){
	(as.function(AHdist_poly)(      AHdist)) * fade_function(cosBAH_fade,     cosBAH) * fade_function(cosAHD_fade,     cosAHD) +
	fade_function(AHdist_short_fade, AHdist) * as.function(cosBAH_short_poly)(cosBAH) * fade_function(cosAHD_fade,     cosAHD) +
	fade_function(AHdist_long_fade,  AHdist) * as.function(cosBAH_long_poly)( cosBAH) * fade_function(cosAHD_fade,     cosAHD) +
	fade_function(AHdist_short_fade, AHdist) * fade_function(cosBAH_fade,     cosBAH) * as.function(cosAHD_short_poly)(cosAHD) +
	fade_function(AHdist_long_fade,  AHdist) * fade_function(cosBAH_fade,     cosBAH) * as.function(cosAHD_long_poly)( cosAHD)
}


output_fname <- paste(
	opt$output_prefix, "AHdist_", opt$parameter_set, "_",
	opt$don_chem_type, "_", opt$acc_chem_type, "_", opt$separation,
	opt$output_suffix, sep="")

cat("Making plot of energy along the AHdist dimension:\n")
cat("	parameter_set: ", opt$parameter_set, "\n")
cat("	don_chem_type: ", opt$don_chem_type, "\n")
cat("	acc_chem_type: ", opt$acc_chem_type, "\n")
cat("	separation:    ", opt$separation, "\n")
cat("	output plot:   ", output_fname, "\n")

d <- evaluate_function(
	evaluation_function,
	seq(polynomials[polynomials$name == eval_type[1,"AHdist"], "xmin"],	polynomials[polynomials$name == eval_type[1,"AHdist"], "xmax"], length.out=500),
	c(.5),c(1))

p <- ggplot(d) + theme_bw() +
	geom_line(aes(x=AHdist, y=energy)) +
	geom_indicator(aes(indicator=interaction(paste("cosBAH=", cosBAH, sep=""), paste("cosAHD=", cosAHD, sep=""), sep=" "))) +
	opts(title=paste("Evaluate HBond Energy Along the AHdist Dimension\n", opt$parameter_set, ",", opt$don_chem_type, ", ", opt$acc_chem_type, ", ", opt$separation, sep=""))
ggsave(output_fname, height=8.3, width=10.8)


output_fname <- paste(
	opt$output_prefix, "cosBAH_", opt$parameter_set, "_",
	opt$don_chem_type, "_", opt$acc_chem_type, "_", opt$separation,
	opt$output_suffix, sep="")

cat("Making plot of energy along the cosBAH dimension:\n")
cat("	parameter_set: ", opt$parameter_set, "\n")
cat("	don_chem_type: ", opt$don_chem_type, "\n")
cat("	acc_chem_type: ", opt$acc_chem_type, "\n")
cat("	separation:    ", opt$separation, "\n")
cat("	output plot:   ", output_fname, "\n")

d <- evaluate_function(
	evaluation_function,
	c(1.8),
	seq(0,1, length.out=500),
	c(1))

p <- ggplot(d) + theme_bw() +
	geom_line(aes(x=cosBAH, y=energy)) + geom_indicator(aes(indicator=interaction(paste("AHdist=", AHdist, sep=""), paste("cosAHD=", cosAHD, sep=""), sep=" "))) +
	opts(title=paste("Evaluate HBond Energy Along the cosBAH Dimension\n", opt$parameter_set, ",", opt$don_chem_type, ", ", opt$acc_chem_type, ", ", opt$separation, sep=""))
ggsave(output_fname, height=8.3, width=10.8)


output_fname <- paste(
	opt$output_prefix, "cosAHD_", opt$parameter_set, "_",
	opt$don_chem_type, "_", opt$acc_chem_type, "_", opt$separation,
	opt$output_suffix, sep="")

cat("Making plot of energy along cosAHD dimension:\n")
cat("	parameter_set: ", opt$parameter_set, "\n")
cat("	don_chem_type: ", opt$don_chem_type, "\n")
cat("	acc_chem_type: ", opt$acc_chem_type, "\n")
cat("	separation:    ", opt$separation, "\n")
cat("	output plot:   ", output_fname, "\n")

d <- evaluate_function(
	evaluation_function,
	c(1.8),
	c(.5),
	seq(0,1, length.out=500))

p <- ggplot(d) + theme_bw() +
	geom_line(aes(x=cosAHD, y=energy)) + geom_indicator(aes(indicator=interaction(paste("AHdist=", AHdist, sep=""), paste("cosBAH=", cosBAH, sep=""), sep=" "))) +
	opts(title=paste("Evaluate HBond Energy Along the cosAHD Dimension\n", opt$parameter_set, ",", opt$don_chem_type, ", ", opt$acc_chem_type, ", ", opt$separation, sep=""))
ggsave(output_fname, height=8.3, width=10.8)
