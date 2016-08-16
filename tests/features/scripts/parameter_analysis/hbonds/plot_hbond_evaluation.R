#! /usr/bin/Rscript --vanilla
# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


####################################################################
#
#plot the hbond evaluation function for given class of hydrogen bonds
#
####################################################################

# Get path to 'rosetta_tests/features/scripts/parameter_analysis', where this script lives
command_args <- commandArgs(trailingOnly = FALSE)
base_dir <- dirname(substring(command_args[grep("--file=", command_args)], 8))

# load_packages() will kindly ask to install the packages if they are missing
source(paste(base_dir, "../../methods/load_packages.R", sep="/"))
libraries <- c(
	"optparse",
	"ggplot2",
	"RSQLite",
	"inline",
	"plyr",
	"Rcpp",
	"polynom")
#	"rbenchmark")
load_packages(libraries)

includes <- c(
	"methods/methods.R",
	"methods/rcpp_hbond_evaluation.R",
	"../methods/methods.R",
	"../../methods/ggplot2_geom_indicator.R")
for(inc in includes) source(paste(base_dir, inc, sep="/"))

##### Setup command line options #########
option_list <- list(
	make_option(c("--database"),
		help="Path to rosetta_database [default: $ROSETTA3_DB]"),
	make_option(c("--parameter_set"), default="standard_params",
		help="Parameter set [default: %default]"),
	make_option(c("--don_chem_type"), default="hbdon_PBA",
		help="Donor chemical type [default: %default]"),
	make_option(c("--acc_chem_type"), default="hbacc_PBA",
		help="Acceptor chemical type [default: %default]"),
	make_option(c("--separation"), default="seq_sep_other",
		help="Sequence separation, used to determine secondary structure [default: %default]"),
	make_option(c("--output_prefix"), default="polynomial_potential_",
		help="Prefix of filename where the plot is saved [default: %default]"),
	make_option(c("--output_suffix"), default=".png",
		help="Suffix of filename where the plot is saved [default: %default]"),
	make_option(c("--list"), action="store_true", default=FALSE,
		help="List all polynomials for specified parameter set and exit"),
	make_option(c("--debug"), action="store_true", default=TRUE,
		help="Print debug information during script execution."))
opt <- parse_args(OptionParser(option_list=option_list))

opt <- parse_rosetta_database_path(opt)
parameter_set_path <- convert_parameter_set_into_database(
	opt$database, opt$parameter_set, opt$debug)

########### Extract evaluation functions ###############
if(opt$debug) cat("Extract evaluation functions\n")

poly_params <- get_polynomial_parameters(
  parameter_set_path,
	opt$don_chem_type, opt$acc_chem_type, opt$separation, opt$debug)

fade_params <- get_fade_parameters(
  parameter_set_path, opt$don_chem_type, opt$acc_chem_type, opt$separation, opt$debug)

if(opt$debug)cat("evaluating polynomials...\n")
AHdist_poly <- evaluate_polynomial(poly_params[1,], "AHdist")
cosBAH_short_poly <- evaluate_polynomial(poly_params[2,],"cosBAH_short")
cosBAH_long_poly <- evaluate_polynomial(poly_params[3,], "cosBAH_long")
cosAHD_short_poly <- evaluate_polynomial(poly_params[4,], "cosAHD_short")
cosAHD_long_poly <- evaluate_polynomial(poly_params[5,], "cosAHD_long")

if(opt$debug) cat("evaluating fade intervals...\n")
AHdist_short_fade <- evaluate_fade_interval(
  fade_params[1,], poly_params[1,2], poly_params[1,3], "AHdist_short_fade")
AHdist_long_fade <- evaluate_fade_interval(
  fade_params[2,], poly_params[1,2], poly_params[1,3], "AHdist_long_fade")
cosBAH_fade <- evaluate_fade_interval(
  fade_params[3,], poly_params[2,2], poly_params[2,3], "cosBAH_fade")
cosAHD_fade <- evaluate_fade_interval(
  fade_params[4,], poly_params[4,2], poly_params[2,3], "cosAHD_fade")

if(opt$debug) cat("Begin generating plots ... \n")
#########  Describe the plot that is going to be made ############
output_fname <- paste(
	opt$output_prefix, opt$parameter_set, "_",
	opt$don_chem_type, "_", opt$acc_chem_type, "_", opt$separation,
	opt$output_suffix, sep="")

cat("Making plot of polynomial parameter:\n")
cat("	parameter_set: ", opt$parameter_set, "\n")
cat("	don_chem_type: ", opt$don_chem_type, "\n")
cat("	acc_chem_type: ", opt$acc_chem_type, "\n")
cat("	separation:    ", opt$separation, "\n")
cat("	output plot:   ", output_fname, "\n")

########### Assemble plot data ####################
if(opt$debug){
	cat("Assemble plot data\n")
	cat("	it has the following columns plus 'term', 'geo_dim', and 'type'\n")
	print(summary(AHdist_poly))
	print(summary(AHdist_short_fade))
}



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
if(opt$debug){
	cat("actually make the plot of the polynomial parameters.\n")
}

p <- ggplot(evaluation_function) + theme_bw() +
	geom_hline(yintercept=0, size=.5) +
	geom_line(aes(x, y, size=type)) +
	geom_indicator(aes(indicator=name), group=1) +
	facet_grid(term ~ geo_dim, scales="free_x") +
	ggtitle(paste("Hydrogen Bond Evaluation Function: ", opt$don_chem_type, ", ", opt$acc_chem_type, ", ", opt$separation, "
1) evaluate each cell  2) multiply across  3) add results", sep="")) +
	scale_x_continuous("Geometric Degree Of Freedom") +
	scale_y_continuous("Energy Term", limit=c(-.7, 1.2)) +
	scale_size_manual(values = c("poly"=1.5, "fade"=.7))
ggsave(output_fname, height=8.3, width=10.8)


################# make plots of energy along each dimension ################

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

dofs <- data.frame(
  AHdist=seq(poly_params[1,2], poly_params[1,3], length.out=500), cosBAH=c(.5), cosAHD=c(1))

d <- adply(dofs, 1, function(d){
  d$density = eval_hbond_energy(data.matrix(d), poly_params, fade_params)
  d
})

summary(d)


p <- ggplot(d) + theme_bw() +
	geom_line(aes(x=AHdist, y=density)) +
	geom_indicator(aes(indicator=interaction(paste("cosBAH=", cosBAH, sep=""), paste("cosAHD=", cosAHD, sep=""), sep=" "))) +
	ggtitle(paste("Evaluate HBond Energy Along the AHdist Dimension\n", opt$parameter_set, ",", opt$don_chem_type, ", ", opt$acc_chem_type, ", ", opt$separation, sep=""))
#        scale_x_continuous(limits=c(poly_params[1,2], poly_params[1,3]))
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

dofs <- expand.grid(AHdist=c(1.8), cosBAH=seq(0, 1, length.out=500), cosAHD=c(1))
d <- adply(dofs, 1, function(d){
  d$density = eval_hbond_energy(data.matrix(d), poly_params, fade_params)
  d
})

p <- ggplot(d) + theme_bw() +
	geom_line(aes(x=cosBAH, y=density)) +
  geom_indicator(aes(indicator=interaction(paste("AHdist=", AHdist, sep=""), paste("cosAHD=", cosAHD, sep=""), sep=" "))) +
	ggtitle(paste("Evaluate HBond Energy Along the cosBAH Dimension\n", opt$parameter_set, ",", opt$don_chem_type, ", ", opt$acc_chem_type, ", ", opt$separation, sep=""))
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

dofs <- expand.grid(
  AHdist=c(1.9),
  cosBAH=c(.5),
  cosAHD=seq(0, 1, length.out=500))
d <- adply(dofs, 1, function(d){
  d$density = eval_hbond_energy(data.matrix(d), poly_params, fade_params)
  d
})

d$AHdist <- factor(d$AHdist)
summary(d)

p <- ggplot(d) + theme_bw() +
  geom_line(aes(x=cosAHD, y=density, color=AHdist)) +
  geom_indicator(aes(indicator=interaction(paste("AHdist=", AHdist, sep=""), paste("cosBAH=", cosBAH, sep=""), sep=" "))) +
	ggtitle(paste("Evaluate HBond Energy Along the cosAHD Dimension\n", opt$parameter_set, ",", opt$don_chem_type, ", ", opt$acc_chem_type, ", ", opt$separation, sep=""))
ggsave(output_fname, height=8.3, width=10.8)



output_fname <- paste(
	opt$output_prefix, "cosAHD_AHdist", opt$parameter_set, "_",
	opt$don_chem_type, "_", opt$acc_chem_type, "_", opt$separation,
	opt$output_suffix, sep="")

cat("Making plot of energy along cosAHD and AHdist dimension:\n")
cat("	parameter_set: ", opt$parameter_set, "\n")
cat("	don_chem_type: ", opt$don_chem_type, "\n")
cat("	acc_chem_type: ", opt$acc_chem_type, "\n")
cat("	separation:    ", opt$separation, "\n")
cat("	output plot:   ", output_fname, "\n")


dofs <- expand.grid(
  AHdist=seq(poly_params[1,2], poly_params[1,3], length.out=200),
  cosBAH=c(.5),
  cosAHD=seq(0, 1, length.out=200))
d <- adply(dofs, 1, function(d){
  d$density = eval_hbond_energy(data.matrix(d), poly_params, fade_params)
  d
})

summary(d)

p <- ggplot(d) + theme_bw() +
	geom_tile(aes(x=cosAHD, y=AHdist, fill=density)) +
  geom_indicator(aes(indicator=paste("cosBAH=", cosBAH, sep=""))) +
	ggtitle(paste("Evaluate HBond Energy Along the cosAHD and AHdist Dimensions\n", opt$parameter_set, ",", opt$don_chem_type, ", ", opt$acc_chem_type, ", ", opt$separation, sep=""))
ggsave(output_fname, height=8.3, width=10.8)
