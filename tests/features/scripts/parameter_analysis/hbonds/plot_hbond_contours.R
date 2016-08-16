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
#plot the predicted probability for each geometric dimension for given
#class of hydrogen bonds
#
#  cd minirosetta_database/scoring/score_functions/hbonds
#  ./scripts/plot_hbond_projection.R --help
#
#  see core::scoring::hbonds::hbond_compute_energy(...)
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
	"polynom",
	"rbenchmark")
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
	make_option(c("--bins"), action="store", default=150,
		help="Number of subdivision for plot[default: %default]"),
	make_option(c("--output_prefix"), default="polynomial_potential_",
		help="Prefix of filename where the plot is saved [default: %default]"),
	make_option(c("--output_suffix"), default=".png",
		help="Suffix of filename where the plot is saved [default: %default]"),
	make_option(c("--debug"), action="store_true", default=TRUE,
		help="Print debug information during script execution."))
opt <- parse_args(OptionParser(option_list=option_list))

opt <- parse_rosetta_database_path(opt)
parameter_set_path <- convert_parameter_set_into_database(
	opt$database, opt$parameter_set, opt$debug)


########### Extract evaluation functions ###############
if(opt$debug){
	cat("Extract evaluation functions\n")
}

poly_params <- get_polynomial_parameters(
  parameter_set_path, opt$don_chem_type, opt$acc_chem_type, opt$separation)

fade_params <- get_fade_parameters(
  parameter_set_path, opt$don_chem_type, opt$acc_chem_type, opt$separation)

############# Project onto each dimension #############
if(opt$debug) cat("Begin plotting contours\n")


#### AHdist vs cosBAH ############
d <- adply(
  expand.grid(
    AHdist=seq(poly_params[1,2], poly_params[1,3], length.out=opt$bins),
    cosBAH=seq(poly_params[2,2], poly_params[2,3], length.out=opt$bins),
    cosAHD=c(1)), 1, function(df){
      df$density =
        eval_hbond_energy(c(df$AHdist[1], df$cosBAH[1], df$cosAHD[1]), poly_params, fade_params)
      df
    })
if(opt$debug){
  print(summary(d))
}

output_fname <- paste(
  opt$output_prefix, "contour_AHdist_cosBAH_", opt$parameter_set, "_",
	opt$don_chem_type, "_", opt$acc_chem_type, "_", opt$separation,
	opt$output_suffix, sep="")

cat("Making plot of contour of probability density along the AHdist and cosBAH dimensions:\n")
cat("	parameter_set: ", opt$parameter_set, "\n")
cat("	don_chem_type: ", opt$don_chem_type, "\n")
cat("	acc_chem_type: ", opt$acc_chem_type, "\n")
cat("	separation:    ", opt$separation, "\n")
cat("	output plot:   ", output_fname, "\n")

p <- ggplot(d) + theme_bw() +
  geom_rect(aes(xmin=1.9, xmax=2.3, ymin=-Inf, ymax=Inf), alpha=.05) +
	geom_indicator(indicator="cosAHD=1") +
  geom_contour(aes(x=AHdist, y=cosBAH, z=density), breaks=seq(1,max(d$density), length.out=20)) +
  ggtitle(paste("Project HBond energy to the AHdist and cosBAH Dimensions\n",
    opt$parameter_set, ",",
    opt$don_chem_type, ", ",
    opt$acc_chem_type, ", ",
    opt$separation, sep=""))
ggsave(output_fname, height=8.3, width=10.8)


########### AHdist vs cosAHD #############
d <- adply(
  expand.grid(
    AHdist=seq(poly_params[1,2], poly_params[1,3], length.out=opt$bins),
    cosBAH=c(.5),
    cosAHD=seq(poly_params[2,2], poly_params[2,3], length.out=opt$bins)),
    1, function(df){
      df$density =
        eval_hbond_energy(c(df$AHdist[1], df$cosBAH[1], df$cosAHD[1]), poly_params, fade_params)
      df
    })
if(opt$debug){
  print(summary(d))
}

output_fname <- paste(
  opt$output_prefix, "contour_AHdist_cosAHD_", opt$parameter_set, "_",
	opt$don_chem_type, "_", opt$acc_chem_type, "_", opt$separation,
	opt$output_suffix, sep="")

cat("Making plot of contour of probability density along the AHdist and cosAHD dimensions:\n")
cat("	parameter_set: ", opt$parameter_set, "\n")
cat("	don_chem_type: ", opt$don_chem_type, "\n")
cat("	acc_chem_type: ", opt$acc_chem_type, "\n")
cat("	separation:    ", opt$separation, "\n")
cat("	output plot:   ", output_fname, "\n")

p <- ggplot(d) + theme_bw() +
  geom_rect(aes(xmin=1.9, xmax=2.3, ymin=-Inf, ymax=Inf), alpha=.05) +
	geom_indicator(indicator="cosBAH=.5") +
  geom_contour(aes(x=AHdist, y=cosAHD, z=density), breaks=seq(1,max(d$density), length.out=20)) +
  ggtitle(paste("Project HBond energy to the AHdist and cosAHD Dimensions\n",
    opt$parameter_set, ",",
    opt$don_chem_type, ", ",
    opt$acc_chem_type, ", ",
    opt$separation, sep=""))
ggsave(output_fname, height=8.3, width=10.8)

