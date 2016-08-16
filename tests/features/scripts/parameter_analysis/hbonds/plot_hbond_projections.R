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
	"cubature",
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
	make_option(c("--tolerance"), action="store", default=1e-3,
		help="The tolerance parameter for the numeric integration. [default: %default]"),
	make_option(c("--bins"), action="store", default=50,
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
if(opt$debug) cat("Extract evaluation functions\n")

poly_params <- get_polynomial_parameters(
  parameter_set_path, opt$don_chem_type, opt$acc_chem_type, opt$separation)

fade_params <- get_fade_parameters(
  parameter_set_path, opt$don_chem_type, opt$acc_chem_type, opt$separation)

############# Project onto each dimension #############
if(opt$debug) cat("Begin plotting projections\n")

d <- adply(data.frame(AHdist=seq(poly_params[1,2], poly_params[1,3], length.out=opt$bins)),
	1, function(AHdist){
  	if(opt$debug){
          cat("AHdist=", AHdist[[1]],"\n")
  	}
  	data.frame(
    	AHdist = AHdist,
    	density = adaptIntegrate(function(x){
          eval_hbond_energy(
            c(AHdist[[1]], x[1], x[2]), poly_params, fade_params)},
          c(0.0, 0.0), c(1.0, 1.0), tol=opt$tolerance)$integral)
})

if(opt$debug){
  print(summary(d))
}

output_fname <- paste(
  opt$output_prefix, "AHdist_projection_", opt$parameter_set, "_",
	opt$don_chem_type, "_", opt$acc_chem_type, "_", opt$separation,
	opt$output_suffix, sep="")

cat("Making plot of projection of probability density along the AHdist dimension:\n")
cat("	parameter_set: ", opt$parameter_set, "\n")
cat("	don_chem_type: ", opt$don_chem_type, "\n")
cat("	acc_chem_type: ", opt$acc_chem_type, "\n")
cat("	separation:    ", opt$separation, "\n")
cat("	output plot:   ", output_fname, "\n")

p <- ggplot(d) + theme_bw() +
  geom_line(aes(x=AHdist, y=density), size=1.2) +
	ggtitle(paste("Project HBond energy to the AHdist Dimension\n",
    opt$parameter_set, ",",
    opt$don_chem_type, ", ",
    opt$acc_chem_type, ", ",
    opt$separation, sep=""))
ggsave(output_fname, height=8.3, width=10.8)


d <-adply(data.frame(cosBAH=seq(0,1, length.out=opt$bins)),
  1, function(cosBAH){
    if(opt$debug){
      cat("cosBAH=", cosBAH[[1]],"\n")
    }
    data.frame(
      cosBAH = cosBAH,
      density = adaptIntegrate(
        function(x){
          eval_hbond_energy(
            c(x[1], cosBAH[[1]], x[2]), poly_params, fade_params)},
        c(poly_params[1,2], 0), c(poly_params[1,3],1), tol=opt$tolerance)$integral)
})

output_fname <- paste(
  opt$output_prefix, "cosBAH_projection_", opt$parameter_set, "_",
  opt$don_chem_type, "_", opt$acc_chem_type, "_", opt$separation,
	opt$output_suffix, sep="")

cat("Making plot of projection of probability density along the cosBAH dimension:\n")
cat("	parameter_set: ", opt$parameter_set, "\n")
cat("	don_chem_type: ", opt$don_chem_type, "\n")
cat("	acc_chem_type: ", opt$acc_chem_type, "\n")
cat("	separation:    ", opt$separation, "\n")
cat("	output plot:   ", output_fname, "\n")

p <- ggplot(d) + theme_bw() +
  geom_line(aes(x=cosBAH, y=density), size=1.2) +
	ggtitle(paste("Project HBond energy to the cosBAH Dimension\n",
    opt$parameter_set, ",",
    opt$don_chem_type, ", ",
    opt$acc_chem_type, ", ",
    opt$separation, sep=""))
ggsave(output_fname, height=8.3, width=10.8)




d <- adply(data.frame(cosAHD=seq(0,1, length.out=opt$bins)),
  1, function(cosAHD){
    if(opt$debug){
      cat("cosAHD=", cosAHD[[1]],"\n")
    }
    data.frame(
      cosAHD = cosAHD,
      density = adaptIntegrate(function(x){
        eval_hbond_energy(c(x[1], x[2], cosAHD[[1]]), poly_params, fade_params)},
        c(1.4, 0), c(3, 1), tol=opt$tolerance)$integral)
  })

output_fname <- paste(
  opt$output_prefix, "cosAHD_projection_", opt$parameter_set, "_",
  opt$don_chem_type, "_", opt$acc_chem_type, "_", opt$separation,
  opt$output_suffix, sep="")

cat("Making plot of projection of probability density along the cosAHD dimension:\n")
cat("	parameter_set: ", opt$parameter_set, "\n")
cat("	don_chem_type: ", opt$don_chem_type, "\n")
cat("	acc_chem_type: ", opt$acc_chem_type, "\n")
cat("	separation:    ", opt$separation, "\n")
cat("	output plot:   ", output_fname, "\n")

p <- ggplot(d) + theme_bw() +
  geom_line(aes(x=cosAHD, y=density), size=1.2) +
	ggtitle(paste("Project HBond energy to the cosAHD Dimension\n",
    opt$parameter_set, ",",
    opt$don_chem_type, ", ",
    opt$acc_chem_type, ", ",
    opt$separation, sep=""))
ggsave(output_fname, height=8.3, width=10.8)



