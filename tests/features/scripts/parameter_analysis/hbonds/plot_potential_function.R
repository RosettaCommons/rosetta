#! /usr/bin/Rscript --vanilla
# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


###################################################################
#
#plot a polynomial from a hydrogen bond parameter set
#
#  cd minirosetta_database/scoring/score_functions/hbonds
#  ./scripts/plot_potential_function.R --help
#
###################################################################


# install packages: with install.packages("<PACKAGE_NAME>")

# to convert parametrizations into database
library("RSQLite")

# represent polynomials
library("polynom")

# plot functionality
library("ggplot2")

# parse command line options
library("optparse")

source("scripts/parameter_analysis/hbonds/methods/methods.R")

##### Setup command line options #########
option_list <- list(
	make_option(c("--database"),
		help="Path to rosetta_database [default: $ROSETTA3_DB]"),
  make_option(c("--parameter_set"), default="standard_params",
		help="Parameter set where the polynomial is defined [default: %default]"),
	make_option(c("--polynomial"), default="poly_AHdist_3",
		help="Name of the polynomial to plot [default: %default]"),
  make_option(c("--output_prefix"), default="build/polynomial_potential_",
		help="Prefix of filename where the plot is saved [default: %default]"),
	make_option(c("--output_suffix"), default=".png",
		help="Suffix of filename where the plot is saved [default: %default]"),
	make_option(c("--list"), action="store_true", default=FALSE,
		help="List all polynomials for specified parameter set and exit"))
opt <- parse_args(OptionParser(option_list=option_list))

opt$database <- valid_rosetta_database_path(opt$database)

convert_parameter_set_into_database(opt$database, opt$parameter_set)
polynomials <- get_polynomials(opt$database, opt$parameter_set)

# If the user request the --list option, make the table and then exit.
if(opt$list){
	display_polynomials(opt$parameter_set, polynomials)
	quit()
}

#########  Describe the plot that is going to be made ############
output_fname <- paste(
  opt$output_prefix, opt$parameter_set, "_", opt$polynomial, opt$output_suffix,
	sep="")

cat("Making plot of polynomial parameter:\n")
cat("    parameter_set: ", opt$parameter_set, "\n")
cat("    polynomial:    ", opt$polynomial, "\n")
cat("    output plot:   ", output_fname, "\n")

evaluated_polynomial <- evaluate_polynomial(polynomials, opt$polynomial)
dimension <- polynomials[polynomials$name==opt$polynomial, "dimension"]

######### Make the plot ################
p <- ggplot(evaluated_polynomial) + theme_bw() +
  geom_hline(yintercept=0, size=.7) +
  geom_line(aes(x, y), size=1.5) +
  ggtitle(paste("Hydrogen Bond Polynomial Evaluation:", opt$polynomial)) +
  scale_x_continuous(dimension) +
  scale_y_continuous("Raw Energy")
ggsave(output_fname)
