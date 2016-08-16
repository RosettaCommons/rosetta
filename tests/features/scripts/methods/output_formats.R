# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

all_output_formats <- read.csv(
	file.path(base_dir, "scripts/methods/output_formats.csv"), header=T, sep="\t")

make_output_formats_options_list <- function(output_formats){
	alply(output_formats, 1, function(output_format){
		make_option(
			c(paste("--", output_format$id, sep="")),
			action="store",
			type="logical",
			default=output_format$use_by_default,
			dest=as.character(output_format$id),
			help=paste("Generate output plots using the", output_format$id, "format.  [Default \"%default\"]"))
	})
}

get_output_formats <- function(options, output_formats){
	adply(output_formats, 1, function(output_format){
	       	if(options[as.character(output_format$id)] == TRUE){
			return(output_format)
		} else {
			return(data.frame())
		}
	})
}

get_output_formats_from_comparison <- function(ss_cmp, output_formats){
	ldply(ss_cmp$output_formats, function(of){
		if(of %in% output_formats$id){
			return(output_formats[output_formats$id == of,])
		} else {
			stop(paste("ERROR: Output format '", of, "' specified in the configuration script is not recognized. Here are the available output formats: \n  ", paste(output_formats$id, collapse="\n  ", sep="")))
		}
	})
}
