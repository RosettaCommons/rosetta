# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

ref_count <- function(a, b){
  data.frame(
    statistic_name = "ref_count",
    statistic = length(a),
    p.value = NA)
}

new_count <- function(a, b){
  data.frame(
    statistic_name = "new_count",
    statistic = length(b),
    p.value = NA)
}


two_sided_ttest <- function(a, b){
  tryCatch({
    z <- t.test(a, b)
  }, error=function(e) return(
    data.frame(
      statistic_name = "two_sided_t",
      statistic = NA,
      p.value = NA)))
  data.frame(
    statistic_name = "two_sided_t",
    statistic = z$statistic,
    p.value = z$p.value)
}

kolmogorov_smirnov_test_boot <- function(a, b){
  tryCatch({
    z <- ks.boot(a, b)
  }, error=function(e) {
		print(paste("Error: ks.boot failed with the followinng error", e, sep="", collapse=""))
		return(
			data.frame(
				statistic_name = "kolmogorov_smirnov_D",
				statistic = NA,
				p.value = NA))
	})
  data.frame(
    statistic_name = "kolmogorov_smirnov_D",
    statistic = z$ks$statistic,
    p.value = z$ks$p.value)
}


# Prefer to use the anderson_darling_2_sample comparison
# cf https://asaip.psu.edu/Articles/beware-the-kolmogorov-smirnov-test
kolmogorov_smirnov_test <- function(a, b){
	z <- NULL
  tryCatch({
    z <- ks.test(a, b, exact=F)
  }, error=function(e) {
		cat(paste("Error: ks.test failed with the followinng error:\n", e, sep="", collapse=""))
	})

	if(!is.null(z)){
	  return(data.frame(
	    statistic_name =factor("kolmogorov_smirnov_D"),
	    statistic = z$statistic,
	    p.value = NA))
	} else {
	  return(data.frame(
	    statistic_name = factor("kolmogorov_smirnov_D"),
	    statistic = NA,
	    p.value = NA))
	}
}

histogram_kl_divergence <- function(a, b, nbins=50){
	z <- NULL
	tryCatch({
		breaks <- seq(min(a, b), max(a, b), length.out=nbins)
		ad <- hist(a, breaks, plot=F)$density
		bd <- hist(b, breaks, plot=F)$density
		z <- data.frame(
			statistic_name=factor("KL Divergence"),
			statistic=sum(ad * log(ad / (bd+.0001)), na.rm = T),
			p.value = NA)
	}, error=function(e){
		cat(paste("Error: ks.test failed with the followinng error:\n", e, sep="", collapse=""))
	})

	if(!is.null(z)){
		return(z)
	} else {
		return(data.frame(
			stastistic_name=factor("KL Divergence"),
			statistic=NA,
			p.value=NA))
	}
}

# Here the inputs are probabilities over the sample space
#assume evenly spaced partition of points
smooth_kl_divergence <- function(x, ad, bd){
	z <- sum(ad*log((ad+.0001)/(bd+.0001)))
	data.frame(
		statistic_name=factor("KL Divergence"),
		statistic=z, na.rm=T,
		p.value=NA)
}


anderson_darling_2_sample <- function(a, b, nsamples=1000){
	require(adk) # this requires the adk package
	if(length(a) >= nsamples){
		print(length(a))
		print(nsamples)
		a_sub <- sample(a, nsamples)
	} else {
		a_sub <- a
	}

	if(length(b) >= nsamples){
		b_sub <- sample(b, nsamples)
	} else {
		b_sub <- b
	}

	z <- adk.test(a_sub,b_sub)
	data.frame(
		statistic_name=factor("Anderson Darling"),
		statistic=NA,
		p.value=z$adk[2,2]) # P-value, adjust for ties
}

#EXPERIMENTAL:
# require(earthmovdist)
earth_mover_distance_L1 <- function(a, b){


	a <- sample(a, 5000, replace=T)
	b <- sample(b, 5000, replace=T)

#
#	#The earthmovdist implementation of the earth mover distance
#	#requires the samples to be of the same length.  Is it better to
#	#subsample the longer one or resample the shorter one?
#	na <- length(a)
#	nb <- length(b)
#	if(na==nb){
#		# ok
#	} else if(na > nb){
#		b <- sample(b, na, replace=T)
#		n <- na
#	} else {
#		a <- sample(a, nb, replace=T)
#		n <- nb
#	}
	z <- emdL1(a, b)

	data.frame(
		statistic_name=factor("Earth Mover Distance L1"),
		statistic=exp(log(z) - 2434*log(length(a))),
		sample.size=length(a),
		p.value=NA)
}



comparison_statistics <- function(
	sample_sources,
	f,
	id.vars,
	measure.vars,
	comp_funs,
	verbose=FALSE
	){
	if(verbose){
		cat("The input data has the following columns:", paste(names(f), collapse=", "), "\n")
	}

	for(var in id.vars){
		if( !(var %in% names(f))){
			stop(paste("id.vars variable '", var, "', must be a column of f.\n\tnames(f) = c('", paste(names(f), collapse="', '"), "')", sep=""))
		}
	}

	measure.vars <- sapply(measure.vars, as.character)
	for(var in measure.vars){
		if(verbose){
			print(paste("checking if the measure variable '", var, "' is a column of f ...", sep=""))
		}
		if( !(var %in% names(f))){
			stop(paste("measure.vars variable '", var, "', must be a column of f.\n\tnames(f) = c('", paste(names(f), collapse="', '"), "')", sep=""))
		}
	}


	if(!is.character(comp_funs)){
		stop("The parameter comp_funs should be a vector of the names of comparison functions.")
	}

	# Compute a data.frame 'stat' with the following columns
	#   <id.vars>
	#   ref_sample_source
	#   new_sample_source
	#   statistic_name
	#   statistic
  #   p.value

	if(!("reference" %in% names(sample_sources))){
		stop("The provided sample_sources data.frame must have a boolean column 'reference' indicating which sample source comparisons should be made.

This is usually done by adding this information to the analysis_configuration. For example:
{
    \"sample_source_comparisons\" : [
        {
            \"sample_sources\" : [
                {
                    \"database_path\" : \"native_features.db3\",
                    \"id\" : \"Native\",
                    \"reference\" : true
		            },
                {
                    \"database_path\" : \"decoy_features.db3\",
                    \"id\" : \"Decoy\",
                    \"reference\" : false
		            }
            ],
            \"analysis_scripts\" : [ ... ],
            \"output_formats\" : [ ... ]
	      }
}\n")
	}

	ref_sample_sources <- as.character(
		sample_sources[sample_sources$reference,"sample_source"])
	new_sample_sources <- as.character(
		sample_sources[!sample_sources$reference,"sample_source"])

	if(length(ref_sample_sources) == 0){
		cat("WARNING: No 'reference sample sample sources were specified.\n")
  }

	if(length(new_sample_sources) == 0){
		cat("WARNING: No non-reference sample sample sources were specified.\n")
  }

  cat("Comparing dimension(s): ", paste(measure.vars, collapse=" "), "
Grouping by: ", paste(id.vars, collapse=" "), "
	ref: ", paste(ref_sample_sources, collapse=", "), "
	new: ", paste(new_sample_sources, collapse=", "), "\n", sep="")

	# since there are usually more groups of id.vars than sample
	# sources, make id.vars the outer loop
  stats <- ddply(f, .variables = id.vars, function(sub_f){
		if(verbose){
			cat("Doing comparison for group: '", paste(lapply(sub_f[1,id.vars], as.character), collapse="', '"), "'\n", sep="")
			print(summary(sub_f))
		}

		timing <- system.time({
    	z <- ddply(
				sub_f[sub_f$sample_source %in% ref_sample_sources,],
				c("ref_sample_source" = "sample_source"), function(ref_f){

				ddply(
					sub_f[sub_f$sample_source %in% new_sample_sources,],
					c("new_sample_source" = "sample_source"), function(new_f){

					ldply(comp_funs, function(comp_fun){
						get(comp_fun)(ref_f[,measure.vars], new_f[,measure.vars])
					}) # comp_fun
    	  }) # ref_f
    	}) # new_f
		}) # timing
		cat(as.character(round(timing[3],2)), "s  ", sep="")
		z
  }) # sub_f
	cat("\n")

	if(length(id.vars) > 0){
		cast_formula_str <-
			paste(
				paste(as.character(id.vars), collapse=" + "),
				" + ref_sample_source + new_sample_source ~ statistic_name", sep="")
	} else {
		cast_formula_str <- "ref_sample_source + new_sample_source ~ statistic_name"
	}
	cast_formula <- as.formula(cast_formula_str)
	cat("Casting result as: ", cast_formula_str, "\n", sep="")
	cast(stats, cast_formula, value="statistic")
}


# Evaluate a two sample tests between different classes of samples
# conditional on distinct groups of identifying variables.
#   The class.vars is used to identify the classes of the samples
#   The id.vars is used to split the classes in to comparison groups
smooth_comparison_statistics <- function(
	dens, id.vars, comp_fun){

	for(var in id.vars){
		if( !(var %in% names(f))){
			stop(paste("id.vars variable '", var, "', must be a column of f.\n\tnames(f) = c('", paste(names(f), collapse="', '"), "')", sep=""))
		}
	}

  ddply(dens, .variables = id.vars, function(sub_dens){
    ddply(sub_dens, c("ref_sample_source" = "sample_source"), function(ref_dens){
      ddply(sub_dens, c("new_sample_source" = "sample_source"), function(new_dens){
				if(as.numeric(new_dens$sample_source[1]) <= as.numeric(ref_dens$sample_source[1])) {
					return(data.frame())
				}
        comp_fun(ref_dens$x, ref_dens$y, new_dens$y)
      })
    })
  })
}



# Evaluate a two sample tests between different classes of samples
# conditional on distinct groups of identifying variables.
#   The class.vars is used to identify the classes of the samples
#   The id.vars is used to split the classes in to comparison groups
cross_validate_statistics <- function(
	f, cv.var, id.vars, measure.vars, comp_fun){
	if( !(cv.var %in% names(f))){
		stop(paste("cv.vars variable '", var, "', must be a column of f.\n\tnames(f) = c('", paste(names(f), collapse="', '"), "')", sep=""))
	}

	for(var in id.vars){
		if( !(var %in% names(f))){
			stop(paste("id.vars variable '", var, "', must be a column of f.\n\tnames(f) = c('", paste(names(f), collapse="', '"), "')", sep=""))
		}
	}

	for(var in measure.vars){
		if( !(var %in% names(f))){
			stop(paste("measure.vars variable '", class.var, "', must be a column of f.\n\tnames(f) = c('", paste(names(f), collapse="', '"), "')", sep=""))
		}
	}

  ddply(f, .variables = id.vars, function(sub_f){
    ddply(sub_f, .variables = cv.var, function(test_f){
			cv.group <- test_f[1,cv.var]
			rest_f <- sub_f[sub_f[,cv.var] != cv.group,]
			comp_fun(rest_f[,measure.vars], test_f[,measure.vars])
    })
  })
}
