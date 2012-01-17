# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

kolmogorov_smirnov_test <- function(a, b){
  tryCatch({
    z <- ks.boot(a, b)
  }, error=function(e) return(
    data.frame(
      statistic_name = "kolmogorov_smirnov_D",
      statistic = NA,
      p.value = NA)))
  data.frame(
    statistic_name = "kolmogorov_smirnov_D",
    statistic = z$ks$statistic,
    p.value = z$ks$p.value)
}

histogram_kl_divergence <- function(a, b, nbins=50){
	breaks <- seq(min(a, b), max(a, b), length.out=nbins)
	ad <- hist(a, breaks, plot=F)$density
	bd <- hist(b, breaks, plot=F)$density
	data.frame(
		statistic_name=factor("KL Divergence"),
		statistic=sum(ad * log(ad / (bd+.0001)), na.rm = T),
		p.value = NA)
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



# Evaluate a two sample tests between different classes of samples
# conditional on distinct groups of identifying variables.
#   The class.vars is used to identify the classes of the samples
#   The id.vars is used to split the classes in to comparison groups
comparison_statistics <- function(
	f, id.vars, measure.vars, comp_fun){

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
    ddply(sub_f, c("ref_sample_source" = "sample_source"), function(ref_f){
      ddply(sub_f, c("new_sample_source" = "sample_source"), function(new_f){
				if(as.numeric(new_f$sample_source[1]) <= as.numeric(ref_f$sample_source[1])) {
					return(data.frame())
				}
        comp_fun(ref_f[,measure.vars], new_f[,measure.vars])
      })
    })
  })
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


# EXAMPLE USAGE  (work in progress)
#
#small_f <- sample_rows(f, 1000)
#comps <- rbind(
#  stats(f, "sample_source", c("don_chem_type", "acc_chem_type"), "ADdist", two_sided_ttest)
#  stats(f, "sample_source", c("don_chem_type", "acc_chem_type"), "ADdist", kolmogorov_smirnov_test)
#)
#
#comps[comps$ref_sample_source != comps$new_sample_source,]
#tt <- comps[comps$ref_sample_source != comps$new_sample_source,]
#
#system.time(z <- stats(f, "sample_source", c("don_chem_type", "acc_chem_type"), "ADdist", kolmogorov_smirnov_test))
#ks <- z[comps$ref_sample_source != comps$new_sample_source,]
#
#tt_vs_tt_sc12_sc12c <- data.frame(
#  don_chem_type    = tt[tt$ref_sample_source == "top4400_1108_apltest4c" & tt$new_sample_source =="t4400_1108_sc12", "don_chem_type"],
#  acc_chem_type    = tt[tt$ref_sample_source == "top4400_1108_apltest4c" & tt$new_sample_source =="t4400_1108_sc12", "acc_chem_type"],
#  tt.p.value.sc12  = tt[tt$ref_sample_source == "top4400_1108_apltest4c" & tt$new_sample_source =="t4400_1108_sc12", "p.value"],
#  tt.p.value.sc12c = tt[tt$ref_sample_source == "top4400_1108_apltest4c" & tt$new_sample_source == "t4400_sc12corr", "p.value"])
#qplot(data=tt_vs_tt_sc12_sc12c, x = interaction(don_chem_type, acc_chem_type), y=tt.p.value.sc12c/tt.p.value.sc12, geom="bar", log="y") + coord_flip()
#ggsave("/tmp/tt_vs_tt_sc12_sc12c.png", height=15, width=7)
#
#df <- rbind(tt, ks)
#
#plot_id <- "ADdist_compare_tt_ks"
#p <- ggplot(df[df$ref_sample_source == "top4400_1108_apltest4c" & df$new_sample_source != "top4400_1108_apltest4c",]) + theme_bw() +
#  geom_bar(aes(x=interaction(don_chem_type, acc_chem_type), y=statistic, fill=log(p.value))) +
#  coord_flip() +
#  facet_grid(statistic_name ~ new_sample_source)
#ggsave("/tmp/tt_vs_tt_sc12_sc12c.png", height=30, width=7)
