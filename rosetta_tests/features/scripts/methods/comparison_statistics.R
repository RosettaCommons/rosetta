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
	data.frame("KL Divergence"=sum(ad * log(ad / (bd+.0001)), na.rm = T))
}

# Evaluate a two sample tests between different classes of samples
# conditional on distinct groups of identifying variables.
#   The class.vars is used to identify the classes of the samples
#   The id.vars is used to split the classes in to comparison groups
comparison_statistics <- function(
	f, class.vars, id.vars, measure.vars, comp_fun){
	for(var in class.vars){
		if( !(var %in% names(f))){
			stop(paste("class.vars variable '", var, "', must be a column of f.\n\tnames(f) = c('", paste(names(f), collapse="', '"), "')", sep=""))
		}
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

	ref.class.vars <- class.vars
  names(ref.class.vars) <- paste("ref_", class.vars, sep="")
  new.class.vars <- class.vars
  names(new.class.vars) <- paste("new_", class.vars, sep="")

  ddply(f, .variables = id.vars, function(sub_f){
    ddply(sub_f, .variables = ref.class.vars, function(ref_f){
      ddply(sub_f, .variables = new.class.vars, function(new_f){
        comp_fun(ref_f[,measure.vars], new_f[,measure.vars])
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
