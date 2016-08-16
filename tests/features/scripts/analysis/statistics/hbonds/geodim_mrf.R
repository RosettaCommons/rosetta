# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

check_setup()
feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "geodim_mrf",
author = "Matthew O'Meara",
brief_description = "",
long_description = "
compute a markov random field a collection of geometric features

This is quite rough, and is a work in progress
",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


# number of folds to do cross validation
set.seed(12345) # for reproducibility

sele <-"
SELECT
	geom.AHdist, geom.cosBAH, geom.cosAHD, geom.chi,
	don.HBChemType AS don_chem_type, acc.HBChemType AS acc_chem_type
FROM
	hbonds AS hb,
	hbond_geom_coords AS geom,
	hbond_sites AS don, hbond_sites AS acc,
	hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb
WHERE
	geom.struct_id = hb.struct_id AND geom.hbond_id = hb.hbond_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30
LIMIT 400;"
f <- query_sample_sources(sample_sources, sele)

discretize_geometric_bins <- function(f, n){
	# number of bins for each geometric dimension
	transform(f,
		AHdist_bin = cut(AHdist, n),
		cosBAH_bin = cut(cosBAH, n),
		cosAHD_bin = cut(cosAHD, n),
		chi_bin = cut(chi, n))
}

assign_cross_validation_folds <- function(f, k){
	f <- f[sample(nrow(f), nrow(f)),]
	f$fold <- seq(nrow(f)) %% k
	f
}


compute_factors <- function(f){
	total = nrow(f)
  list(
    ddply(f, .(AHdist_bin, cosBAH_bin), function(df){data.frame(p_raw=nrow(df)/total)}),
    ddply(f, .(AHdist_bin, cosAHD_bin), function(df){data.frame(p_raw=nrow(df)/total)}),
    ddply(f, .(AHdist_bin, chi_bin   ), function(df){data.frame(p_raw=nrow(df)/total)}),
    ddply(f, .(cosBAH_bin, cosAHD_bin), function(df){data.frame(p_raw=nrow(df)/total)}),
    ddply(f, .(cosBAH_bin, chi_bin   ), function(df){data.frame(p_raw=nrow(df)/total)}),
    ddply(f, .(don_chem_type, AHdist_bin), function(df){data.frame(p_raw=nrow(df)/total)}),
    ddply(f, .(don_chem_type, cosBAH_bin), function(df){data.frame(p_raw=nrow(df)/total)}),
    ddply(f, .(don_chem_type, cosAHD_bin), function(df){data.frame(p_raw=nrow(df)/total)}),
    ddply(f, .(don_chem_type, chi_bin   ), function(df){data.frame(p_raw=nrow(df)/total)}),
    ddply(f, .(don_chem_type, acc_chem_type), function(df){data.frame(p_raw=nrow(df)/total)}),
    ddply(f, .(acc_chem_type, AHdist_bin), function(df){data.frame(p_raw=nrow(df)/total)}),
    ddply(f, .(acc_chem_type, cosBAH_bin), function(df){data.frame(p_raw=nrow(df)/total)}),
    ddply(f, .(acc_chem_type, cosAHD_bin), function(df){data.frame(p_raw=nrow(df)/total)}),
    ddply(f, .(acc_chem_type, chi_bin   ), function(df){data.frame(p_raw=nrow(df)/total)}))
}


	compute_p_raw <- function(factors, AHdist_bin, cosBAH_bin, cosAHD_bin, chi_bin, don_chem_type, acc_chem_type){
		10^-15 +
	  	factors[[1 ]][factors[[1 ]]$AHdist_bin    == AHdist_bin    & factors[[1 ]]$cosBAH_bin    == cosBAH_bin   , "p_raw"]*
	  	factors[[2 ]][factors[[2 ]]$AHdist_bin    == AHdist_bin    & factors[[2 ]]$cosAHD_bin    == cosAHD_bin   , "p_raw"]*
	  	factors[[3 ]][factors[[3 ]]$AHdist_bin    == AHdist_bin    & factors[[3 ]]$chi_bin       == chi_bin      , "p_raw"]*
	  	factors[[4 ]][factors[[4 ]]$cosBAH_bin    == cosBAH_bin    & factors[[4 ]]$cosAHD_bin    == cosAHD_bin   , "p_raw"]*
	  	factors[[5 ]][factors[[5 ]]$cosBAH_bin    == cosBAH_bin    & factors[[5 ]]$chi_bin       == chi_bin      , "p_raw"]*
	  	factors[[6 ]][factors[[6 ]]$don_chem_type == don_chem_type & factors[[6 ]]$AHdist_bin    == AHdist_bin   , "p_raw"]*
	  	factors[[7 ]][factors[[7 ]]$don_chem_type == don_chem_type & factors[[7 ]]$cosBAH_bin    == cosBAH_bin   , "p_raw"]*
	  	factors[[8 ]][factors[[8 ]]$don_chem_type == don_chem_type & factors[[8 ]]$cosAHD_bin    == cosAHD_bin   , "p_raw"]*
	  	factors[[9 ]][factors[[9 ]]$don_chem_type == don_chem_type & factors[[9 ]]$chi_bin       == chi_bin      , "p_raw"]*
	  	factors[[10]][factors[[10]]$don_chem_type == don_chem_type & factors[[10]]$acc_chem_type == acc_chem_type, "p_raw"]*
	  	factors[[11]][factors[[11]]$acc_chem_type == acc_chem_type & factors[[11]]$AHdist_bin    == AHdist_bin   , "p_raw"]*
	  	factors[[12]][factors[[12]]$acc_chem_type == acc_chem_type & factors[[12]]$cosBAH_bin    == cosBAH_bin   , "p_raw"]*
	  	factors[[13]][factors[[13]]$acc_chem_type == acc_chem_type & factors[[13]]$cosAHD_bin    == cosAHD_bin   , "p_raw"]
	}

train_model_brute_force <- function(factors, f, debug=F){
	model <- expand.grid(
	  AHdist_bin = levels(f$AHdist_bin),
	  cosBAH_bin = levels(f$cosBAH_bin),
	  cosAHD_bin = levels(f$cosAHD_bin),
	  chi_bin = levels(f$chi_bin),
		don_chem_type = levels(f$don_chem_type),
		acc_chem_type = levels(f$acc_chem_type))


	model <- adply(model,1,
	  function(row){
			data.frame(
				p_raw=compute_p_raw(
					factors,
					row$AHdist,
					row$cosBAH,
					row$cosAHD,
					row$chi,
					row$don_chem_type,
					row$acc_chem_type))})

	# compute normalizing constant
	Z = sum(model$p_raw)
	model$p <- model$p_raw/Z

	#check normalization, should be 1
	#print(sum(d$p))

	model
}

compute_log_likelihood <- function(f, model){
	e <- merge(f, model)
	sum(log(e$p))
}

train_and_test <- function(f, fold){
	print(paste("begin cross validation", fold))
	f_train <- f[f$fold != fold,]
	f_test <- f[f$fold == fold,]

	# train set clique potentials based on the training data
	print("compute_factors...")
	factors <- compute_factors(f_train)

	# estimate raw probabilities by brute force
	print("train_model...")
	model <- train_model_brute_force(factors, f_train)

	print("compute_log_likelihood...")
	results <- data.frame(
		fold = fold,
		ll_train = compute_log_likelihood(f_train, model),
		ll_test = compute_log_likelihood(f_test, model))

	print("results:")
	print(results)
	results
}



run <- function(nbins,k){

	print(paste("begin run with", nbins, "bins and", k, "fold cross validation."))
	f <- discretize_geometric_bins(f, nbins)
	f <- assign_cross_validation_folds(f, k)
	ldply(0:(k-1), function(fold){
		data.frame(
			nbins=nbins,
			train_and_test(f, fold))

	})

}



#
##Example
#phi1 = data.frame(
#  a=c(0,0,1,1), b=c(0,1,0,1), value=c(30, 5, 1, 10))
#
#phi2 = data.frame(
#  b=c(0,0,1,1), c=c(0,1,0,1), value=c(100, 1, 1, 100))
#
#phi3 = data.frame(
#  c=c(0,0,1,1), d=c(0,1,0,1), value=c(1, 100, 100, 1))
#
#phi4 = data.frame(
#  d=c(0,0,1,1), a=c(0,1,0,1), value=c(100, 1, 1, 100))
#
#
#P_unnormalized <- function(a,b,c,d){
#  phi1[phi1$a == a & phi1$b == b, "value"] *
#    phi2[phi2$b == b & phi2$c == c, "value"] *
#    phi3[phi3$c == c & phi3$d == d, "value"] *
#    phi4[phi4$d == d & phi4$a == a, "value"]
#}
#
#d <- expand.grid(a=c(0,1), b=c(0,1), c=c(0,1), d=c(0,1))
#
##d <- transform(d, p_raw = P_unnormalized(a,b,c,d))
#d <- adply(d,1,
#  function(row){ data.frame(p_raw=P_unnormalized(row$a, row$b, row$c, row$d))})
#
#Z = sum(d$p_raw)
#
#d$p <- d$p_raw/Z



})) # end FeaturesAnalysis
