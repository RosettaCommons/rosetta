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
id = "hbond_cosAHD_regression",
author = "Matthew O'Meara",
brief_description = "",
long_description = "
The cosAHD angle in a hydrogen bond is the cosine of the angle at
the hydrogen.  Usually the hydrogen is placed close to the line
connecting the acceptor and donor atoms. To model this, the beta
function is used. The beta function goes from 0 to 1 where the first
parameter conceptually controls the angle on the left and the second
parameter controls the angle on the right.

This script generates 4 plots grouping the hydrogen bonds different ways:
  1) All the hydrogen bonds together
  2) Conditional on the donor chemical type
  3) Conditional on the acceptor chemical type
  4) Conditional on both the donor and acceptor chemical types

Currently this is only setup for a single sample source at a time.",

feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <-"
SELECT
  geom.cosAHD,
  acc_site.HBChemType AS acc_chem_type,
  don_site.HBChemType AS don_chem_type
FROM
  hbond_geom_coords AS geom,
  hbonds AS hbond,
  hbond_sites AS don_site,
  hbond_sites AS acc_site
WHERE
  hbond.struct_id = geom.struct_id AND
  hbond.hbond_id =  geom.hbond_id AND
  hbond.struct_id = don_site.struct_id AND
  hbond.don_id = don_site.site_id AND
  hbond.struct_id = acc_site.struct_id AND
  hbond.acc_id = acc_site.site_id;"
f <- query_sample_sources(sample_sources, sele)

# This is deprecated please use the hbond_chem_types table for the lables instead
# Order the plots better and give more descriptive labels
f$don_chem_type <- factor(f$don_chem_type,
	levels = c("hbdon_IMD", "hbdon_IME", "hbdon_GDE", "hbdon_GDH",
		"hbdon_AHX", "hbdon_HXL", "hbdon_IND", "hbdon_AMO", "hbdon_CXA", "hbdon_PBA"),
	labels = c("dIMD: h", "dIME: h", "dGDE: r", "dGDH: r",
		"dAHX: y", "dHXL: s,t", "dIND: w", "dAMO: k", "dCXA: n,q", "dPBA: bb"))

# This is deprecated please use the hbond_chem_types table for the lables instead
# Order the plots better and give more descriptive labels
f$acc_chem_type <- factor(f$acc_chem_type,
	levels = c("hbacc_IMD", "hbacc_IME", "hbacc_AHX", "hbacc_HXL",
		"hbacc_CXA", "hbacc_CXL", "hbacc_PBA"),
	labels = c("aIMD: h", "aIME: h", "aAHX: y", "aHXL: s,t",
		"aCXA: n,q", "aCXL: d,e", "aPBA: bb"))


all_dens <- estimate_density_1d(
  f, c("sample_source"), "cosAHD", histogram=TRUE)

all_m <- nls(y ~ c*dbeta(x, a, b), all_dens, start=list(a=11, b=1.2, c=.04), algorithm="port", trace=TRUE, weight=y)
all_dens$fitted <- predict(all_m)
all_dens$param_string <- paste(c("a:", "b:", "c:"), round(all_m$m$getPars(),3), collapse=", ")

plot_id <- "hbond_cosAHD_regression_all"
ggplot(data=all_dens, aes(x=x)) + theme_bw() +
  geom_line(aes(y=y), size=1.4) +
  geom_line(aes(y=fitted), colour="blue", size=1.4) +
  geom_indicator(aes(indicator=param_string)) +
  ggtitle("Hydrogen Bonds AHD Angle Fit with Beta Function\n(normalized for equal volume per unit distance)") +
  labs(x=expression(paste('Acceptor -- Hydrogen -- Donor (degrees)')),
       y="FeatureDensity")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)



don_dens <- estimate_density_1d(
  f, c("sample_source", "don_chem_type"), "cosAHD", histogram=TRUE)

don_dens <- ddply(don_dens, .variables=c("sample_source", "don_chem_type"), function(df){
  cat("sample_source:", df$sample_source[1], "don_chem_type:", as.character(df$don_chem_type[1]), "\n")
  success <- try({
    params.grid <- expand.grid(a=c(2), b=c(.7), c=c(.5))
    m <- nls(y ~ c*dbeta(x, a, b), df, start=params.grid, algorithm="port", trace=TRUE, control=nls.control(maxiter=100), weight=y)})
  if(class(success)=="try-error"){
    df$fitted <- NA
  } else {
    df$fitted <- predict(m)
    df$param_string <- paste(c("a:", "b:", "c:"), round(m$m$getPars(),3), collapse=", ")
  }
  df
})


plot_id <- "hbond_cosAHD_regression_don"
ggplot(data=don_dens, aes(x=x)) + theme_bw() + facet_wrap( ~ don_chem_type ) +
  geom_line(aes(x,y)) +
  geom_line(aes(x,fitted), colour="blue") +
  geom_indicator(aes(indicator=param_string)) +
  ggtitle("Hydrogen Bonds AHD Angle by Don Chemical Type fit with Beta Function\n(normalized for equal volume per unit distance)") +
  labs(x=expression(paste('Acceptor -- Hydrogen -- Donor (degrees)')),
       y="FeatureDensity")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)






acc_dens <- estimate_density_1d(
  f, c("sample_source", "acc_chem_type"), "cosAHD", histogram=TRUE)

acc_dens <- ddply(acc_dens, .variables=c("sample_source", "acc_chem_type"), function(df){
  cat("sample_source:", df$sample_source[1], "acc_chem_type:", as.character(df$acc_chem_type[1]), "\n")
  success <- try({
    params.grid <- expand.grid(a=c(2), b=c(.7), c=c(.5))
    m <- nls(y ~ c*dbeta(x, a, b), df, start=params.grid, algorithm="port", trace=TRUE, control=nls.control(maxiter=100), weight=y)})
  if(class(success)=="try-error"){
    df$fitted <- NA
  } else {
    df$fitted <- predict(m)
    df$param_string <- paste(c("a:", "b:", "c:"), round(m$m$getPars(),3), collapse=", ")
  }
  df
})


plot_id <- "hbond_cosAHD_regression_acc"
ggplot(data=acc_dens, aes(x=x)) + theme_bw() + facet_wrap( ~ acc_chem_type ) +
  geom_line(aes(x,y)) +
  geom_line(aes(x,fitted), colour="blue") +
  geom_indicator(aes(indicator=param_string)) +
  ggtitle("Hydrogen Bonds AHD Angle by Acceptor Chemical Type fit with Beta Function\n(normalized for equal volume per unit distance)") +
  labs(x=expression(paste('Acceptor -- Hydrogen -- Donor (degrees)')),
       y="FeatureDensity")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)




each_dens <- estimate_density_1d(
  f, c("sample_source", "don_chem_type", "acc_chem_type"), "cosAHD", histogram=TRUE)

each_dens <- ddply(each_dens, .variables=c("sample_source", "don_chem_type", "acc_chem_type"), function(df){
  cat("sample_source:", df$sample_source[1], "acc_chem_type:", as.character(df$acc_chem_type[1]), "\n")
  cat("sample_source:", df$sample_source[1], "don_chem_type:", as.character(df$don_chem_type[1]), "\n")
  success <- try({
    params.grid <- expand.grid(a=c(2), b=c(.7), c=c(.5))
#			a=seq(2, 10, length.out=3),
#			b=seq(.01, 2, length.out=3),
#		  c=seq(.01,.5, length.out=3))
    m <- nls(y ~ c*dbeta(x, a, b), df, start=params.grid, algorithm="port", trace=TRUE, control=nls.control(maxiter=100), weight=sqrt(y))})
  if(class(success)=="try-error"){
    df$fitted <- NA
  } else {
    df$fitted <- predict(m)
    df$param_string <- paste(c("a:", "b:", "c:"), round(m$m$getPars(),2), collapse=", ")
  }
  df
})


plot_id <- "hbond_cosAHD_regression_each"
ggplot(data=each_dens, aes(x=x)) + theme_bw() + facet_grid( don_chem_type ~ acc_chem_type ) +
  geom_line(aes(x,y)) +
  geom_line(aes(x,fitted), colour="blue") +
  geom_indicator(aes(indicator=param_string)) +
  ggtitle("Hydrogen Bonds AHD Angle by Chemical Type fit with Beta Function\n(normalized for equal volume per unit distance)") +
  labs(x=expression(paste('Acceptor -- Hydrogen -- Donor (degrees)')),
       y="FeatureDensity")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
