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
id = "hbond_regression",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <-"
SELECT
  geom.AHdist,
  acc_site.HBChemType AS acc_chem_type,
  don_site.HBChemType AS don_chem_type
FROM
  hbond_geom_coords AS geom,
  hbonds AS hbond,
  hbond_sites AS don_site,
  hbond_sites AS acc_site
WHERE
  hbond.struct_id = geom.struct_id AND
  hbond.hbond_id = geom.hbond_id AND
  hbond.struct_id = don_site.struct_id AND
  hbond.don_id = don_site.site_id AND
  hbond.struct_id = acc_site.struct_id AND
  hbond.acc_id = acc_site.site_id AND
  geom.AHdist < 2.6;";

all_geom <- query_sample_sources(sample_sources, sele)

morse_fn <- function(x, D_a, a, r_0, min_e){
	return(D_a*(1+exp(-2*a*(x-r_0))-2*exp(-a*(x-r_0)))+min_e)
}


######## Join all types of hydrogen bonds together in a single distribution
plot_id <- "HBond_AHdist_morse_regression_unified"

all_dens <- estimate_density_1d(
  data = all_geom,
  ids = c("sample_source"),
  variable = "AHdist",
  weight_fun = radial_3d_normalization,
  adjust=1
)
all_dens$neg_log_y = -log(all_dens$y)

all_dens <- ddply(all_dens, c("sample_source"), function(df){
  m1 <- nls(neg_log_y ~ morse_fn(x, D_a, a, r_0, min_e), df, start = list(D_a=7, a=1.8, r_0=1.8, min_e=-1), algorithm="port", trace=TRUE, weights=y)
	print(m1)
  df$fitted <- predict(m1,df$x)
	df
})

dens0.5 <- estimate_density_1d(
  data = all_geom,
  ids = c("sample_source"),
  variable = "AHdist",
  weight_fun = radial_3d_normalization,
  adjust=.5
)
dens0.5$neg_log_y = -log(dens0.5$y)


ggplot(all_dens) + theme_bw() +
  geom_line(aes(x=x, y=fitted, colour=sample_source), size=1.3) +
  geom_line(aes(x=x, y=neg_log_y, colour=sample_source)) +
  geom_line(data=dens0.5, aes(x=x, y=neg_log_y, colour=sample_source)) +
  ggtitle("Hydrogen Bonds A-H Distance Fit with Morse Potential\nnormalized for equal weight per unit distance") +
  labs(x=expression(paste('Acceptor -- Proton Distance (', ring(A), ')')),
       y="-log(FeatureDensity)") +
  scale_y_continuous(limits=c(-2,5), breaks=((0:30)/4-2.5))
  scale_x_continuous(limits=c(1.5,3), breaks=c(1.6, 1.9, 2.2, 2.6))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)




#####  Fit morse potential to each chemical type ####

dens <- estimate_density_1d(
  data = all_geom,
  ids = c("sample_source", "don_chem_type", "acc_chem_type"),
  variable = "AHdist",
  weight_fun = radial_3d_normalization)
dens$neg_log_y = -log(dens$y)

dens <- ddply(dens, .variables=c("don_chem_type", "acc_chem_type"),
	function(df){
	  print(df$don_chem_type[1])
		print(df$acc_chem_type[1])
		print(df$count[1])
		success <- try(
		m <- nls(
			 neg_log_y ~ morse_fn(x, D_a, a, r_0, min_e),
			 df,
			 start=list(D_a=4, a=3, r_0=1.8, min_e=-0.8),
			 algorithm="port",
			 trace=TRUE,
			 weights=y,
			 control=nls.control(maxiter=100)))
		if(class(success)=="try-error"){
				return(data.frame())
		} else {
				return(data.frame(t(m$m$getPars())))
		}
		df$fitted <- predict(m,df$x)
		df
})

print(summary(dens))

plot_id <- "HBond_AHdist_regression_by_chem_type"
ggplot(dens) + theme_bw() +
  geom_line(aes(x=x, y=neg_log_y, colour=sample_source)) +
  geom_indicator(aes(indicator=count, colour=sample_source, group=sample_source)) +
  facet_grid(don_chem_type ~ acc_chem_type) +
  ggtitle("Hydrogen Bonds A-H Distance by Chemical Type Fitted with Morse Function\nnormalized for equal weight per unit distance") +
  labs(x=expression(paste('Acceptor -- Proton Distance (', ring(A), ')')),
       y="log(FeatureDensity)") +
  scale_y_continuous(limits=c(-2,5), breaks=c(-2,0,2,4)) +
  scale_x_continuous(limits=c(1.5,3), breaks=c(1.6, 1.9, 2.2, 2.6)) +
save_plots(self, plot_id, sample_sources, output_dir, output_formats)




})) # end FeaturesAnalysis
