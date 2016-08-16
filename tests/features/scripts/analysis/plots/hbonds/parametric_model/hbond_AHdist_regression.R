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
id = "hbond_AHdist_regression",
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
  hbond.hbond_id =  geom.hbond_id AND
  hbond.struct_id = don_site.struct_id AND
  hbond.don_id = don_site.site_id AND
  hbond.struct_id = acc_site.struct_id AND
  hbond.acc_id = acc_site.site_id AND
	geom.AHdist < 2.6;"
f <- query_sample_sources(sample_sources, sele)

sele <-"
SELECT DISTINCT
	ev.don_chem_type,
	ev.acc_chem_type,
	p.degree,
	p.c_a, p.c_b, p.c_c, p.c_d, p.c_e, p.c_f, p.c_g, p.c_h, p.c_i, p.c_j, p.c_k
FROM
	hbond_evaluation_types AS ev,
	hbond_polynomial_1d AS p
WHERE
	ev.separation = 'seq_sep_other' AND
	ev.database_tag = p.database_tag AND ev.AHdist = p.name;"
polynomials <- query_sample_sources(sample_sources, sele)

## This is deprecated please use the hbond_chem_types table for the lables instead
## Order the plots better and give more descriptive labels
#f$don_chem_type <- factor(f$don_chem_type,
#	levels = c("hbdon_IMD", "hbdon_IME", "hbdon_GDE", "hbdon_GDH",
#		"hbdon_AHX", "hbdon_HXL", "hbdon_IND", "hbdon_AMO", "hbdon_CXA", "hbdon_PBA"),
#	labels = c("dIMD: h", "dIME: h", "dGDE: r", "dGDH: r",
#		"dAHX: y", "dHXL: s,t", "dIND: w", "dAMO: k", "dCXA: n,q", "dPBA: bb"))
#
## This is deprecated please use the hbond_chem_types table for the lables instead
## Order the plots better and give more descriptive labels
#f$acc_chem_type <- factor(f$acc_chem_type,
#	levels = c("hbacc_IMD", "hbacc_IME", "hbacc_AHX", "hbacc_HXL",
#		"hbacc_CXA", "hbacc_CXL", "hbacc_PBA"),
#	labels = c("aIMD: h", "aIME: h", "aAHX: y", "aHXL: s,t",
#		"aCXA: n,q", "aCXL: d,e", "aPBA: bb"))
#
## This is deprecated please use the hbond_chem_types table for the lables instead
## Order the plots better and give more descriptive labels
#polynomials$don_chem_type <- factor(polynomials$don_chem_type,
#	levels = c("hbdon_IMD", "hbdon_IME", "hbdon_GDE", "hbdon_GDH",
#		"hbdon_AHX", "hbdon_HXL", "hbdon_IND", "hbdon_AMO", "hbdon_CXA", "hbdon_PBA"),
#	labels = c("dIMD: h", "dIME: h", "dGDE: r", "dGDH: r",
#		"dAHX: y", "dHXL: s,t", "dIND: w", "dAMO: k", "dCXA: n,q", "dPBA: bb"))
#
## This is deprecated please use the hbond_chem_types table for the lables instead
## Order the plots better and give more descriptive labels
#polynomials$acc_chem_type <- factor(polynomials$acc_chem_type,
#	levels = c("hbacc_IMD", "hbacc_IME", "hbacc_AHX", "hbacc_HXL",
#		"hbacc_CXA", "hbacc_CXL", "hbacc_PBA"),
#	labels = c("aIMD: h", "aIME: h", "aAHX: y", "aHXL: s,t",
#		"aCXA: n,q", "aCXL: d,e", "aPBA: bb"))
#
get_polynomial <- function(sample_source, don_chem_type, acc_chem_type){
	df <- polynomials[
		polynomials$sample_source == sample_source &
		polynomials$don_chem_type == don_chem_type &
		polynomials$acc_chem_type == acc_chem_type,]
	if(nrow(df) != 1){
		print("ERROR: in get_polynomial")
		print(paste(don_chem_type, acc_chem_type))
		print(df)
	}
	switch(df$degree,
		p <- df[1,c("c_a")],
		p <- df[1,c("c_b","c_a")],
		p <- df[1,c("c_c","c_b","c_a")],
		p <- df[1,c("c_d","c_c","c_b","c_a",)],
		p <- df[1,c("c_e","c_d","c_c","c_b","c_a")],
		p <- df[1,c("c_f","c_e","c_d","c_c","c_b","c_a")],
		p <- df[1,c("c_g","c_f","c_e","c_d","c_c","c_b","c_a")],
		p <- df[1,c("c_h","c_g","c_f","c_e","c_d","c_c","c_b","c_a")],
		p <- df[1,c("c_i","c_h","c_g","c_f","c_e","c_d","c_c","c_b","c_a")],
		p <- df[1,c("c_j","c_i","c_h","c_g","c_f","c_e","c_d","c_c","c_b","c_a")],
		p <- df[1,c("c_k","c_j","c_i","c_h","c_g","c_f","c_e","c_d","c_c","c_b","c_a")])
	polynomial(p)
}


morse_fn <- function(x, d, a, r, m){
	return(d*(1+exp(-2*a*(x-r))-2*exp(-a*(x-r)))+m)
}


all_dens <- estimate_density_1d(
	f, c("sample_source"), "AHdist", radial_3d_normalization)

by_don_dens <- estimate_density_1d(
  f, c("sample_source", "don_chem_type"), "AHdist", radial_3d_normalization)

by_acc_dens <- estimate_density_1d(
  f, c("sample_source", "acc_chem_type"), "AHdist", radial_3d_normalization)

by_chem_dens <- estimate_density_1d(
  f, c("sample_source", "don_chem_type", "acc_chem_type"),
	"AHdist", radial_3d_normalization)


#plot_id <- "AHdist_morse_regression"
#ggplot() + theme_bw() +
#	stat_smooth(aes(x,-log(y), weight=y), data=all_dens, se=FALSE, size=1.2,
#		method=nls, formula=y ~ morse_fn(x, d, a, r, m),
#		start = list(d=7, a=1.8, r=1.97, m=-1), trace=TRUE) +
#	geom_line(aes(x=x, y=-log(y), colour=don_chem_type), data=by_don_dens) +
#	ggtitle("Hydrogen Bonds A-H Distance by Chemical Type\nnormalized for equal weight per unit distance") +
#	facet_wrap(~sample_source, ncol=1) +
#	scale_y_continuous("-log(FeatureDensity)", limits=c(-2,5), breaks=((0:15)/2-2.5))
#	scale_x_continuous(expression(paste('Acceptor -- Proton Distance (', ring(A), ')')), limits=c(1.5,2.6), breaks=c(1.6, 1.9, 2.2, 2.5))
#save_plots(self, plot_id, sample_sources, output_dir, output_formats)


param.grid <- expand.grid(list(
	d=seq(2, 6, by=3), a=seq(2, 6, by=3), r=1.7, m=seq(-1.8, -.4, by=.4)))

#
#m2 <- nls(
# -log(y) ~ morse_fn(x, d, a, r, m),
# all_dens, start=param.grid, algorithm="port", trace=TRUE, weights=y)

m.chem_type <- ddply(by_chem_dens, c("sample_source", "don_chem_type", "acc_chem_type"),
	function(df){
		print(df[1, c("don_chem_type", "acc_chem_type", "counts")])
		success <- try(
			m2 <- nls(data=df, -log(y) ~ morse_fn(x, d, a, r, m),
				 start=list(d=4, a=3, r=1.8, m=-0.8),
				 algorithm="port",
				 trace=TRUE,
				 weights=y,
				 control=nls.control(maxiter=100)))
		if(class(success)=="try-error"){return(data.frame())
		} else {
				return(data.frame(t(m2$m$getPars())))
		}
})

dens.morse_fitted <- ddply(m.chem_type, c("sample_source", "don_chem_type", "acc_chem_type"), function(v){
	df <- by_chem_dens[
		by_chem_dens$don_chem_type==v$don_chem_type & by_chem_dens$acc_chem_type==v$acc_chem_type,]
	transform(df,
		neg_log_y_fitted = morse_fn(df$x, v$d, v$a, v$r, v$m))
})


print(summary(m.chem_type))
m.chem_type_standardized <- transform(m.chem_type,
	m=-.8,
	d=1.5,
	a=2.4)

dens.morse_fitted_standardized <- ddply(m.chem_type_standardized, c("sample_source", "don_chem_type", "acc_chem_type"), function(v){
	df <- by_chem_dens[
		by_chem_dens$don_chem_type==v$don_chem_type & by_chem_dens$acc_chem_type==v$acc_chem_type,]
	transform(df,
		neg_log_y_fitted = morse_fn(df$x, v$d, v$a, v$r, v$m))
})



plot_id = "morse_potential_parameter_plotmatrix"
plotmatrix(data=m.chem_type[,c("d", "a", "r", "m")], aes(colour=don_chem_type)) +
	theme_bw() +
	ggtitle("Parameters for Morse Potential fit to Hydrogen Bond AHdist Densities\nGrouped by Donor and Acceptor Chemical Types\nd*(1 + exp(-2*a*(x - r)) - 2*exp(-a*(x - r))) + m") +
	coord_equal(ratio=1)
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

print(summary(dens.morse_fitted))
plot_id <- "AHdist_morse_regression_residuals_by_chem_type"
a_ply(sample_sources, 1, function(sample_source){
	ggplot(dens.morse_fitted) + theme_bw() +
		geom_line(aes(x=x, y=-log(y) - neg_log_y_fitted)) +
		geom_indicator(aes(indicator=counts)) +
		facet_grid(don_chem_type ~ acc_chem_type) +
		ggtitle(paste("Hydrogen Bonds A-H Distance by Chemical Type\nnormalized for equal weight per unit distance  Sample Source:", sample_source$sample_source[1])) +
		labs(x=expression(paste('Acceptor -- Proton Distance (', ring(A), ')')),
		     y="log(FeatureDensity + 1)") +
		scale_y_continuous(limits=c(-5,10), breaks=c(0,10)) +
		scale_x_continuous(limits=c(1.5,2.6), breaks=c(1.6, 1.9, 2.2, 2.5))
	save_plots(self, plot_id, sample_source, output_dir, output_formats)
})

dens.poly_fitted <- ddply(by_chem_dens,
	c("sample_source", "don_chem_type", "acc_chem_type"), function(df){
		data.frame(
			x=df$x,
			y=predict(get_polynomial(
				df$sample_source[1],
				as.character(df$don_chem_type[1]),
				as.character(df$acc_chem_type[1])),
				df$x))
	})

## This is deprecated please use the hbond_chem_types table for the lables instead
## Order the plots better and give more descriptive labels
#by_chem_dens$don_chem_type <- factor(by_chem_dens$don_chem_type,
#	levels = c("hbdon_IMD", "hbdon_IME", "hbdon_GDE", "hbdon_GDH",
#		"hbdon_AHX", "hbdon_HXL", "hbdon_IND", "hbdon_AMO", "hbdon_CXA", "hbdon_PBA"),
#	labels = c("dIMD: h", "dIME: h", "dGDE: r", "dGDH: r",
#		"dAHX: y", "dHXL: s,t", "dIND: w", "dAMO: k", "dCXA: n,q", "dPBA: bb"))
#
## This is deprecated please use the hbond_chem_types table for the lables instead
## Order the plots better and give more descriptive labels
#by_chem_dens$acc_chem_type <- factor(by_chem_dens$acc_chem_type,
#	levels = c("hbacc_IMD", "hbacc_IME", "hbacc_AHX", "hbacc_HXL",
#		"hbacc_CXA", "hbacc_CXL", "hbacc_PBA"),
#	labels = c("aIMD: h", "aIME: h", "aAHX: y", "aHXL: s,t",
#		"aCXA: n,q", "aCXL: d,e", "aPBA: bb"))
#
kT <- .4
plot_id <- "AHdist_morse_regression_by_chem_type"
a_ply(sample_sources, 1, function(sample_source){
	ss_id <- sample_source$sample_source[1]
	ggplot(data=subset(by_chem_dens, sample_source==ss_id)) + theme_bw() +
		geom_line(aes(x, neg_log_y_fitted*kT, group=sample_source), color="gray",
			data=subset(dens.morse_fitted, sample_source==ss_id)) +
		geom_line(aes(x, neg_log_y_fitted, group=sample_source), color="blue",
			data=subset(dens.morse_fitted_standardized, sample_source==ss_id)) +
		geom_line(aes(x,y), data=dens.poly_fitted, color="purple") +
		geom_line(aes(x, -log(y)*kT)) +
		geom_indicator(aes(indicator=counts)) +
		facet_grid(don_chem_type ~ acc_chem_type) +
		ggtitle(paste("Hydrogen Bonds A-H Distance by Chemical Type Fitted with Morse Function\nnormalized for equal weight per unit distance  Sample Source:", sample_source$sample_source[1])) +
		labs(x=expression(paste('Acceptor -- Proton Distance (', ring(A), ')')),
		     y="log(FeatureDensity)") +
		scale_y_continuous(limits=c(-.8,1.3)) +
	#	scale_y_continuous(limits=c(-2,5), breaks=c(-2,0,2,4)) +
		scale_x_continuous(limits=c(1.5,2.6), breaks=c(1.6, 1.9, 2.2, 2.5))
	save_plots(self, plot_id, sample_source, output_dir, output_formats)
})

#l_ply(by_chem_dens, aes("don_chem_type", "acc_chem_type"), function(dens){
#	plot_id <- paste("AHdist_morse_polynomial_regression", dens$don_chem_type[1], dens$acc_chem_type[1], sep="_")
#	ggplot(dens) + theme_bw() +
#		geom_line(aes(x, y, colour=sample_source)) +
#		geom_indicator(aes(indicator=counts, colour=sample_source)) +
#		ggtitle(paste("HBond AHdist Fitted with Morse Function and Polynomial\nnormalized for equal weight per unit distance Donor:", don_chem_type, " Acceptor:", acc_chem_type, sep="")) +
#		scale_x_continuous(expression(paste('Acceptor -- Proton Distance (', ring(A), ')')), limits=c(1.5,2.6), breaks=c(1.6, 1.9, 2.2, 2.5)) +
#		scale_y_continuous("log(FeatureDensity)")
#	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
#})


})) # end FeaturesAnalysis
