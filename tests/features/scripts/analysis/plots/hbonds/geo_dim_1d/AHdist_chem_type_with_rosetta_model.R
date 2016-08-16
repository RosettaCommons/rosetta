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
id = "AHdist_chem_type_with_rosetta_model",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures", "HBondParameterFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

source("scripts/analysis/plots/hbonds/hbond_geo_dim_scales.R")
source("scripts/methods/polynomial_methods.R")


sele <-"
SELECT
	geom.AHdist,
	don.HBChemType AS don_chem_type, acc.HBChemType AS acc_chem_type
FROM
	hbonds AS hb,
	hbond_geom_coords AS geom,
	hbond_sites AS don, hbond_sites AS acc
WHERE
	geom.struct_id = hb.struct_id AND geom.hbond_id = hb.hbond_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id;"
f <- query_sample_sources(sample_sources, sele)

sele <-"
SELECT DISTINCT
	ev.don_chem_type,
	ev.acc_chem_type,
	ev.separation AS seq_sep,
	p.dimension,
	p.xmin, p.xmax,
	p.degree,
	p.c_a, p.c_b, p.c_c, p.c_d, p.c_e, p.c_f, p.c_g, p.c_h, p.c_i, p.c_j, p.c_k
FROM
	hbond_evaluation_types AS ev,
	hbond_polynomial_1d AS p
WHERE
	ev.separation = 'seq_sep_other' AND
	ev.database_tag = p.database_tag AND ev.AHdist = p.name;"
polynomials <- query_sample_sources(sample_sources, sele)


xmin <- min(f$AHdist, polynomials$xmin)
xmax <- max(f$AHdist, polynomials$xmax)

dens <- estimate_density_1d(
	f, c("sample_source", "acc_chem_type", "don_chem_type"),
	"AHdist", weight_fun = radial_3d_normalization)


dens$y_energy <- -log(dens$y)*.35

dens.model <-
  expand.grid(
    sample_source = sample_sources$sample_source[2],
    acc_chem_type = levels(f$acc_chem_type),
    don_chem_type = levels(f$don_chem_type))

dens.model <- ddply(dens.model,
  .(sample_source, acc_chem_type, don_chem_type), function(df){

  sample_source <- as.character(df$sample_source[1])
  acc_chem_type <- as.character(df$acc_chem_type[1])
  don_chem_type <- as.character(df$don_chem_type[1])

  poly <- get_1d_polynomial(
    polynomials, don_chem_type, acc_chem_type)

  x <- seq(xmin, xmax, length.out=100)

	normalization <- .03
	temperature <- 10
	z <- data.frame(
		x=x,
		y_energy=predict(poly$poly, x),
		sample_source=factor("Rosetta Model"))
	z$y_energy <- ifelse(z$x >= poly$xmin & z$x <= poly$xmax, z$y, NA)
	z$y <- normalization * exp(-z$y_energy * temperature)
	z
})

dens.model$counts <- NA
dens <- rbind(dens, dens.model)


# Order the plots better and give more descriptive labels
dens$don_chem_type_name <- factor(dens$don_chem_type,
	levels = c("hbdon_IMD", "hbdon_IME", "hbdon_GDE", "hbdon_GDH",
		"hbdon_AHX", "hbdon_HXL", "hbdon_IND", "hbdon_AMO", "hbdon_CXA", "hbdon_PBA"),
	labels = c("dIMD: h", "dIME: h", "dGDE: r", "dGDH: r",
		"dAHX: y", "dHXL: s,t", "dIND: w", "dAMO: k", "dCXA: n,q", "dPBA: bb"))

# Order the plots better and give more descriptive labels
dens$acc_chem_type_name <- factor(dens$acc_chem_type,
	levels = c("hbacc_IMD", "hbacc_IME", "hbacc_AHX", "hbacc_HXL",
		"hbacc_CXA", "hbacc_CXL", "hbacc_PBA"),
	labels = c("aIMD: h", "aIME: h", "aAHX: y", "aHXL: s,t",
		"aCXA: n,q", "aCXL: d,e", "aPBA: bb"))

dens <- dens[!is.na(dens$acc_chem_type_name) & !is.na(dens$don_chem_type_name),]

plot_id <- "hbond_AHdist_chem_type_with_parameters_energy"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x, y_energy, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	facet_grid(don_chem_type_name ~ acc_chem_type_name) +
	ggtitle("HBond A-H Distance by Chemical Type, B-Factor < 30\nnormalized for equal weight per unit distance in density estimation") +
	scale_y_continuous("Energy (arbitrary units)", limits=c(-.6,1.1), breaks=c(-.5,0, .5)) +
	scale_x_continuous(expression(paste('Acceptor -- Proton Distance (', ring(A), ')')), limits=c(.5,3), breaks=c(1, 1.5, 2, 2.5))

if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


plot_id <- "hbond_AHdist_chem_type_with_parameters_probability"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x, y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	facet_grid(don_chem_type_name ~ acc_chem_type_name) +
	ggtitle("HBond A-H Distance by Chemical Type, B-Factor < 30\nnormalized for equal weight per unit distance in density estimation") +
	scale_y_continuous("Feature Density") +
	scale_x_continuous(expression(paste('Acceptor -- Proton Distance (', ring(A), ')')), limits=c(.5,3), breaks=c(1, 1.5, 2, 2.5))

if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
