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
id = "OHdonor",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

source("scripts/methods/polynomial_methods.R")


extract_features <- function(sample_sources, don_chem_type, acc_chem_type, xmin, xmax){
	sele <-paste("
SELECT
	geom.AHdist
FROM
	hbond_geom_coords AS geom,
	hbonds AS hb,
	hbond_sites AS don, hbond_sites AS acc,
	hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb
WHERE
	geom.struct_id = hb.struct_id AND geom.hbond_id = hb.hbond_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	don.HBChemType = '", don_chem_type, "' AND
	acc.HBChemType = '", acc_chem_type, "' AND
	", xmin, " <= geom.AHdist AND
	geom.AHdist <= ", xmax, " AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30;", sep="")
	query_sample_sources(sample_sources, sele)
}

estimate_density <- function(f){
	estimate_density_1d(
		data = f,
		ids = c("sample_source"),
		variable = "AHdist",
		weight_fun = radial_3d_normalization)
}


extract_model <- function(sample_sources, don_chem_type, acc_chem_type){
	sele <-paste("
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
		ev.don_chem_type = '", don_chem_type, "' AND
		ev.acc_chem_type = '", acc_chem_type, "' AND
		ev.separation = 'seq_sep_other' AND
		ev.database_tag = p.database_tag AND ev.AHdist = p.name;", sep="")
	query_sample_sources(sample_sources, sele)

}

extract_evaluate_models <- function(
	sample_sources,
	don_chem_type,
	acc_chem_type,
	xmin,
	xmax,
	sample_source="A-H Distance Term",
	temperature = 10,
	normalization = .03
) {
	m_polys <- extract_model(
		subset(sample_sources, !is.na(model)),
		don_chem_type, acc_chem_type)

	ddply(m_polys, .(sample_source), function(m){
		poly <- get_1d_polynomial(m, don_chem_type, acc_chem_type)

		x <- seq(xmin, xmax, length.out=300)

		boltzmann_transform <- function(x, e, xmin, xmax) {
			ifelse(x >= xmin & x <= xmax, normalization * exp(-e*temperature), NA)
		}

		data.frame(
			x = x,
			y = boltzmann_transform(x, predict(poly$poly, x), poly$xmin, poly$xmax),
			sample_source = sample_source,
			model = sample_sources[sample_sources$sample_source == m$sample_source, "model",],
			counts = NA)
	})
}

generate_plot <- function(
	dens, model,
	don_chem_type, acc_chem_type,
	sample_sources, output_dir, output_formats) {

	plot_id <- paste("OHdonor_AHdist", don_chem_type, acc_chem_type, sep="_")

	ref_dens <- transform(subset(dens, reference), sample_source=NULL)
	new_dens <- subset(dens, !reference)

#	facet_labels <- data.frame(
#		sample_source = sample_sources$sample_source,
#		label=toupper(letters[1:nrow(sample_sources)]))
#
#	facet_labels <- rbind(
#		data.frame(sample_source = unique(model$sample_source)),
#		data.frame(sample_source = unique(new_dens$sample_source)))
#	facet_labels$label <- toupper(letters[1:nrow(facet_labels)])

	#hack!
	facet_labels <- data.frame(
		sample_source = c("A-H Distance Term", "Relaxed Native Score12", "Relaxed Native NewHB", "Relaxed Native NewHB LJcorr"),
		label = c("A", "B", "C", "D"))

	p <- ggplot() +
		geom_line(data=ref_dens, aes(x=x, y=y), colour="black", size=2) +
		geom_line(data=model, aes(x=x, y=y, colour=model), size=1) +
		geom_line(data=new_dens, aes(x=x, y=y, colour=model), size=2) +
		geom_indicator(data=facet_labels, aes(indicator=label),
			group=1, color="black", xpos=.03, ypos=.97, size=12) +
		facet_wrap(~sample_source, ncol=2) +
		scale_y_continuous("Feature Density") +
		scale_x_continuous(
			expression(paste('Acceptor -- Hydrogen Distance (', ring(A), ')'))) +
		scale_colour_grey("HBond Potential", start=.4, end=.7) +
		theme_bw() +
		theme(legend.position = c(.36, .85))

	save_plots(self, plot_id=plot_id, sample_sources, output_dir, output_formats)
}


##########################

xmin <- 1.5; xmax <- 2.2
don_chem_type <- "hbdon_HXL"; acc_chem_type <- "hbacc_PBA"

f <- extract_features(
	sample_sources,
	don_chem_type, acc_chem_type,
	xmin, xmax)

#counts, e.g. to add to caption
print(ddply(f, .(sample_source), nrow))

dens <- estimate_density(f)
dens <- merge(dens, sample_sources[, c("sample_source", "reference", "model")])

model <- extract_evaluate_models(
	sample_sources,
	don_chem_type, acc_chem_type,
	xmin, xmax)

generate_plot(
	dens, model,
	don_chem_type, acc_chem_type,
	sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
