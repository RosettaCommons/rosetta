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
id = "cosAHD_chem_type_with_rosetta_model",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

source("scripts/analysis/plots/hbonds/hbond_geo_dim_scales.R")
source("scripts/methods/polynomial_methods.R")


sele <-"
SELECT
	geom.cosAHD,
	don.HBChemType AS don_chem_type, acc.HBChemType AS acc_chem_type
FROM
	hbonds AS hb,
	hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb,
	hbond_geom_coords AS geom,
	hbond_sites AS don, hbond_sites AS acc
WHERE
	geom.struct_id = hb.struct_id AND geom.hbond_id = hb.hbond_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	don_pdb.struct_id = don.struct_id AND don_pdb.site_id = don.site_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = acc.struct_id AND acc_pdb.site_id = acc.site_id AND
	acc_pdb.heavy_atom_temperature < 30 AND
	abs(don.resNum - acc.resNum) > 5;"


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
	ev.database_tag = p.database_tag AND ev.cosAHD_short = p.name;"
polynomials <- query_sample_sources(sample_sources, sele)



xmin <- min(f$AHdist, polynomials$xmin)
xmax <- max(f$AHdist, polynomials$xmax)


#coAHD goes from 0 to 1, where 1 is linear since there is significant
#density at 1, to accurately model a discontinuity, reflect the data
#across the right boundary, in computing the density esitmation
dens <- estimate_density_1d_reflect_boundary(
 data=f,
 ids = c("sample_source", "acc_chem_type", "don_chem_type"),
 variable = "cosAHD",
 reflect_right=TRUE,
 right_boundary=1, adjust=.7)


dens$y <- -log(dens$y)*.2

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
	z <- data.frame(x=x, y=predict(poly$poly, x),
        	sample_source=factor("Rosetta Model"))
	z$y <- ifelse(z$x >= poly$xmin & z$x <= poly$xmax, z$y, NA)
	z
})

dens.model$counts <- NA

dens <- rbind(dens, dens.model)


dens$don_chem_type_name <- don_chem_type_name_linear(dens$don_chem_type)
dens$acc_chem_type_name <- acc_chem_type_name_linear(dens$acc_chem_type)

dens <- dens[!is.na(dens$acc_chem_type_name) & !is.na(dens$don_chem_type_name),]

plot_id <- "hbond_cosAHD_chem_type_with_parameters"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(180-acos(x)*180/pi, y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	facet_grid(don_chem_type_name ~ acc_chem_type_name) +
	ggtitle("HBond AHD Angle by Chemical Type, SeqSep > 5, B-Factor < 30\nnormalized for equal weight per unit distance in density estimation") +
	scale_y_continuous("Energy (arbitrary units)", limits=c(-.6,1.1), breaks=c(-.5,0,.5)) +
	scale_x_continuous(
		"Acceptor -- Hydrogen -- Donor (degrees)", trans="reverse",
		limits=c(180, 120), breaks=c(180, 160, 140))

if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
