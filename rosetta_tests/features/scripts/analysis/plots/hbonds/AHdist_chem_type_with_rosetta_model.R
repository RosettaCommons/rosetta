# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

check_setup()

source("scripts/analysis/plots/hbonds/hbond_geo_dim_scales.R")
source("scripts/parameter_analysis/hbonds/methods/methods.R")

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
	ev.separation,
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
polynomials <- query_sample_sources(sample_sources[2,], sele)


xmin <- min(f$AHdist)
xmax <- max(f$AHdist)

dens <- estimate_density_1d(
	f, c("sample_source", "acc_chem_type", "don_chem_type"),
	"AHdist", weight_fun = radial_3d_normalization, adjust=.6)


dens$y <- -log(dens$y)*.35

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
  separation <- "seq_sep_other"
  dimension <- "hbgd_AHdist"

  poly <- get_1d_polynomial_from_types(
    polynomials,
    sample_source,
    acc_chem_type, don_chem_type, separation, dimension)
  x <- seq(xmin, xmax, length.out=100)
  data.frame(x=x, y=aaply(x, 1, as.function(poly)), sample_source=factor("Rosetta Model"))
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

dens$sample_source <- factor(dens$sample_source,
	levels = levels(dens$sample_source),
	labels = c("top8000 (native)", "Rosetta (simulated)", "AHdist dimension of Rosetta HBond module"))

plot_id <- "hbond_AHdist_chem_type_with_parameters"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x, y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source)) +
	facet_grid(don_chem_type_name ~ acc_chem_type_name) +
	opts(title = "HBond A-H Distance by Chemical Type, B-Factor < 30\nnormalized for equal weight per unit distance in density estimation") +
	scale_y_continuous("Energy (arbitrary units)", limits=c(-.7,.2), breaks=c(-.5,-.25,0)) +
	scale_x_continuous(expression(paste('Acceptor -- Proton Distance (', ring(A), ')')), limits=c(1.4,2.7), breaks=c(1.6, 1.9, 2.2, 2.6))

if(nrow(sample_sources) <= 3){
	p <- p + opts(legend.position="bottom", legend.direction="horizontal")
}

save_plots(plot_id, sample_sources, output_dir, output_formats)
