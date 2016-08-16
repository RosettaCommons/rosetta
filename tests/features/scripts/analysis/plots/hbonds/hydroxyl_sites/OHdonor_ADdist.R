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
id = "OHdonor_ADdist",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

source("scripts/methods/polynomial_methods.R")


extract_features <- function(sample_sources, don_chem_type, acc_chem_type, xmin, xmax){
	sele <-paste("
SELECT
	don_atoms.base_x AS dx, don_atoms.base_y AS dy, don_atoms.base_z AS dz,
	acc_atoms.atm_x  AS ax, acc_atoms.atm_y  AS ay, acc_atoms.atm_z  AS az
FROM
	hbonds AS hb,
	hbond_sites AS don, hbond_sites AS acc,
	hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb,
	hbond_site_atoms AS don_atoms, hbond_site_atoms AS acc_atoms
WHERE
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30 AND
	don_atoms.struct_id = hb.struct_id AND don_atoms.site_id = hb.don_id AND
	acc_atoms.struct_id = hb.struct_id AND acc_atoms.site_id = hb.acc_id AND
	don.HBChemType = '", don_chem_type, "' AND
	acc.HBChemType = '", acc_chem_type, "';", sep="")
	f <- query_sample_sources(sample_sources, sele)

        # A-D distance is not stored directly in the features database,
        # however it can be computed from the coordinates of the hydrogen
        # bonding atoms.
        transform(f,
        	ADdist = vector_distance(cbind(dx, dy, dz), cbind(ax, ay, az)))

}

estimate_density <- function(f){
	estimate_density_1d(
		data = f,
		ids = c("sample_source"),
		variable = "ADdist",
		weight_fun = radial_3d_normalization)
}


generate_plot <- function(
	dens,
	don_chem_type, acc_chem_type,
	sample_sources, output_dir, output_formats) {

	plot_id <- paste("OHdonor_ADdist", don_chem_type, acc_chem_type, sep="_")

	ref_dens <- transform(subset(dens, reference), sample_source=NULL)
	new_dens <- subset(dens, !reference)

	# hack!
	facet_labels <- data.frame(
		sample_source = factor(c("Relaxed Native Score12", "Relaxed Native NewHB", "Relaxed Native NewHB LJcorr")),
		label = c("A", "B", "C"))

	p <- ggplot() +
		geom_line(data=ref_dens, aes(x=x, y=y), colour="black", size=2) +
		geom_line(data=new_dens, aes(x=x, y=y, colour=model), size=2) +
		geom_indicator(data=facet_labels, aes(indicator=label),
			group=1, color="black", xpos=.03, ypos=.97, size=20) +
		facet_wrap(~sample_source, ncol=2) +
		scale_y_continuous("Feature Density") +
		scale_x_continuous(
			expression(paste('Acceptor -- Donor Distance (', ring(A), ')'))) +
		scale_colour_discrete("HBond Potential") +
		theme_bw() +
		theme(legend.position = c(.36, .85))

	save_plots(self, plot_id=plot_id, sample_sources, output_dir, output_formats)
}


##########################

xmin <- 1.5
xmax <- 2.2
don_chem_type <- "hbdon_HXL"
acc_chem_type <- "hbacc_PBA"

print(don_chem_type)
print(acc_chem_type)

f <- extract_features(
	sample_sources,
	don_chem_type,
	acc_chem_type,
	xmin,
	xmax)

#counts, e.g. to add to caption
print(ddply(f, .(sample_source), nrow))

dens <- estimate_density(f)

if(!("reference" %in% names(sample_sources))){
	stop("Please add a boolean 'reference' field to each sample source block in the analsis configuration.")

}

if(!("model" %in% names(sample_sources))){
	stop("Please add a string 'model' field to each sample source block in the analsis configuration.")
}



dens <- merge(dens, sample_sources[, c("sample_source", "reference", "model")])


generate_plot(
	dens,
	don_chem_type, acc_chem_type,
	sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
