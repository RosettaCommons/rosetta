# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "ADdist_chem_type",
author = "Matthew O'Meara",

brief_description = "Measure the distance between acceptor heavy atoms and the donor heavy atoms for hydrogen bonds.",


long_description = "
Hydrogen Bonds are chemical contacts between a negatively charged
acceptor group and a positively charged donor group. This script measures
the distance between pairs of donor and acceptor atoms. The A-D distance
can be computed from the coordinates of the atoms, which are stored in the
hbond_site_atoms table. The ADchem_type plot groups A-D distance observations
by the chemical type of the donor group and the chemical type of the acceptor
group displaying them as a lattice of plots.

The null model assumes that in the coordinate frame of the acceptor,
the position of the donor atom is uniformaly distributed in cartesian space.
To compute the density in the radial dimension requires normalizing by
the 1/r^2 so the plotted density of the null distribution is uniform. See
the hbond_null analysis script for more details.

Plotting the log(FeatureDensity + 1) allows high density and low density
aspects of the distribution to be included in roughly the same scale. In other
words, the subplot with the highest density sets the scale for all the plots.
If there is a highly concentrated region this can set the scale large to see
low density aspects of the distribution. The '+ 1' acts a pseudo-count,
making the logarithm of never go below zero. Other choices for the value of the
pseudo count would control minimum value.

The indicator in each subplot counts the number of hydrogen bonds in
each group for each sample source. For example, the red number in the
('aAHX: y', 'dPBA: bb') cell is the number of hydrogen bonds in the
sample source corresponding to the red distribution with an aromatic
hydroxyl acceptor (e.g., in a tyrosine) a protein backbone donor
(i.e., the N-H group of a protein backbone).",

feature_reporter_dependencies = c("HBondFeatures"),

run=function(self, sample_sources, output_dir, output_formats){

source("scripts/analysis/plots/hbonds/hbond_geo_dim_scales.R")

sele <- "
SELECT
	don_atoms.base_x AS dx, don_atoms.base_y AS dy, don_atoms.base_z AS dz,
	acc_atoms.atm_x  AS ax, acc_atoms.atm_y  AS ay, acc_atoms.atm_z  AS az,
	don.HBChemType AS don_chem_type, acc.HBChemType AS acc_chem_type
FROM
	hbonds AS hb,
	hbond_sites AS don, hbond_sites AS acc,
	hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb,
	hbond_site_atoms AS don_atoms, hbond_site_atoms AS acc_atoms
WHERE
	don.struct_id = hb.struct_id AND don.site_id =hb.don_id AND
	acc.struct_id = hb.struct_id AND acc.site_id =hb.acc_id AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30 AND
	don_atoms.struct_id = hb.struct_id AND don_atoms.site_id = hb.don_id AND
	acc_atoms.struct_id = hb.struct_id AND acc_atoms.site_id = hb.acc_id;"


# Execute the SQL query on each sample source.
f <- query_sample_sources(sample_sources, sele)

# A-D distance is not stored directly in the features database,
# however it can be computed from the coordinates of the hydrogen
# bonding atoms.
f <- transform(f,
	ADdist = vector_distance(cbind(dx, dy, dz), cbind(ax, ay, az)))

# This shouldn't happend, but if it does get rid of them
f <- f[f$ADdist < 5,]

f$don_chem_type_name <- don_chem_type_name_linear(f$don_chem_type)
f$acc_chem_type_name <- acc_chem_type_name_linear(f$acc_chem_type)

f <- na.omit(f, method="r")

# Compute density estimation for over the A-D distance grouping by the
# donor type, acceptor type and sample source. Apply the radial 3d
# normalization. This corrects for the fact that there is more volume
# in a spherical shell at a farther distance then a closer distance.
dens <- estimate_density_1d(f,
	c("sample_source", "don_chem_type_name", "acc_chem_type_name"),
	"ADdist", radial_3d_normalization)

# Generate a lattice of density plots for each donor and acceptor type
plot_id <- "hbond_ADdist_chem_type"
p <- ggplot(dens) + theme_bw() +
	geom_line(aes(x, y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	facet_grid(don_chem_type_name ~ acc_chem_type_name) +
	ggtitle("H-Bond A-D Distance by Chemical Type, B-Factor < 30\n(normalized for equal volume per unit distance)") +
	scale_x_ADdist +
	scale_y_continuous("FeatureDensity", limits=c(0,8.5), breaks=c(1,3,5,7))
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


#summary_stats <- compute_summary_statistics_1d(f,
#	c("sample_source", "don_chem_type_name", "acc_chem_type_name"),
#	"ADdist", radial_3d_normalization, c("primary_mode", "count")
#
#comp_stats <- compute_comparison_statistics_1d(self, f, sample_sources,
#	c("don_chem_type_name", "acc_chem_type_name"),
#	"ADdist", c("primary_mode_diff", "two_sided_ttest", "kolmogorov_smirnov_test")
#
#table_id <- "hbond_ADdist_chem_type_comp_stats"
#table_title <- "H-Bond A-D Distance by Chemical Type, B-Factor < 30 Comparison Statistics"
#save_tables(self, f, table_id, output_dir, output_formats,
#	caption=table_title, caption.placement="top")
#
#save_stats(self, "hbond_ADdist_chem_type", output_dir, output_formats, comp_stats)

})) # end FeaturesAnalysis
