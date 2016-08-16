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
id = "AHdist_by_rank",
author = "Matthew O'Meara",
brief_description = "",
long_description = "
The rank of a hydrogen bond at donor site or acceptor site is rank of the
relative Rosetta HBond energy of the hydrogen bond at the site. A rank of
0 indicates the hydrogen bond is the only hydrogen bond at the site.

The indicator in each subplot counts the number of hbonds in each group for
each sample sample source.",

feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

source("scripts/analysis/plots/hbonds/hbond_geo_dim_scales.R")

sele <-"
SELECT
	geom.AHdist,
	acc.HBChemType AS acc_chem_type, don.HBChemType AS don_chem_type,
	hb.accRank AS accRank, hb.donRank AS donRank
FROM
	hbond_geom_coords AS geom,
	hbonds AS hb,
	hbond_sites AS don, hbond_sites AS acc,
	hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb
WHERE
	geom.struct_id = hb.struct_id AND geom.hbond_id = hb.hbond_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30;"

f <- query_sample_sources(sample_sources, sele)
f$accRank <- factor(f$accRank)
f$donRank <- factor(f$donRank)

f$don_chem_type_name <- don_chem_type_name_linear(f$don_chem_type)
f$acc_chem_type_name <- acc_chem_type_name_linear(f$acc_chem_type)

plot_parts <- list(
  geom_line(aes(x, y)),
  geom_indicator(aes(indicator=counts)),
  labs(x=expression(paste('Acceptor -- Proton Distance (', ring(A), ')')),
       y="FeatureDensity"),
  theme_bw(),
	theme(legend.position="bottom", legend.direction="horizontal"))

if(nrow(sample_sources) <= 3){
	plot_parts_all <- c(plot_parts,
		list(
			facet_wrap( ~ sample_source, nco=1),
			list(theme(legend.position="bottom", legend.direction="horizontal"))))
} else {
	plot_parts_all <- c(plot_parts,
		list(facet_wrap( ~ sample_source)))
}

plot_id <- "hbond_AHdist_by_donRank"
d <- estimate_density_1d(f, c("sample_source", "donRank"), "AHdist", radial_3d_normalization)
ggplot(d, aes(colour=donRank, group=donRank)) + plot_parts_all +
  ggtitle("Hydrogen Bonds A-H Distance by Donor Rank, B-Factor < 30\nnormalized for equal weight per unit distance")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id <- "hbond_AHdist_by_accRank"
d <- estimate_density_1d(f, c("sample_source", "accRank"), "AHdist", radial_3d_normalization)
ggplot(d, aes(colour=accRank, group=accRank)) + plot_parts_all +
  ggtitle("HBonds A-H Distance by Acceptor Rank, B-Factor < 30\nnormalized for equal weight per unit distance")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

f1 <- f[f$donRank==1,]
plot_id <- "hbond_AHdist_donRank_is_1"
d <- estimate_density_1d(f1, c("sample_source"), "AHdist", radial_3d_normalization)
ggplot(d, aes(colour=sample_source, group=sample_source)) + plot_parts +
  ggtitle("HBonds A-H Distance, Primary HBond at Bifuracted Donor, B-Factor < 30\nnormalized for equal weight per unit distance")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

f1 <- f[f$donRank==1,]
plot_id <- "hbond_AHdist_chem_type_donRank_is_1"
d <- estimate_density_1d(f1, c("sample_source", "don_chem_type_name", "acc_chem_type_name"), "AHdist", radial_3d_normalization)
ggplot(d, aes(colour=sample_source, group=sample_source)) + plot_parts +
  facet_grid( don_chem_type_name ~ acc_chem_type_name) +
  ggtitle("HBonds A-H Distance, Donor Rank=1, B-Factor < 30\nnormalized for equal weight per unit distance")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

f1 <- f[f$accRank==1,]
plot_id <- "hbond_AHdist_accRank_is_1"
d <- estimate_density_1d(f1, c("sample_source"), "AHdist", radial_3d_normalization)
ggplot(d, aes(colour=sample_source, group=sample_source)) + plot_parts +
  ggtitle("HBonds A-H Distance, Primary HBond at Bifuracted Acceptor, B-Factor < 30\nnormalized for equal weight per unit distance")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

f1 <- f[f$accRank==1,]
plot_id <- "hbond_AHdist_chem_type_accRank_is_1"
d <- estimate_density_1d(f1, c("sample_source", "acc_chem_type_name", "acc_chem_type_name"), "AHdist", radial_3d_normalization)
ggplot(d, aes(colour=sample_source, group=sample_source)) + plot_parts +
  facet_grid( don_chem_type_name ~ acc_chem_type_name) +
  ggtitle("HBonds A-H Distance, Acceptor Rank=1, B-Factor < 30\nnormalized for equal weight per unit distance")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
