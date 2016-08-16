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
id = "HAH_HXL_as_acceptor.R",
author = "Andrew Leaver-Fay",
brief_description = "Measure the angle of the tyrosine, serine, and threonine HAH angle when acting as an acceptor",

long_description = "
Measure the H-A-H angle for hydroxyl acceptors.  Since the flags that are active during extraction can decide whether
the cosBAH angle reported in the features database is either the CB-OG-H angle OR the HG-OG-H angle, its worth simply
measuring the angle directly from the cooridinates. This script creates HAH angle plots for SER/THR and for TYR
acceptors -- both aggregate plots, and plots broken down by donor-type.


                  D
                 /
                /
               H
   HG        '
     \     '
       \ '
        OG
        |
        |
     ---CB
       /|
      / |

",


feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	don.HBChemType AS don_chem_type,
	acc.HBChemType AS acc_chem_type,
	acc.resNum AS acc_residue_number,
	don.resNum AS don_residue_number,
	acc_atoms.base2_x AS Bx, acc_atoms.base2_y AS By, acc_atoms.base2_z AS Bz,
	acc_atoms.atm_x  AS Ax, acc_atoms.atm_y AS Ay,  acc_atoms.atm_z AS Az,
	don_atoms.atm_x  AS Hx, don_atoms.atm_y AS Hy,  don_atoms.atm_z AS Hz
FROM
	hbonds AS hb,
	hbond_sites AS acc,
  hbond_sites AS don,
	hbond_site_atoms AS acc_atoms,
	hbond_site_atoms AS don_atoms,
  hbond_sites_pdb AS acc_pdb,
  hbond_sites_pdb AS don_pdb
WHERE
	acc_atoms.struct_id = hb.struct_id AND acc_atoms.site_id = hb.acc_id AND
	don_atoms.struct_id = hb.struct_id AND don_atoms.site_id = hb.don_id AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	(acc.HBChemType = 'hbacc_HXL' OR acc.HBChemType = 'hbacc_AHX') AND
	abs(acc.resNum - don.resNum) > 5 AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.heavy_atom_temperature < 30;"

f <- query_sample_sources(sample_sources, sele)

f <- transform(f,
  cosBAH = vector_dotprod(
    vector_normalize(cbind(Ax-Bx, Ay-By, Az-Bz)),
		vector_normalize(cbind(Hx-Ax, Hy-Ay, Hz-Az))))

dens <- estimate_density_1d(
  data = f,
  ids = c("sample_source", "acc_chem_type"),
  variable = "cosBAH")

plot_id = "cosHAH_hydroxyl_acceptors"
p <- ggplot(data=dens) + theme_bw() +
  geom_line(aes(x=180-acos(x)*180/pi, y=y, colour=sample_source)) +
  geom_indicator(aes(colour=sample_source, indicator=counts, group=sample_source)) +
  facet_wrap( ~ acc_chem_type) +
  ggtitle("Hydrogen Bonds HAH Angle for Hydroxyl Acceptors, Measured from Coordinates\nB-Factors < 30; Sequence Separation > 5") +
  scale_x_continuous(paste('Hydroxyl Hydrogen -- Acceptor -- Donor Hydrogen (degrees)')) +
  scale_y_continuous("FeatureDensity")
if(nrow(sample_sources) <= 3){
  p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

dens2 <- estimate_density_1d(
  data = f,
  ids = c("sample_source", "acc_chem_type", "don_chem_type"),
  variable = "cosBAH")

plot_id = "cosHAH_hydroxyl_acceptors_chem_type"
p <- ggplot(data=dens2) + theme_bw() +
  geom_line(aes(x=180-acos(x)*180/pi, y=y, colour=sample_source)) +
  geom_indicator(aes(colour=sample_source, indicator=counts, group=sample_source)) +
  facet_grid( don_chem_type ~ acc_chem_type) +
  ggtitle("Hydrogen Bonds HAH Angle for Hydroxyl Acceptors, Measured from Coordinates\nB-Factors < 30; Sequence Separation > 5") +
  scale_x_continuous(paste('Hydroxyl Hydrogen -- Acceptor -- Donor Hydrogen (degrees)')) +
  scale_y_continuous("FeatureDensity")
if(nrow(sample_sources) <= 3){
  p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
