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
id = "RealBAH_HXL_as_acceptor.R",
author = "Andrew Leaver-Fay",
brief_description = "Measure the angle of the tyrosine, serine, and threonine BAH angle -- where the base is the heavy-atom base of the hydroxyl oxygen -- when acting as an acceptor",

long_description = "
Currently the BAH angle extracted from the database is measured wrt to the hydroxyl hydrogen instead of the hydroxyl oxygen's heavyatom base.
As shown for serine below, the angle HG-OG-H will be reported as the BAH angle.  This is a tough angle to work with, though, because
the coordinate of the HG atom may only be inferred and not directly observed.  Instead, it would be good to look at the CB-OG-H angle.


                  A
                 '
                '
               H
   HG        /
     \     /
       \ /
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
	acc_atoms.base_x AS Bx, acc_atoms.base_y AS By, acc_atoms.base_z AS Bz,
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

plot_id = "cosBAH_hydroxyl_acceptors"
p <- ggplot(data=dens) + theme_bw() +
  geom_line(aes(x=180-acos(x)*180/pi, y=y, colour=sample_source)) +
  geom_indicator(aes(colour=sample_source, indicator=counts, group=sample_source)) +
  facet_wrap( ~ acc_chem_type) +
  ggtitle("Hydrogen Bonds BAH Angle for Hydroxyl Acceptors, Measured from the Heavy-Atom Base\nB-Factors < 30; Sequence Separation > 5") +
  scale_x_continuous(paste('Heavyatom Base -- Acceptor -- Hydrogen (degrees)')) +
  scale_y_continuous("FeatureDensity")
if(nrow(sample_sources) <= 3){
  p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

dens2 <- estimate_density_1d(
  data = f,
  ids = c("sample_source", "acc_chem_type", "don_chem_type"),
  variable = "cosBAH")

plot_id = "cosBAH_hydroxyl_acceptors_chem_type"
p <- ggplot(data=dens2) + theme_bw() +
  geom_line(aes(x=180-acos(x)*180/pi, y=y, colour=sample_source)) +
  geom_indicator(aes(colour=sample_source, indicator=counts, group=sample_source)) +
  facet_grid( don_chem_type ~ acc_chem_type) +
  ggtitle("Hydrogen Bonds BAH Angle for Hydroxyl Acceptors, Measured from the Heavy-Atom Base\nB-Factors < 30; Sequence Separation > 5") +
  scale_x_continuous(paste('Heavyatom Base -- Acceptor -- Hydrogen (degrees)')) +
  scale_y_continuous("FeatureDensity")
if(nrow(sample_sources) <= 3){
  p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
