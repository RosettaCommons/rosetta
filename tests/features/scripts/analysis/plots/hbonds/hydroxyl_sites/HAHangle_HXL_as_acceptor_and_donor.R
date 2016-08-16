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
id = "HAH_HXL_as_acceptor_and_donor.R",
author = "Andrew Leaver-Fay",
brief_description = "Measure the angle of the tyrosine, serine, and threonine HAH angle when acting as an acceptor",

long_description = "
Look for the distribution of HAH angles -- measured from the hydroxyl hydrogen -- when the hydroxyl group
is acting as both a donor and an acceptor.  In such cases, the expectation would be that the hydrogen's location
can be more accurately inferred than in cases where the hydroxyl is only acting as an acceptor and there's little
other information in the vicinity of the hydroxyl group that indicates where the hydroxyl hydrogen might be.


                  D
                 /
                /
               H
   HG         '
     \      '
       \  '
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
CREATE TEMPORARY TABLE hydroxyl_acceptor_sites_that_also_donate (
  struct_id INTEGER,
  site_id INTEGER,
  resNum INTEGER,
  HBChemType TEXT,
  PRIMARY KEY (struct_id,site_id));
INSERT INTO hydroxyl_acceptor_sites_that_also_donate
SELECT DISTINCT
  hxl_acc.struct_id AS struct_id,
  hxl_acc.site_id AS site_id,
  hxl_acc.resNum as resNum,
  hxl_acc.HBChemType as HBChemType
FROM
  hbonds AS hb1,
  hbonds AS hb2,
  hbond_sites AS don1,
  hbond_sites AS hxl_acc,
  hbond_sites AS hxl_don,
  hbond_sites AS acc2
WHERE
  hb1.struct_id = hb2.struct_id AND
  hb1.struct_id = don1.struct_id AND hb1.don_id = don1.site_id AND
  hb1.struct_id = hxl_acc.struct_id AND hb1.acc_id = hxl_acc.site_id AND
  (hxl_acc.HBChemType = 'hbacc_HXL' OR hxl_acc.HBChemType = 'hbacc_AHX') AND
  hb2.struct_id = hxl_don.struct_id AND hb2.don_id = hxl_don.site_id AND
  hb2.struct_id = acc2.struct_id AND hb2.acc_id = acc2.site_id AND
  (hxl_don.HBChemType = 'hbdon_HXL' OR hxl_don.HBChemType = 'hbdon_AHX' ) AND
  hxl_don.resNum = hxl_acc.resNum;
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
	hydroxyl_acceptor_sites_that_also_donate AS acc,
  hbond_sites AS don,
  hbond_sites_pdb AS acc_pdb,
  hbond_sites_pdb AS don_pdb,
  hbond_site_atoms AS acc_atoms,
  hbond_site_atoms AS don_atoms
WHERE
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	acc_atoms.struct_id = hb.struct_id AND acc_atoms.site_id = hb.acc_id AND
	don_atoms.struct_id = hb.struct_id AND don_atoms.site_id = hb.don_id AND
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

plot_id = "cosHAH_hydroxyl_acceptors_that_also_donate"
p <- ggplot(data=dens) + theme_bw() +
  geom_line(aes(x=180-acos(x)*180/pi, y=y, colour=sample_source)) +
  geom_indicator(aes(colour=sample_source, indicator=counts, group=sample_source)) +
  facet_wrap( ~ acc_chem_type) +
  ggtitle("Hydrogen Bonds HAH Angle for Hydroxyl Acceptors where the Hydroxyl also Donates\nB-Factors < 30; Sequence Separation > 5") +
  scale_x_continuous(paste('Acceptor Hydrogen -- Acceptor -- Donor Hydrogen (degrees)')) +
  scale_y_continuous("FeatureDensity")
if(nrow(sample_sources) <= 3){
  p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

dens2 <- estimate_density_1d(
  data = f,
  ids = c("sample_source", "acc_chem_type", "don_chem_type"),
  variable = "cosBAH")

plot_id = "cosHAH_hydroxyl_acceptors_that_also_donate_chem_type"
p <- ggplot(data=dens2) + theme_bw() +
  geom_line(aes(x=acos(x)*180/pi, y=y, colour=sample_source)) +
  geom_indicator(aes(colour=sample_source, indicator=counts, group=sample_source)) +
  facet_grid( don_chem_type ~ acc_chem_type) +
  ggtitle("Hydrogen Bonds BAH Angle for Hydroxyl Acceptors where the Hydroxyl also Donates\nB-Factors < 30; Sequence Separation > 5") +
  scale_x_continuous(paste('Acceptor Hydrogen -- Acceptor -- Donor Hydrogen (degrees)')) +
  scale_y_continuous("FeatureDensity")
if(nrow(sample_sources) <= 3){
  p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
