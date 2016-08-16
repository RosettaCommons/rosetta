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
id = "tyr_proton_chi",
author = "Matthew O'Meara",
brief_description = "Measure the angle of the tyrosine dihedral angle of the hydroxyl hydrogen and the best hydrogen bonded hydrogen relative to the plane of the aromatic ring",

long_description = "
The dihedral angle in tyrosine residues defined by the atoms (CE2, C, OH, HH) may be rotameric because the C-OH bond is has partial sp2 character because of the resonance involving the aromatic ring. It is challening to definitively place the HH atom. Therefore use the dihedral angle (CE2, C, OH, H) where the final H is the hydrogen atom of the best hydrogen bond donor making a hydrogen bond with the OH group.

Currently the atoms extracted for hydrogen bond sites are the atm, base, bbase, base2. For the OH acceptors on tyrosines have acceptor chemical type hbacc_AHX (for aromatic hydroxyl) and have the following atoms:

  atm:  OH
  base: CZ
  bbase: CE2
  base2: HH

  H         HH
   `       /
     `    /
       `OH
        |
        |
        CZ
       /  \
      /    \
     |      CE2
     |  ()  |
     |      |
      \    /
       \ _/

  proton_chi = DIHEDRAL(CE2, CZ, OH, HH)
  hbond_chi = DIHEDRAL(CE2, CZ, OH, H)
",


feature_reporter_dependencies = c("ResidueFeatures", "HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <- "
SELECT
	res_conf.chi3 AS proton_chi
FROM
	residues AS res,
	protein_residue_conformation AS res_conf
WHERE
	res.name3 = 'TYR' AND
	res_conf.struct_id = res.struct_id AND
	res_conf.seqpos = res.resNum;"
f <- query_sample_sources(sample_sources, sele)

f <- na.omit(f, method="r")
dens <- estimate_density_1d_wrap(
  f, c("sample_source"), "proton_chi", xlim=c(-180, 180), adjust=.1)


plot_id = "tyr_proton_chi"
ggplot(data=dens) + theme_bw() +
  geom_line(aes(x=x, y=y, colour=sample_source)) +
  geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
  ggtitle("Proton Chi Torsional angle on Tyrosine Residues") +
  scale_x_continuous('Protein Chi (degrees)') +
  scale_y_continuous('Feature Density') +
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
