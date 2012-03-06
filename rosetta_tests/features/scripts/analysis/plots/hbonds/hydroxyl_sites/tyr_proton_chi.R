# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

                                            
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self){

#sele <-"
#SELECT
#	atoms.bbase_x AS CE2x, atoms.bbase_y AS CE2y, atoms.bbase_z AS CE2z,
# 	atoms.base_x  AS CZx,  atoms.base_y  AS CZy,  atoms.base_z  AS CZz,
#	atoms.atm_x   AS OHx,  atoms.atm_y   AS OHy,  atoms.atm_z   AS OHz,
#	atoms.base2_x AS HHx,  atoms.atm_y   AS HHy,  atoms.atm_z   AS HHz
#FROM
#	structures AS struct,
#	hbond_sites AS site,
#	hbond_site_atoms AS atoms,
#	hbond_sites_pdb AS pdb
#WHERE
#	site.struct_id = struct.struct_id AND
#	atoms.struct_id = site.struct_id AND atoms.site_id = site.site_id AND
#	site.HBChemType = 'hbacc_AHX' AND
#	pdb.struct_id = site.struct_id AND pdb.site_id = site.site_id AND
#	pdb.heavy_atom_temperature < 30;"
#f <- query_sample_sources(sample_sources, sele)
#
#f$proton_chi <- with(f, 
#	vector_dihedral(cbind(CE2x, CE2y, CE2z), cbind(CZx, CZy, CZz),
#		cbind(OHx, OHy, OHz), cbind(HHx, HHy, HHz)))
#
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
  geom_indicator(aes(indicator=counts, colour=sample_source)) +
  opts(title = "Proton Chi Torsional angle on Tyrosine Residues") +
  scale_x_continuous('Protein Chi (degrees)') +
  scale_y_continuous('Feature Density') +
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
