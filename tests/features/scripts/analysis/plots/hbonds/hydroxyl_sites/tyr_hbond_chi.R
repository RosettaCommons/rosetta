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
id = "tyr_hbond_chi",
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
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	don_site.HBChemType AS don_chem_type,
	CASE don_site.HBChemType
		WHEN 'hbdon_IMD' THEN 'imidazole: h'
		WHEN 'hbdon_IME' THEN 'imidazole: h'
		WHEN 'hbdon_IND' THEN 'indol: w'
		WHEN 'hbdon_AHX' THEN 'hydroxyl: s,t,y'
		WHEN 'hbdon_HXL' THEN 'hydroxyl: s,t,y'
		WHEN 'hbdon_GDE' THEN 'other'
		WHEN 'hbdon_GDH' THEN 'other'
		WHEN 'hbdon_PBA' THEN 'other'
		WHEN 'hbdon_CXA' THEN 'other'
		WHEN 'hbdon_AMO' THEN 'other'
			END AS don_hybrid,
	ABS(don_site.resNum - acc_site.resNum) AS seq_sep,
	acc_atoms.bbase_x AS CE2x, acc_atoms.bbase_y AS CE2y, acc_atoms.bbase_z AS CE2z,
 	acc_atoms.base_x  AS CZx,  acc_atoms.base_y  AS CZy,  acc_atoms.base_z  AS CZz,
	acc_atoms.atm_x   AS OHx,  acc_atoms.atm_y   AS OHy,  acc_atoms.atm_z   AS OHz,
	don_atoms.atm_x   AS Hx,   don_atoms.atm_y   AS Hy,   don_atoms.atm_z   AS Hz
FROM
	hbonds AS hb,
	hbond_sites AS don_site, hbond_sites AS acc_site,
	hbond_site_atoms AS don_atoms, hbond_site_atoms AS acc_atoms,
	hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb
WHERE
	hb.struct_id = acc_site.struct_id AND hb.acc_id = acc_site.site_id AND
	hb.struct_id = don_site.struct_id AND hb.don_id = don_site.site_id AND
	abs( acc_site.resNum - don_site.resNum ) > 5 AND
	don_atoms.struct_id = hb.struct_id AND don_atoms.site_id = hb.don_id AND
	acc_atoms.struct_id = hb.struct_id AND acc_atoms.site_id = hb.acc_id AND
	acc_site.HBChemType = 'hbacc_AHX' AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30;"
f <- query_sample_sources(sample_sources, sele)

f$hbond_chi <- with(f,
	vector_dihedral(cbind(CE2x, CE2y, CE2z), cbind(CZx, CZy, CZz),
	cbind(OHx, OHy, OHz), cbind(Hx, Hy, Hz)))


f$don_chem_type <- factor(f$don_chem_type,
	levels=c("hbdon_IMD", "hbdon_GDE", "hbdon_HXL", "hbdon_PBA",
		"hbdon_IME", "hbdon_GDH", "hbdon_AHX", "hbdon_IND",
		"hbdon_CXA", "hbdon_AMO"),
	labels = c("dIMD: h", "dGDE: r", "dHXL: s,t", "dPBA: bb",
		"dIME: h", "dGDH: r", "dAHX: y", "dIND: w",
		"dCXA: n,q", "dAMO: k"))

dens <- estimate_density_1d_wrap(f,
  c("sample_source", "don_chem_type"), "hbond_chi", xlim=c(-pi, pi), adjust=.3)

plot_id = "tyr_hbond_chi"
ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x*180/pi, y=y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	facet_wrap( ~ don_chem_type ) +
	ggtitle("HBond Torsional Angle for Tyrosine Acceptors") +
	scale_x_continuous('HBond Torsional Angle (degrees)') +
	scale_y_continuous('Feature Density') +
	theme(legend.position=c(.58,.35)) +
	theme(legend.justification=c("left", "top"))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

dens <- estimate_density_1d_wrap(f,
  c("sample_source", "don_hybrid"), "hbond_chi", xlim=c(-pi, pi), adjust=.2)

plot_id = "tyr_hbond_chi_by_hybrid"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x*180/pi, y=y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	facet_wrap( ~ don_hybrid ) +
	ggtitle("HBond Torsional Angle for Tyrosine Acceptors") +
	scale_x_continuous('HBond Torsional Angle (degrees)') +
	scale_y_continuous('Feature Density')
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


dens <- estimate_density_1d_wrap(f[f$seq_sep > 5,],
  c("sample_source", "don_hybrid"), "hbond_chi", xlim=c(-pi, pi), adjust=.2)

plot_id = "tyr_hbond_chi_by_hybrid_seq_sep_5"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x*180/pi, y=y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	facet_wrap( ~ don_hybrid ) +
	ggtitle("HBond Torsional Angle for Tyrosine Acceptors SeqSep > 5") +
	scale_x_continuous('HBond Torsional Angle (degrees)') +
	scale_y_continuous('Feature Density')
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)



})) # end FeaturesAnalysis
