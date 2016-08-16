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
id = "HXL_as_donor_hbond_chi",
author = "Matthew O'Meara",
brief_description = "Measure the angle of the serine and threonine dihedral angle of the hydroxyl hydrogen and the best hydrogen bonded hydrogen relative to the plane of the aromatic ring",

long_description = "
The dihedral angle in tyrosine residues defined by the atoms (CE2, C, OH, HH) may be rotameric because the C-OH bond is has partial sp2 character because of the resonance involving the aromatic ring. It is challening to definitively place the HH atom. Therefore use the dihedral angle (CE2, C, OH, H) where the final H is the hydrogen atom of the best hydrogen bond donor making a hydrogen bond with the OH group.

Currently the atoms extracted for hydrogen bond sites are the atm, base, bbase, base2. For the OH acceptors on tyrosines have acceptor chemical type hbacc_AHX (for aromatic hydroxyl) and have the following atoms:


                 A
                '
               '
              '
            HG
           /
          /
        OG
        |
        |
     ---CB
       /|
      / |


  proton_chi = DIHEDRAL(CA, CB, OG, HG)
  hbond_chi = DIHEDRAL(CA, CB, OG, A)
",


feature_reporter_dependencies = c("StructureFeatures", "HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
CREATE TEMPORARY TABLE HXL_sites(
	struct_id INTEGER,
	don_id INTEGER,
	resNum INTEGER,
	CAx REAL, CAy REAL, CAz REAL,
	CBx REAL, CBy REAL, CBz REAL,
	OGx REAL, OGy REAL, OGz REAL,
	HGx REAL, HGy REAL, HGz REAL,
	PRIMARY KEY(struct_id, don_id));

INSERT INTO HXL_sites
SELECT
	don.struct_id,
	don.site_id AS don_id,
	don.resNum,
	acc_atoms.bbase_x AS CAx, acc_atoms.bbase_y AS CAy, acc_atoms.bbase_z AS CAz,
	acc_atoms.base_x  AS CBx,  acc_atoms.base_y AS CBy, acc_atoms.base_z  AS CBz,
 	don_atoms.base_x  AS OGx,  don_atoms.base_y AS OGy, don_atoms.base_z  AS OGz,
	don_atoms.atm_x   AS HGx,  don_atoms.atm_y  AS HGy, don_atoms.atm_z   AS HGz
FROM
	hbond_sites AS don, hbond_sites AS acc,
	hbond_site_atoms AS don_atoms, hbond_site_atoms AS acc_atoms,
	hbond_sites_pdb AS pdb
WHERE
	acc.struct_id = don.struct_id AND acc.resNum = don.resNum AND

	don_atoms.struct_id = don.struct_id AND
	don_atoms.site_id = don.site_id AND

	acc_atoms.struct_id = acc.struct_id AND
	acc_atoms.site_id = acc.site_id AND

	pdb.struct_id = don.struct_id AND pdb.site_id = don.site_id AND
	pdb.heavy_atom_temperature < 30 AND
	don.HBChemType = 'hbdon_HXL' AND acc.HBChemType = 'hbacc_HXL';

SELECT
	struct.tag,
	acc.HBChemType AS acc_chem_type,
	acc.resNum AS acc_residue_number,
	hxl.resNum AS hxl_residue_number,
	ABS(hxl.resNum - acc.resNum) AS seq_sep,
	CASE acc.HBChemType
			WHEN 'hbacc_IMD' THEN 'ring' WHEN 'hbacc_IME' THEN 'ring'
			WHEN 'hbacc_AHX' THEN 'sp3'  WHEN 'hbacc_HXL' THEN 'sp3'
			WHEN 'hbacc_CXA' THEN 'sp2'  WHEN 'hbacc_CXL' THEN 'sp2'
			WHEN 'hbacc_PBA' THEN 'sp2'  END AS acc_hybrid,
	hb.donRank AS don_rank,
	hxl.CAx, hxl.CAy, hxl.CAz,
	hxl.CBx, hxl.CBy, hxl.CBz,
	hxl.OGx, hxl.OGy, hxl.OGz,
	hxl.HGx, hxl.HGy, hxl.HGz,
	acc_atoms.atm_x AS Ax, acc_atoms.atm_y AS Ay, acc_atoms.atm_z AS Az
FROM
	structures AS struct,
	hbonds AS hb,
	hxl_sites AS hxl, hbond_sites AS acc,
	hbond_site_atoms AS acc_atoms,
	hbond_sites_pdb AS acc_pdb
WHERE
	hb.struct_id = struct.struct_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	hxl.struct_id = hb.struct_id AND hxl.don_id = hb.don_id AND
	abs(acc.resNum - hxl.resNum) > 5 AND

	acc_atoms.struct_id = acc.struct_id AND
	acc_atoms.site_id = acc.site_id AND

	acc_pdb.struct_id = acc.struct_id AND
	acc_pdb.site_id = acc.site_id AND
	acc_pdb.heavy_atom_temperature < 30;"

f <- query_sample_sources(sample_sources, sele)

f$hbond_chi <- with(f,
	vector_dihedral(cbind(CAx, CAy, CAz), cbind(CBx, CBy, CBz),
	cbind(OGx, OGy, OGz), cbind(Ax, Ay, Az)))

f$proton_chi <- with(f,
	vector_dihedral(cbind(CAx, CAy, CAz), cbind(CBx, CBy, CBz),
	cbind(OGx, OGy, OGz), cbind(HGx, HGy, HGz)))

f$chi_diff <- f$hbond_chi - f$proton_chi
f$chi_diff <- with(f, ifelse(chi_diff < pi, chi_diff, chi_diff - 2*pi))
f$chi_diff <- with(f, ifelse(chi_diff > -pi, chi_diff, chi_diff + 2*pi))

###since tyrosines are symmetric, don't distinguish points points close
###to 0 different from points close to pi
#f$hbond_chi <- abs(f$raw_hbond_chi)
#f$hbond_chi <- ifelse(f$hbond_chi < pi/2, f$hbond_chi, pi - f$hbond_chi)


adjust = .1
adjust_hybrid = .1
adjust_diff = .3
adjust_diff_hybrid = .3

f <- na.omit(f, method="r")

# Order the plots better and give more descriptive labels
f$acc_chem_type <- factor(f$acc_chem_type,
	levels = c("hbacc_IMD", "hbacc_IME", "hbacc_AHX", "hbacc_HXL",
		"hbacc_CXA", "hbacc_CXL", "hbacc_PBA"),
	labels = c("aIMD: h", "aIME: h", "aAHX: y", "aHXL: s,t",
		"aCXA: n,q", "aCXL: d,e", "aPBA: bb"))

dens <- estimate_density_1d_wrap(f,
  c("sample_source", "acc_chem_type"), "hbond_chi", xlim=c(-pi, pi),
	adjust=adjust)


plot_id = "hxl_as_donor_hbond_chi"
ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x*180/pi, y=y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	facet_wrap( ~ acc_chem_type ) +
	ggtitle("HBond Torsional Angle for Hydroxyl As Donors: (CA, CB, OG, A) SeqSep > 5, BFact < 30") +
	scale_x_continuous('HBond Torsional Angle (degrees)') +
	scale_y_continuous('Feature Density') +
	theme(legend.position=c(.58,.35)) +
	theme(legend.justification=c("left", "top"))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


dens <- estimate_density_1d_wrap(f,
  c("sample_source", "acc_hybrid"), "hbond_chi", xlim=c(-pi, pi),
	adjust=adjust_hybrid)

plot_id = "hxl_as_donor_hbond_chi_by_hybrid"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x*180/pi, y=y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	facet_wrap( ~ acc_hybrid ) +
	ggtitle("HBond Torsional Angle for Hydroxyl (s,t) As Donors: (CA, CB, OG, A) SeqSep > 5, BFact < 30") +
	scale_x_continuous('HBond Torsional Angle (degrees)') +
	scale_y_continuous('Feature Density') +
	theme(legend.position=c(.58,.35)) +
	theme(legend.justification=c("left", "top"))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


dens <- estimate_density_1d_wrap(f,
  c("sample_source", "acc_chem_type"), "proton_chi", xlim=c(-pi, pi),
	adjust=adjust)

plot_id = "hxl_as_donor_proton_chi"
ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x*180/pi, y=y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	facet_wrap( ~ acc_chem_type ) +
	ggtitle("Proton Torsional Angle for Hydroxyl (s,t): (CA, CB, OG, HG) SeqSep > 5, BFact < 30") +
	scale_x_continuous('HBond Torsional Angle (degrees)') +
	scale_y_continuous('Feature Density') +
	theme(legend.position=c(.58,.35)) +
	theme(legend.justification=c("left", "top"))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

dens <- estimate_density_1d_wrap(f,
  c("sample_source", "acc_hybrid"), "proton_chi", xlim=c(-pi, pi),
	adjust=adjust_hybrid)

plot_id = "hxl_as_donor_proton_chi_by_hybrid"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x*180/pi, y=y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	facet_wrap( ~ acc_hybrid ) +
	ggtitle("Proton Torsional Angle for Hydroxyl (s,t): (CA, CB, OG, HG) SeqSep > 5, BFact < 30") +
	scale_x_continuous('HBond Torsional Angle (degrees)') +
	scale_y_continuous('Feature Density') +
	theme(legend.position=c(.58,.35)) +
	theme(legend.justification=c("left", "top"))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

dens <- estimate_density_1d(f,
  c("sample_source", "acc_chem_type"), "chi_diff", adjust=adjust_diff)

plot_id = "hxl_as_donor_hbond_proton_chi_diff"
ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x*180/pi, y=y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	facet_wrap( ~ acc_chem_type ) +
	ggtitle("HBond Torsion Angle - Proton Torsional Angle for Hydroxyl (s,t) Donors\n SeqSep > 5, BFact < 30") +
	scale_x_continuous('HBond Torsional Angle (degrees)') +
	scale_y_continuous('Feature Density') +
	theme(legend.position=c(.58,.35)) +
	theme(legend.justification=c("left", "top"))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

dens <- estimate_density_1d(f,
  c("sample_source", "acc_hybrid"), "chi_diff", adjust=adjust_diff_hybrid)

plot_id = "hxl_as_donor_proton_chi_diff_by_hybrid"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x*180/pi, y=y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	facet_wrap( ~ acc_hybrid ) +
	ggtitle("HBond Torsion Angle - Proton Torsional Angle for Hydroxyl (s,t) Donors\nSeqSep > 5, BFact < 30") +
	scale_x_continuous('HBond Torsional Angle (degrees)') +
	scale_y_continuous('Feature Density') +
	theme(legend.position=c(.58,.35)) +
	theme(legend.justification=c("left", "top"))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)





})) # end FeaturesAnalysis
