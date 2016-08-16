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
id = "SSdists_Sheet_O_HA",
author = "Matthew O'Meara",
brief_description = "",
long_description = "

        RESi
  ---ACC1-DON2 -----<
 /    |    |
|    HB1  HB2             ANTI-PARALLEL HBonding
 \    |    |
  ---DON1-ACC2------>
       RESj


          ||              ||
          ||      Bres1   Ares1
  ------- Ri+1 -- Ri ---- Ri-1 -------<
 /               /||    /               CLOSE Carbon HBond
|            O_Ca || Ca_O               Geometry in anti-parallel sheets
 \           /    || /
  ------- Rj-1 -- Rj ---- Rj+1 -------->
          Bres2   Ares2   ||
          ||              ||


         |               |
         |      RESi     |
----Ca---|---DON1--ACC2--|--Ca-------->
     \      /      /  \      \
      \   HB1     /   HB2     \         PARALLEL HBonding
       \  /      /      \      \
--DON--ACC1--|--Ca---|---DON2--ACC---->
   RESj-1    |  RESj |     RESj+1
             |       |
",

feature_reporter_dependencies = c("ResidueFeatures", "HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
CREATE TEMPORARY TABLE ee_bb_bb_hbonds AS
SELECT
  hb.struct_id,
  acc_site.resNum AS acc_resNum,
  don_site.resNum AS don_resNum
FROM
  hbonds AS hb,
  hbond_sites AS acc_site,
  hbond_sites AS don_site,
  residue_secondary_structure AS r1ss,
  residue_secondary_structure AS r2ss
WHERE
  acc_site.HBChemType = 'hbacc_PBA' AND
  acc_site.struct_id = hb.struct_id AND
  acc_site.site_id   = hb.acc_id AND
  don_site.site_id   = hb.don_id AND
  don_site.HBChemType = 'hbdon_PBA' AND
  don_site.struct_id = acc_site.struct_id AND
  r1ss.struct_id     = acc_site.struct_id AND
  r1ss.resNum        = acc_site.resNum AND
  r1ss.dssp          = 'E' AND
  r2ss.struct_id     = don_site.struct_id AND
  r2ss.resNum        = don_site.resNum AND
  r2ss.dssp          = 'E';

CREATE UNIQUE INDEX IF NOT EXISTS
	residue_type_atom_residue_type_set_name_residue_type_name_atom_name ON
	residue_type_atom ( residue_type_set_name, residue_type_name, atom_name );

CREATE TEMPORARY TABLE O_HA_atom_coords AS
SELECT
	res.struct_id AS struct_id, res.resNum AS resNum,
	O.x AS Ox, O.y AS Oy, O.z AS Oz,
	HA.x AS HAx, HA.y AS HAy, HA.z AS HAz
FROM
	residues AS res,
	residue_type_atom AS O_type, residue_type_atom AS HA_type,
	residue_atom_coords AS O, residue_atom_coords AS HA
WHERE
	O_type.residue_type_set_name = 'fa_standard' AND
	O_type.residue_type_name = res.res_type AND O_type.atom_name = ' O  ' AND
	HA_type.residue_type_set_name = 'fa_standard' AND
	HA_type.residue_type_name = res.res_type AND HA_type.atom_name = ' HA ' AND
	O.struct_id = res.struct_id AND O.seqpos = res.resNum AND
	O.atomno = O_type.atom_index AND
	HA.struct_id = res.struct_id AND HA.seqpos = res.resNum AND
	HA.atomno = HA_type.atom_index;

CREATE UNIQUE INDEX O_HA_atom_coords_struct_id_resNum ON
	O_HA_atom_coords (struct_id, resNum);

CREATE TEMPORARY TABLE antiparallel_close_contact_residue_pairs AS
SELECT
  hb1.struct_id,
  hb1.acc_resNum AS resNum_i,
  hb1.don_resNum AS resNum_j
FROM
	ee_bb_bb_hbonds AS hb1,
  ee_bb_bb_hbonds AS hb2
WHERE
  hb1.struct_id  = hb2.struct_id AND
  hb1.acc_resNum = hb2.don_resNum AND
  hb1.don_resNum = hb2.acc_resNum AND
  hb1.acc_resNum < hb1.don_resNum;

CREATE TEMPORARY TABLE antiparallel_HA_O AS
SELECT (O.Ox-HA.HAx)*(O.Ox-HA.HAx) +
       (O.Oy-HA.HAy)*(O.Oy-HA.HAy) +
	     (O.Oz-HA.HAz)*(O.Oz-HA.HAz) AS dist_sq,
	     'antiparallel' AS strand_orientation
FROM   O_HA_atom_coords AS O, O_HA_atom_coords AS HA,
       antiparallel_close_contact_residue_pairs AS atpairs
WHERE  HA.struct_id = atpairs.struct_id AND HA.resNum = atpairs.resNum_i - 1 AND
	     O.struct_id = atpairs.struct_id AND O.resNum = atpairs.resNum_j;

CREATE TEMPORARY TABLE antiparallel_O_HA AS
SELECT (O.Ox-HA.HAx)*(O.Ox-HA.HAx) +
       (O.Oy-HA.HAy)*(O.Oy-HA.HAy) +
	     (O.Oz-HA.HAz)*(O.Oz-HA.HAz) AS dist_sq,
	     'antiparallel' AS strand_orientation
FROM   O_HA_atom_coords AS O, O_HA_atom_coords AS HA,
       antiparallel_close_contact_residue_pairs AS atpairs
WHERE  O.struct_id = atpairs.struct_id AND O.resNum = atpairs.resNum_i AND
	     HA.struct_id = atpairs.struct_id AND HA.resNum = atpairs.resNum_j - 1;

CREATE TEMPORARY TABLE parallel_close_contact_residue_pairs AS
SELECT
  hb1.struct_id,
  hb1.acc_resNum AS resNum_i,
  hb1.don_resNum AS resNum_j
FROM
  ee_bb_bb_hbonds AS hb1,
  ee_bb_bb_hbonds AS hb2
WHERE
  hb1.struct_id   = hb2.struct_id AND
  hb1.don_resNum  = hb2.acc_resNum AND
  hb1.acc_resNum +2 = hb2.don_resNum;

CREATE TEMPORARY TABLE parallel_O_HA AS
SELECT
	(O.Ox-HA.HAx)*(O.Ox-HA.HAx) +
	(O.Oy-HA.HAy)*(O.Oy-HA.HAy) +
	(O.Oz-HA.HAz)*(O.Oz-HA.HAz) AS dist_sq,
	'parallel' AS strand_orientation
FROM
	O_HA_atom_coords AS O, O_HA_atom_coords AS HA,
	parallel_close_contact_residue_pairs AS apairs
WHERE
	O.struct_id = apairs.struct_id AND O.resNum = apairs.resNum_i AND
	HA.struct_id = apairs.struct_id AND HA.resNum = apairs.resNum_j;

SELECT * FROM antiparallel_HA_O UNION
SELECT * FROM antiparallel_O_HA UNION
SELECT * FROM parallel_O_HA;"

f <- query_sample_sources(sample_sources, sele)
f$dist <- sqrt(f$dist_sq)

dens <- estimate_density_1d(
  f, c("sample_source", "strand_orientation"),
	"dist", weight_fun = radial_3d_normalization, n_pts=1000)

dp <- dens[dens$strand_orientation == "parallel",]
dap <- dens[dens$strand_orientation == "antiparallel",]

plot_id <- "SSDists_Sheet_O_HA_parallel_anti_parallel"
p <- ggplot() + theme_bw() +
	geom_line(data=dp, aes(x=x, y=y, colour=sample_source)) +
	geom_indicator(data=dp, aes(indicator=counts, colour=sample_source, group=sample_source)) +
	geom_line(data=dap, aes(x=x, y=y, colour=sample_source)) +
	geom_indicator(data=dap, aes(indicator=counts, colour=sample_source, group=sample_source, xpos="left")) +
	ggtitle("O--HA atom atom distances involving beta-sheet residues\nnormalized for equal weight per unit distance\n(antiparallel-counts left, parallel-counts right)") +
	scale_y_continuous("Feature Density") +
	scale_x_continuous(expression(paste('O--HA Distances (', ring(A), ')')),
		limit=c(1.9, 6), breaks=c(2, 3, 4, 5, 6))

if(nrow(sample_sources) <= 3){
  p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
