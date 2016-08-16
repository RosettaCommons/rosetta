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
id = "SSDists_Sheet_gt2bonds",
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
feature_reporter_dependencies = c("ResidueSecondaryStructureFeatures", "ProteinBackboneAtomPairFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

# new idiom: union select
sele <-"
CREATE TEMPORARY TABLE bbdists AS
SELECT
  dist.resNum1,
	dist.resNum2,
  dist.N_N_dist,   dist.N_Ca_dist,  dist.N_C_dist,   dist.N_O_dist,
  dist.Ca_N_dist,  dist.Ca_Ca_dist, dist.Ca_C_dist,  dist.Ca_O_dist,
  dist.C_N_dist,   dist.C_Ca_dist,  dist.C_C_dist,   dist.C_O_dist,
  dist.O_N_dist,   dist.O_Ca_dist,  dist.O_C_dist,   dist.O_O_dist
FROM
  protein_backbone_atom_atom_pairs AS dist,
	residue_secondary_structure AS r1ss,
	residue_secondary_structure AS r2ss
WHERE
  dist.struct_id = r1ss.struct_id AND
  dist.struct_id = r2ss.struct_id AND
  dist.resNum1 = r1ss.resNum AND
  dist.resNum2 = r2ss.resNum AND
  r1ss.dssp = 'H' AND r2ss.dssp = 'H'
LIMIT 1000000;
SELECT
  'N' AS at1, 'N' AS at2, N_N_dist AS dist FROM bbdists
UNION
SELECT
  'N' AS at1, 'CA' AS at2, N_Ca_dist AS dist FROM bbdists
UNION
SELECT
  'N' AS at1, 'CA' AS at2, Ca_N_dist AS dist FROM bbdists WHERE bbdists.resNum1 + 1 != bbdists.resNum2
UNION
SELECT
  'N' AS at1, 'C' AS at2, N_C_dist AS dist FROM bbdists
UNION
SELECT
  'N' AS at1, 'C' as at2, C_N_dist AS dist FROM bbdists WHERE bbdists.resNum1 + 1 != bbdists.resNum2
UNION
SELECT
  'N' AS at1, 'O' AS at2, N_O_dist AS dist FROM bbdists
UNION
SELECT
  'N' AS at1, 'O' AS at2, O_N_dist AS dist FROM bbdists WHERE bbdists.resNum1 + 1 != bbdists.resNum2
UNION
SELECT
  'CA' AS at1, 'CA' AS at2, Ca_Ca_dist AS dist FROM bbdists
UNION
SELECT
  'CA' AS at1, 'C' AS at2, Ca_C_dist AS dist FROM bbdists
UNION
SELECT
  'CA' AS at1, 'C' AS at2, C_Ca_dist AS dist FROM bbdists WHERE bbdists.resNum1 + 1 != bbdists.resNum2
UNION
SELECT
  'CA' AS at1, 'O' AS at2, Ca_O_dist AS dist FROM bbdists
UNION
SELECT
  'CA' AS at1, 'O' AS at2, O_Ca_dist AS dist FROM bbdists
UNION
SELECT
  'C' AS at1, 'C' AS at2, C_C_dist AS dist FROM bbdists
UNION
SELECT
  'C' AS at1, 'O' AS at2, C_O_dist AS dist FROM bbdists
UNION
SELECT
  'C' AS at1, 'O' AS at2, O_C_dist AS dist FROM bbdists
UNION
SELECT
  'O' AS at1, 'O' AS at2, O_O_dist AS dist FROM bbdists;"

f <- query_sample_sources(sample_sources, sele)


f$at1 <- factor(f$at1,
  levels = c("N", "CA", "C", "O") )
f$at2 <- factor(f$at2,
  levels = c("N", "CA", "C", "O") )
print(summary(f))

dens <- estimate_density_1d(
  f, c("sample_source", "at1", "at2"),
  "dist", weight_fun = radial_3d_normalization)

plot_id <- "SSDists_Sheet_gt2bonds"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	facet_grid(at1 ~ at2) +
	ggtitle("Backbone atom atom distances involving sheet residues (w/ > 2 bond separation)\nnormalized for equal weight per unit distance") +
	scale_y_log10("FeatureDensity", limits=c(1e-3,1e0)) +
	scale_x_continuous(expression(paste('Atom Atom Distances (', ring(A), ')')), limits=c(1.5,7), breaks=2:7)

if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


#sele <-
#SELECT
#  r1conf.secstruct, r2conf.secstruct,
#  dist.N_N_dist,   dist.N_Ca_dist,  dist.N_C_dist,   dist.N_O_dist,
#  dist.Ca_N_dist,  dist.Ca_Ca_dist, dist.Ca_C_dist,  dist.Ca_O_dist,
#  dist.C_N_dist,   dist.C_Ca_dist,  dist.C_C_dist,   dist.C_O_dist,
#  dist.O_N_dist,   dist.O_Ca_dist,  dist.O_C_dist,   dist.O_O_dist
#FROM
#  hbond_sites AS acc_site1,
#  hbond_sites AS don_site1,
#  hbonds AS hb1,
#  hbond_sites AS acc_site2,
#  hbond_sites AS don_site2,
#  hbonds AS hb2,
#  protein_backbone_atom_atom_pairs AS dist,
#	protein_residue_conformation AS r1conf,
#	protein_residue_conformation AS r2conf
#WHERE
#  acc_site1.struct_id = don_site1.struct_id AND
#  hb1.struct_id = acc_site1.struct_id AND
#  hb1.acc_id = acc_site1.site_id AND
#  hb1.don_id = don_site1.site_id AND
#  acc_site1.resNum = don_site1.resNum + 4 AND
#
#  acc_site2.struct_id = don_site2.struct_id AND
#  hb2.struct_id = acc_site2.struct_id AND
#  hb2.acc_id = acc_site2.site_id AND
#  hb2.don_id = don_site2.site_id AND
#  acc_site2.resNum = don_site2.resNum + 4 AND
#
#  don_site1.resNum = acc_site2.resnum AND
#  hb1.struct_id = hb2.struct_id AND
#
#  dist.struct_id = hb1.struct_id AND
#  ( dist.resNum1 = don_site1.resNum OR dist.resNum2 = don_site1.resNum )
#
#LIMIT 10000;

})) # end FeaturesAnalysis
