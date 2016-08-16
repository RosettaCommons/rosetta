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
id = "SSDists_Sheet_O_N",
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

sele <-"
CREATE TEMPORARY TABLE ee_atpair_dists AS
SELECT
  dist.O_N_dist, dist.N_O_dist
FROM
  protein_backbone_atom_atom_pairs as dist,
  residue_secondary_structure as r1ss,
  residue_secondary_structure as r2ss
WHERE
  dist.struct_id = r1ss.struct_id AND
  dist.struct_id = r2ss.struct_id AND
  dist.resNum1 = r1ss.resNum AND
  dist.resNum2 = r2ss.resNum AND
  dist.resNum1 + 1 != dist.resNum2 AND
  r1ss.dssp = 'E' AND r2ss.dssp = 'E'
LIMIT 500000;

SELECT O_N_dist as dist FROM ee_atpair_dists UNION
SELECT N_O_dist as dist FROM ee_atpair_dists;"
f <- query_sample_sources(sample_sources, sele)

dens <- estimate_density_1d(
  f, c("sample_source"),
  "dist", weight_fun = radial_3d_normalization)

plot_id <- "SSDists_Sheet_O_N_seqsep_gt1"
p <- ggplot(data=dens) + theme_bw() +
  geom_line(aes(x=x, y=y, colour=sample_source)) +
  geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
  ggtitle("O--N atom atom distances involving beta-sheet residues (seq. sep. > 1)\nnormalized for equal weight per unit distance") +
  scale_y_log10("FeatureDensity", limits=c(1e-3,1e0)) +
  scale_x_continuous(expression(paste('Atom Atom Distances (', ring(A), ')')), limits=c(1.5,7), breaks=2:7)

if(nrow(sample_sources) <= 3){
  p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)



lims <- c(10^-1.5, 10^-.5)
plot_id <- "SSDists_Sheet_O_N_seqsep_gt1_zoom_4_6p5"
p <- ggplot(data=dens) + theme_bw() +
  geom_line(aes(x=x, y=y, colour=sample_source)) +
  geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
  ggtitle("O--N atom atom distances involving beta-sheet residues (seq. sep. > 1)\nnormalized for equal weight per unit distance") +
  scale_y_log10("FeatureDensity", limits=lims) +
  scale_x_continuous(expression(paste('Atom Atom Distances (', ring(A), ')')), limits=c(4,6.5), breaks=c(4, 4.5, 5, 5.6, 6, 6.5) )

if(nrow(sample_sources) <= 3){
  p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
