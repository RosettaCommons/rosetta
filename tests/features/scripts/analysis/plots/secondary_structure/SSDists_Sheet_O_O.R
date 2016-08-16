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
id = "SSDists_Sheet_O_O",
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
feature_reporter_dependencies = c("ResidueSecondaryStructureFeatures", "HBondFeatures"),
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
  r2ss.dssp          = 'E'
ORDER BY RANDOM()
LIMIT 100000;

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

CREATE TEMPORARY TABLE antiparallel_O_O AS
SELECT dist.O_O_dist AS dist, 'antiparallel' AS strand_orientation
FROM   protein_backbone_atom_atom_pairs AS dist,
       antiparallel_close_contact_residue_pairs AS atpairs
WHERE  dist.struct_id =  atpairs.struct_id AND
       dist.resNum1 = atpairs.resNum_i AND dist.resNum2 = atpairs.resNum_j;

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

CREATE TEMPORARY TABLE parallel_O_O_minus AS
SELECT dist.O_O_dist AS dist, 'parallel_minus' AS strand_orientation
FROM   protein_backbone_atom_atom_pairs AS dist,
       parallel_close_contact_residue_pairs AS apairs
WHERE  dist.struct_id =  apairs.struct_id AND
       apairs.resNum_i < apairs.resNum_j AND
       dist.resNum1 = apairs.resNum_i AND dist.resNum2 = apairs.resNum_j - 1;

CREATE TEMPORARY TABLE parallel_O_O_plus AS
SELECT dist.O_O_dist AS dist, 'parallel_plus' AS strand_orientation
FROM   protein_backbone_atom_atom_pairs AS dist,
       parallel_close_contact_residue_pairs AS apairs
WHERE  dist.struct_id =  apairs.struct_id AND
       apairs.resNum_i < apairs.resNum_j AND
       dist.resNum1 = apairs.resNum_i AND dist.resNum2 = apairs.resNum_j + 1;

SELECT * FROM antiparallel_O_O    UNION
SELECT * FROM parallel_O_O_minus UNION
SELECT * FROM parallel_O_O_plus;"

f <- query_sample_sources(sample_sources, sele)

dens <- estimate_density_1d(
  f, c("sample_source", "strand_orientation"),
  "dist", weight_fun = radial_3d_normalization)

dpM <- dens[dens$strand_orientation == "parallel_minus",]
dpP <- dens[dens$strand_orientation == "parallel_plus",]
dap <- dens[dens$strand_orientation == "antiparallel",]

plot_id <- "SSDists_Sheet_O_O_parallel_anti_parallel"
p <- ggplot() + theme_bw() +
  geom_line(data=dpP, aes(x=x, y=y, colour=sample_source)) +
  geom_line(data=dpM, aes(x=x, y=y, colour=sample_source)) +
  geom_line(data=dap, aes(x=x, y=y, colour=sample_source)) +
  geom_indicator(data=dpP, aes(indicator=counts, colour=sample_source, group=sample_source, xpos="right")) +
  geom_indicator(data=dpM, aes(indicator=counts, colour=sample_source, group=sample_source, xpos="center")) +
  geom_indicator(data=dap, aes(indicator=counts, colour=sample_source, group=sample_source, xpos="left")) +
  ggtitle("O--O atom atom distances involving beta-sheet residues\nnormalized for equal weight per unit distance\n(antiparallel                    parallel-minus                    parallel-plus)") +
  scale_y_log10("FeatureDensity", limits=c(1e-2,5))+
  scale_x_continuous(expression(paste('Atom Atom Distances (', ring(A), ')')), limits=c(2.5,9), breaks=c(3, 4, 5, 6, 7, 8, 9))
if(nrow(sample_sources) <= 3){
  p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


# zooom in on the anti_parallel peak
lims <- c(.1, 3)
plot_id <- "SSDists_Sheet_O_O_anti_parallel_zoom"
p <- ggplot(data=dap) + theme_bw() +
  geom_line(aes(x=x, y=y, colour=sample_source)) +
  geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
  ggtitle("O--O atom atom distances involving antiparallel beta-sheet residues\nnormalized for equal weight per unit distance") +
  scale_y_log10("FeatureDensity", limits=lims) +
  scale_x_continuous(expression(paste('Atom Atom Distances (', ring(A), ')')), limits=c(3,4), breaks=c(3, 3.25, 3.5, 3.75, 4) )

if(nrow(sample_sources) <= 3){
  p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


# zooom in on the parallel-minus peak
lims <- c(.1, 5)
plot_id <- "SSDists_Sheet_O_O_parallel_minus_zoom"
p <- ggplot(data=dpM) + theme_bw() +
  geom_line(aes(x=x, y=y, colour=sample_source)) +
  geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
  ggtitle("O--O atom atom distances involving parallel beta-sheet residues i to j-1\nnormalized for equal weight per unit distance") +
  scale_y_log10("FeatureDensity", limits=lims) +
  scale_x_continuous(expression(paste('Atom Atom Distances (', ring(A), ')')), limits=c(4.5,5.5), breaks=c(4.5, 4.75, 5, 5.25, 5.5) )

if(nrow(sample_sources) <= 3){
  p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


# zooom in on the parallel-plus peak
lims <- c(.1, 1.5)
plot_id <- "SSDists_Sheet_O_O_parallel_plus_zoom"
p <- ggplot(data=dpP) + theme_bw() +
  geom_line(aes(x=x, y=y, colour=sample_source)) +
  geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
  ggtitle("O--O atom atom distances involving parallel beta-sheet residues i to j+1\nnormalized for equal weight per unit distance") +
  scale_y_log10("FeatureDensity", limits=lims) +
  scale_x_continuous(expression(paste('Atom Atom Distances (', ring(A), ')')), limits=c(7.5,8.5), breaks=c(7.5, 7.75, 8, 8.25, 8.5) )

if(nrow(sample_sources) <= 3){
  p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
