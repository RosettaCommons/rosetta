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
id = "SSDist_Sheet_cosBAH",
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
	hb.hbond_id,
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

CREATE TEMPORARY TABLE antiparallel_close_contact_hbonds AS
SELECT
  hb1.struct_id,
	hb1.hbond_id,
  hb1.acc_resNum AS resNum_i,
  hb1.don_resNum AS resNum_j
FROM
	ee_bb_bb_hbonds AS hb1,
  ee_bb_bb_hbonds AS hb2
WHERE
  hb1.struct_id  = hb2.struct_id AND
  hb1.acc_resNum = hb2.don_resNum AND
  hb1.don_resNum = hb2.acc_resNum;

CREATE TEMPORARY TABLE parallel_close_contact_hbonds AS
SELECT
  hb1.struct_id,
	hb1.hbond_id,
  hb1.acc_resNum AS resNum_i,
  hb1.don_resNum AS resNum_j
FROM
  ee_bb_bb_hbonds AS hb1,
  ee_bb_bb_hbonds AS hb2
WHERE
  hb1.struct_id   = hb2.struct_id AND
  (( hb1.don_resNum  = hb2.acc_resNum AND hb1.acc_resNum +2 = hb2.don_resNum     ) OR
	 ( hb1.acc_resNum  = hb2.don_resNum AND hb1.don_resNum    = hb2.acc_resNum + 2 ));

CREATE TEMPORARY TABLE antiparallel_cosBAH AS
SELECT geom.cosBAH AS cosBAH, 'antiparallel' AS strand_orientation
FROM	 hbond_geom_coords AS geom,
			 antiparallel_close_contact_hbonds AS aphbond
WHERE	 aphbond.struct_id = geom.struct_id AND aphbond.hbond_id = geom.hbond_id;


CREATE TEMPORARY TABLE parallel_cosBAH AS
SELECT geom.cosBAH AS cosBAH, 'parallel' AS strand_orientation
FROM   hbond_geom_coords as geom,
			 parallel_close_contact_hbonds AS phbond
WHERE	 geom.struct_id =  phbond.struct_id AND geom.hbond_id = phbond.hbond_id;

SELECT * FROM antiparallel_cosBAH UNION
SELECT * FROM     parallel_cosBAH;"

f <- query_sample_sources(sample_sources, sele)

dens <- estimate_density_1d(
  f, c("sample_source", "strand_orientation"),
	"cosBAH")

dp <- dens[dens$strand_orientation == "parallel",]
dap <- dens[dens$strand_orientation == "antiparallel",]

ylims <- c(0,8)
plot_id <- "SSDists_Sheet_cosBAH_parallel_anti_parallel"
p <- ggplot() + theme_bw() +
  geom_line(data=dp, aes(x=180-180/pi*acos(x), y=y, colour=sample_source)) +
  geom_indicator(data=dp, aes(indicator=counts, colour=sample_source, group=sample_source)) +
  geom_line(data=dap, aes(x=180-180/pi*acos(x), y=y, colour=sample_source)) +
  geom_indicator(data=dap, aes(indicator=counts, colour=sample_source, group=sample_source, xpos="left")) +
  ggtitle("BAH angles for beta-sheet hydrogen bonds\nnormalized for equal weight per unit distance\n(antiparallel-counts left, parallel-counts right)") +
  scale_y_continuous("FeatureDensity", ylims) +
  scale_x_continuous("BAH interior angle (degrees)",
		limits=c(80,180), breaks=c(90, 120, 150, 180))

if(nrow(sample_sources) <= 3){
  p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


#ylims <- c(.1,5)
#plot_id <- "SSDists_Sheet_O_CA_parallel_anti_parallel_zoom_3_4"
#p <- ggplot(data=dens) + theme_bw() +
#  geom_line(data=dp, aes(x=x, y=y, colour=sample_source)) +
#  geom_indicator(data=dp, aes(indicator=counts, colour=sample_source)) +
#	geom_line(data=dap, aes(x=x, y=y, colour=sample_source)) +
#  geom_indicator(data=dap, aes(indicator=counts, colour=sample_source, xpos="left")) +
#  ggtitle("O--Ca atom atom distances involving beta-sheet residues (seq. sep. > 1)\nnormalized for equal weight per unit distance") +
#  scale_y_log10("FeatureDensity", limits=ylims) +
#  scale_x_continuous(expression(paste('Atom Atom Distances (', ring(A), ')')), limits=c(3,4.5), breaks=c(3, 3.25, 3.5, 3.75, 4, 4.25, 4.5) )
#
#if(nrow(sample_sources) <= 3){
#  p <- p + theme(legend.position="bottom", legend.direction="horizontal")
#}
#
#save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
