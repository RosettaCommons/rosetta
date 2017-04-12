# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

# check_setup()


library(ggplot2)
library(plyr)
library("scales")
library(grid)
library(viridis)

feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "rama_zoom_helix_stbb_hbonds",
author = "Andrew Leaver-Fay",
brief_description = "Scatter the second sample source on the phi/psi distribution of the first sample source",
feature_reporter_dependencies = c("ResidueFeatures", "ProteinBackboneTorsionAngleFeatures", "ResidueSecondaryStructureFeatures", "PdbDataFeatures", "ProteinResidueConformationFeatures","HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

# TEMP ! THIS WORKS! sele <-"
# TEMP ! THIS WORKS! CREATE TEMPORARY TABLE bif_hbonds_to_bb3 AS SELECT
# TEMP ! THIS WORKS!   hbond.struct_id AS struct_id,
# TEMP ! THIS WORKS!   hbond.don_id AS don_id,
# TEMP ! THIS WORKS!   hbond.acc_id AS acc_id,
# TEMP ! THIS WORKS!   hbond.energy AS energy
# TEMP ! THIS WORKS! FROM
# TEMP ! THIS WORKS!   hbonds as hbond,
# TEMP ! THIS WORKS!   hbond_sites AS acc_site,
# TEMP ! THIS WORKS!   hbond_site_environment as acc_env
# TEMP ! THIS WORKS! WHERE
# TEMP ! THIS WORKS!   hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
# TEMP ! THIS WORKS!   acc_site.HBChemType = 'hbacc_PBA' AND
# TEMP ! THIS WORKS!   acc_site.struct_id = acc_env.struct_id AND acc_site.site_id = acc_env.site_id AND
# TEMP ! THIS WORKS!   acc_env.num_hbonds >= 2;
# TEMP ! THIS WORKS! 
# TEMP ! THIS WORKS! SELECT
# TEMP ! THIS WORKS! 	CASE WHEN ss1.dssp='H' AND ss2.dssp='H' AND ss3.dssp='H' AND ss4.dssp='H'
# TEMP ! THIS WORKS!     THEN 1
# TEMP ! THIS WORKS!     ELSE 0
# TEMP ! THIS WORKS!     END
# TEMP ! THIS WORKS!   AS all_helix,
# TEMP ! THIS WORKS!   substr( don_site.resType, 1, 3 ) as don_res_type,
# TEMP ! THIS WORKS!   ( acc_site.resNum - don_site.resNum ) as seq_sep
# TEMP ! THIS WORKS! FROM
# TEMP ! THIS WORKS!   bif_hbonds_to_bb3 AS hbond,
# TEMP ! THIS WORKS!   bif_hbonds_to_bb3 AS hbond2,
# TEMP ! THIS WORKS!   hbond_sites AS acc_site,
# TEMP ! THIS WORKS!   hbond_sites AS don_site,
# TEMP ! THIS WORKS!   hbond_sites AS don_site2,
# TEMP ! THIS WORKS! 	protein_backbone_torsion_angles AS bb,
# TEMP ! THIS WORKS! 	residue_secondary_structure AS ss1,
# TEMP ! THIS WORKS! 	residue_secondary_structure AS ss2,
# TEMP ! THIS WORKS! 	residue_secondary_structure AS ss3,
# TEMP ! THIS WORKS! 	residue_secondary_structure AS ss4
# TEMP ! THIS WORKS! WHERE
# TEMP ! THIS WORKS!   hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
# TEMP ! THIS WORKS!   hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
# TEMP ! THIS WORKS!   hbond2.struct_id = hbond.struct_id AND hbond2.acc_id = hbond.acc_id AND
# TEMP ! THIS WORKS!   hbond2.struct_id = don_site2.struct_id AND hbond2.don_id = don_site2.site_id AND
# TEMP ! THIS WORKS!   don_site.HBChemType = 'hbdon_HXL' AND
# TEMP ! THIS WORKS!   don_site2.HBChemType = 'hbdon_PBA' AND
# TEMP ! THIS WORKS! 	( acc_site.resNum - don_site.resNum == -3 OR acc_site.resNum - don_site.resNum == -4 ) AND
# TEMP ! THIS WORKS! 	bb.struct_id = acc_site.struct_id AND bb.resNum == acc_site.resNum AND
# TEMP ! THIS WORKS! 	ss1.struct_id = acc_site.struct_id AND ss1.resNum == acc_site.resNum AND
# TEMP ! THIS WORKS! 	ss2.struct_id = acc_site.struct_id AND ss2.resNum == acc_site.resNum+1 AND
# TEMP ! THIS WORKS! 	ss3.struct_id = acc_site.struct_id AND ss3.resNum == acc_site.resNum+2 AND
# TEMP ! THIS WORKS! 	ss4.struct_id = acc_site.struct_id AND ss4.resNum == acc_site.resNum+3
# TEMP ! THIS WORKS! ;"
# TEMP ! THIS WORKS! ss2 = data.frame(sample_sources[2,])
# TEMP ! THIS WORKS! g <- query_sample_sources( ss2, sele)
# TEMP ! THIS WORKS! 
# TEMP ! THIS WORKS! print(table(g))

sele <-"
SELECT
  substr( don_site.resType, 1, 3 ) as don_res_type
FROM
  hbond_sites AS don_site,
	residue_secondary_structure AS ss
WHERE
  don_site.HBChemType = 'hbdon_HXL' AND
	ss.struct_id = don_site.struct_id AND ss.resNum == don_site.resNum AND
	ss.dssp = 'H'
;"
ss2 = data.frame(sample_sources[2,])
g <- query_sample_sources( ss2, sele)
print(table(g))

})) # end FeaturesAnalysis
