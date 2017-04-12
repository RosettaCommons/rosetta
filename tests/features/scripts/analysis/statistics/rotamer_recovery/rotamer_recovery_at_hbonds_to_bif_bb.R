# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


library(plyr)

feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "rotamer_recovery_summary",
author = "Matthew O'Meara",
brief_description = "",
long_description = "
This features analysis requires the comparison to be the RRComparerRotBins
",

feature_reporter_dependencies = c("ResidueFeatures", "RotamerRecoveryFeatures", "PdbDataFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

# Rotamer recovery rate of first chi angle
sele <-"
CREATE TEMPORARY TABLE bif_hbonds_to_bb AS SELECT
  hbond.struct_id AS struct_id,
  hbond.don_id AS don_id,
  hbond.acc_id AS acc_id,
  hbond.energy AS energy
FROM
  hbonds as hbond,
  hbond_sites AS acc_site,
  hbond_site_environment as acc_env
WHERE
  hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
  acc_site.HBChemType = 'hbacc_PBA' AND
  acc_site.struct_id = acc_env.struct_id AND acc_site.site_id = acc_env.site_id AND
  acc_env.num_hbonds >= 2;

SELECT
  don_site.resType AS res_type,
  rr.divergence < 20 AS recovered
FROM
  bif_hbonds_to_bb AS hbond,
  bif_hbonds_to_bb AS hbond2,
  hbond_sites AS acc_site,
  hbond_sites AS don_site,
  hbond_sites AS don_site2,
	protein_backbone_torsion_angles AS bb,
  rotamer_recovery AS rr
WHERE
  hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
  hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
  hbond2.struct_id = hbond.struct_id AND hbond2.acc_id = hbond.acc_id AND
  hbond2.struct_id = don_site2.struct_id AND hbond2.don_id = don_site2.site_id AND
  don_site.HBChemType != 'hbdon_PBA' AND
  don_site2.HBChemType = 'hbdon_PBA' AND
	bb.struct_id = acc_site.struct_id AND bb.resNum == acc_site.resNum AND
  rr.struct_id = hbond.struct_id AND rr.resNum == don_site.resNum
;
"

f <- query_sample_sources(sample_sources, sele)
#print(f)

#reorder the residue types by functional group similarity
f$res_type <- factor(f$res_type,
	levels=c("TRP", "PHE", "GLN", "GLU", "TYR", "SER", "ARG", "HIS", "LEU", "MET", "CYS", "THR", "ASP", "ILE", "VAL", "LYS", "ASN", "PRO"),
	labels=c("TRP", "PHE", "GLN", "GLU", "TYR", "SER", "ARG", "HIS", "LEU", "MET", "CYS", "THR", "ASP", "ILE", "VAL", "LYS", "ASN", "PRO"))

tf <- table( f )
print( "Rotamer recovery broken down by amino acid type" );
print( tf )

print( "Rotamer recovery broken down by AA and sequence separation of the donor and acceptor" );

sele <-"
SELECT
  don_site.resType AS res_type,
--  rr.divergence AS divergence,
  rr.divergence < 20 AS recovered,
--	bb.phi,
--  bb.psi,
  CASE ( acc_site.resNum - don_site.resNum )
    WHEN -1 THEN '-1' WHEN -2 THEN '-2' WHEN -3 THEN '-3' WHEN -4 THEN '-4'
    WHEN 1 THEN '1' WHEN 2 THEN '2' WHEN 3 THEN '3' WHEN 4 THEN '4'
    ELSE 'long' END AS seq_sep
FROM
  bif_hbonds_to_bb AS hbond,
  bif_hbonds_to_bb AS hbond2,
  hbond_sites AS acc_site,
  hbond_sites AS don_site,
  hbond_sites AS don_site2,
	protein_backbone_torsion_angles AS bb,
  rotamer_recovery AS rr
WHERE
  hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
  hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
  hbond2.struct_id = hbond.struct_id AND hbond2.acc_id = hbond.acc_id AND
  hbond2.struct_id = don_site2.struct_id AND hbond2.don_id = don_site2.site_id AND
  don_site.HBChemType != 'hbdon_PBA' AND
  don_site2.HBChemType = 'hbdon_PBA' AND
	bb.struct_id = acc_site.struct_id AND bb.resNum == acc_site.resNum AND
  rr.struct_id = hbond.struct_id AND rr.resNum == don_site.resNum
;
"


f <- query_sample_sources(sample_sources, sele)
#print(f)

#reorder the residue types by functional group similarity
f$res_type <- factor(f$res_type,
	levels=c("TRP", "PHE", "GLN", "GLU", "TYR", "SER", "ARG", "HIS", "LEU", "MET", "CYS", "THR", "ASP", "ILE", "VAL", "LYS", "ASN", "PRO"),
	labels=c("TRP", "PHE", "GLN", "GLU", "TYR", "SER", "ARG", "HIS", "LEU", "MET", "CYS", "THR", "ASP", "ILE", "VAL", "LYS", "ASN", "PRO"))

tf <- table( f )
print( tf )

#sort the table first by res_type then by sample_source
# f <- f[order(f$res_type, f$sample_source),]
#
# #make recovery values two digit numbers
# f[,c(4,5,6,7)] <- round(colwise(as.numeric)(f[,c(4,5,6,7)]),2)
#
#
# table_id <- "rotamer_recovery_summary"
# table_title <- "Recovery of Chi Angle Rotamers after applying MinPack by Residue type, BFactor < 30\n"
# save_tables(self,
# 	f, table_id,
# 	sample_sources, output_dir, output_formats,
# 	caption=table_title, caption.placement="top")


})) # end FeaturesAnalysis
