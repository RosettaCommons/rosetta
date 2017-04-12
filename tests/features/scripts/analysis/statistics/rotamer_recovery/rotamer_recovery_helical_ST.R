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
SELECT
  res.name3 AS res_type,
  rr.divergence < 20 AS recovered
FROM
	residues AS res,
  rotamer_recovery AS rr,
	residue_secondary_structure AS ss
WHERE
  rr.struct_id = res.struct_id AND rr.resNum == res.resNum AND
	res.struct_id = ss.struct_id AND res.resNum == ss.resNum AND
	ss.dssp = 'H' AND
	( res.name3 == 'SER' OR res.name3 == 'THR' )
;
"


f <- query_sample_sources(sample_sources, sele)
#print(f)

#reorder the residue types by functional group similarity
#f$res_type <- factor(f$res_type,
#	levels=c("TRP", "PHE", "GLN", "GLU", "TYR", "SER", "ARG", "HIS", "LEU", "MET", "CYS", "THR", "ASP", "ILE", "VAL", "LYS", "ASN", "PRO"),
#	labels=c("TRP", "PHE", "GLN", "GLU", "TYR", "SER", "ARG", "HIS", "LEU", "MET", "CYS", "THR", "ASP", "ILE", "VAL", "LYS", "ASN", "PRO"))

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
