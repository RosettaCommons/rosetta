# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

check_setup()
feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "hbond_dAMOaPBA",
filename = "scripts/analysis/instances/hbonds/dAMOaPBA.R",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self){

sele <-"
SELECT DISTINCT
	struct.tag,
	struct.struct_id || '_' || acc.resNum || '_' || don.resNum AS id,
	'' AS chain,
	'NZ' AS don_atom1, 'CE' AS don_atom2, 'CD' AS don_atom3,
	don.resNum AS don_resNum,
 	'CA' AS acc_atom1, 'C'  AS acc_atom2, 'O'  AS acc_atom3,
	acc.resNum AS acc_resNum
FROM
	structures as struct,
	hbonds AS hb,
	hbond_sites AS don, hbond_sites AS acc,
	hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb
WHERE
	hb.struct_id = struct.struct_id AND
	don.struct_id = struct.struct_id AND don.site_id = hb.don_id AND
	acc.struct_id = struct.struct_id AND acc.site_id = hb.acc_id AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30 AND
	acc.HBChemType == 'hbacc_PBA' AND don.HBChemType == 'hbdon_AMO' AND
	ABS(don.resNum - acc.resNum) > 10
ORDER BY RANDOM()
LIMIT 15;"

f <- query_sample_sources(sample_sources, sele)

if(nrow(f) == 0){
	cat("WARNING: Query returned no rows. Skipping rest of features analysis.\n")
	return()
}


ss_ids <- as.character(unique(f$sample_source))


# f:
#
#           sample_source, tag, id, chain, CA, C, O, acc_resNum, don_resNum
# instance1
#   ...
#

# instance_atoms:
#
#                  sample_source, tag, id, chain, atom
# instance1, atom1
# instance1, atom2
#   ...



don_atoms <-
	melt(f,
		id.vars=c("id", "sample_source", "tag", "chain", "don_resNum"),
		measure.vars=c("don_atom1", "don_atom2", "don_atom3"),
		variable_name = "atom_name")
names(don_atoms)[5] <- "resNum"
names(don_atoms)[7] <- "atom"

acc_atoms <-
	melt(f,
		id.vars=c("id", "sample_source", "tag", "chain", "acc_resNum"),
		measure.vars=c("acc_atom1", "acc_atom2", "acc_atom3"),
		variable_name = "atom_name")
names(acc_atoms)[5] <- "resNum"
names(acc_atoms)[7] <- "atom"

instance_atoms <- rbind(don_atoms, acc_atoms)

instances_id <- "hbond_dAMOaPBA"
prepare_feature_instances(instances_id, sample_sources, instance_atoms)

})) # end FeaturesAnalysis
