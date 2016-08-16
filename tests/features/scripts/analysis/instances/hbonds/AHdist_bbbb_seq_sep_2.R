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
id = "AHdist_bbbb_seq_sep_2",
author = "Matthew O'Meara",
brief_description = "",
long_description = "
Backbone-backbone hbonds with sequence separation
(don.resNum-acc.resNum) 2, are predominately in turns and have the
residue in between with phi~-80 and psi ~(50, 80). It is a little
strained but it there are plenty of examples in the natives.",
feature_reporter_dependencies = c("StructureFeatures", "HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <-"
SELECT
	struct.tag,
	acc.site_id AS id,
	'' AS chain,
	acc.resNum,
	'CA' AS CA, 'C' AS C, 'O' AS O
FROM
	structures as struct,
	hbonds AS hbond,
	hbond_sites AS don, hbond_sites AS acc
WHERE
	hbond.struct_id = struct.struct_id AND
	don.struct_id = struct.struct_id AND don.site_id = hbond.don_id AND
	acc.struct_id = struct.struct_id AND acc.site_id = hbond.acc_id AND
	acc.HBChemType == 'hbacc_PBA' AND don.HBChemType == 'hbdon_PBA' AND
	don.resNum - acc.resNum = 2
ORDER BY RANDOM()
LIMIT 30;"
f <- query_sample_sources(sample_sources, sele)

if(nrow(f) == 0){
	cat("WARNING: Query returned no rows. Skipping rest of features analysis.\n")
	return()
}

f <- melt(f,
	id.vars=c("sample_source", "tag", "id", "chain", "resNum"),
	measure.vars=c("CA", "C", "O"),
	variable_name = "atom")


instances_id <- "AHdist_bbbb_seq_sep_2"

prepare_feature_instances(instances_id, sample_sources, f, output_dir)

})) # end FeaturesAnalysis
