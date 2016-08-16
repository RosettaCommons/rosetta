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
id = "tyr",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("StructureFeatures", "ResidueFeatures", "PdbDataFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT DISTINCT
	struct.tag,
	struct.struct_id || '_' || res.resNum AS id,
	'' AS chain,
	res.resNum AS resNum,
	'' AS atom
FROM
	structures AS struct,
	residues AS res,
	residue_pdb_confidence AS res_conf
WHERE
	res.struct_id = struct.struct_id AND
	res.name3 = 'TYR' AND
	res_conf.struct_id = res.struct_id AND
	res_conf.residue_number = res.resNum AND
	res_conf.max_temperature < 30
ORDER BY RANDOM()
LIMIT 15;"

f <- query_sample_sources(sample_sources, sele)

if(nrow(f) == 0){
	cat("WARNING: Query returned no rows. Skipping rest of features analysis.\n")
	return()
}

instances_id <- "residues_tyr"
prepare_feature_instances(instances_id, sample_sources, f)

})) # end FeaturesAnalysis
