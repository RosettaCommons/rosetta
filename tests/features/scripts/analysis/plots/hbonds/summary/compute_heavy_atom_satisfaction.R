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
id = "compute_heavy_atom_satisfaction",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <- "
DROP TABLE IF EXISTS hbond_heavy_sites;
CREATE TABLE IF NOT EXISTS hbond_heavy_sites (
	site_id INTEGER PRIMARY KEY AUTOINCREMENT,
	struct_id BLOB,
	residue_number INTEGER,
	HBChemType TEXT,
	is_donor INTEGER,
	res_type TEXT,
	satisfied_label TEXT,
	satisfied INTEGER,
	buried_label TEXT,
	buried INTEGER);

INSERT INTO hbond_heavy_sites
SELECT
	NULL,
	site.struct_id,
	site.resNum,
	site.HBChemType,
	site.is_donor,
	res.name3 AS res_type,
	CASE WHEN env.num_hbonds == 0 THEN 'acc' ELSE 'ACC' END AS satisfied_label,
	env.num_hbonds != 0 AS satisfied,
	CASE WHEN env.sasa_r100 == 0 THEN 'Buried' ELSE 'Exposed' END AS buried_label,
	env.sasa_r100 == 0 AS buried
FROM
	residues AS res,
	hbond_sites AS site,
	hbond_site_environment AS env
WHERE
	site.struct_id = res.struct_id AND
	site.resNum = res.resNum AND
	env.struct_id = site.struct_id AND
	env.site_id = site.site_id AND
	site.is_donor = 0;

INSERT INTO hbond_heavy_sites
SELECT
	NULL,
	site.struct_id,
	site.resNum,
	site.HBChemType,
	site.is_donor,
	res.name3 AS res_type,
	CASE WHEN SUM(env.num_hbonds) == 0 THEN 'don' ELSE 'DON' END AS satisfied_label,
	SUM(env.num_hbonds) != 0 AS satisfied,
	CASE WHEN SUM(env.sasa_r100) == 0 THEN 'Buried' ELSE 'Exposed' END AS buried_label,
	SUM(env.sasa_r100) == 0 AS buried
FROM
	residues AS res,
	hbond_sites AS site,
	hbond_site_environment AS env,
	hbond_site_atoms AS atms
WHERE
	site.struct_id = res.struct_id AND
	site.resNum = res.resNum AND
	env.struct_id = site.struct_id AND
	env.site_id = site.site_id AND
	site.is_donor = 1 AND
	atms.struct_id = site.struct_id AND
	atms.site_id = site.site_id
GROUP BY
	atms.base_x, atms.base_y, atms.base_z;

DROP INDEX IF EXISTS hbond_heavy_sites_struct_id_residue_number;
CREATE INDEX hbond_heavy_sites_struct_id_residue_number ON
	hbond_heavy_sites ( struct_id, residue_number );

SELECT count(*) FROM hbond_heavy_sites;"
heavy_sites <- query_sample_sources(sample_sources, sele);

})) # end FeaturesAnalysis
