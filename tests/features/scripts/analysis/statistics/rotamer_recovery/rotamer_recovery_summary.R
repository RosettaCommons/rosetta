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
CREATE TABLE IF NOT EXISTS nchi (
        name3 TEXT,
        nchi INTEGER,
        PRIMARY KEY (name3));
INSERT OR IGNORE INTO nchi VALUES('ARG', 4);
INSERT OR IGNORE INTO nchi VALUES('LYS', 4);
INSERT OR IGNORE INTO nchi VALUES('MET', 3);
INSERT OR IGNORE INTO nchi VALUES('GLN', 3);
INSERT OR IGNORE INTO nchi VALUES('GLU', 3);
INSERT OR IGNORE INTO nchi VALUES('TYR', 2);
INSERT OR IGNORE INTO nchi VALUES('ILE', 2);
INSERT OR IGNORE INTO nchi VALUES('ASP', 2);
INSERT OR IGNORE INTO nchi VALUES('TRP', 2);
INSERT OR IGNORE INTO nchi VALUES('PHE', 2);
INSERT OR IGNORE INTO nchi VALUES('HIS', 2);
INSERT OR IGNORE INTO nchi VALUES('ASN', 2);
INSERT OR IGNORE INTO nchi VALUES('THR', 1);
INSERT OR IGNORE INTO nchi VALUES('SER', 1);
INSERT OR IGNORE INTO nchi VALUES('PRO', 1);
INSERT OR IGNORE INTO nchi VALUES('CYS', 1);
INSERT OR IGNORE INTO nchi VALUES('VAL', 1);
INSERT OR IGNORE INTO nchi VALUES('LEU', 1);
INSERT OR IGNORE INTO nchi VALUES('ALA', 0);
INSERT OR IGNORE INTO nchi VALUES('GLY', 0);

CREATE TEMPORARY TABLE IF NOT EXISTS chi_recovered AS
SELECT
	res.name3 AS res_type,
	nchi.nchi,
	nchi.nchi-rr.divergence > 0 AS first,
	nchi.nchi-rr.divergence > 1 AS second,
	nchi.nchi-rr.divergence > 2 AS third,
	nchi.nchi-rr.divergence > 3 AS fourth
FROM
	residues AS res,
	rotamer_recovery AS rr,
	residue_pdb_confidence AS res_conf,
	nchi
WHERE
	rr.struct_id = res.struct_id AND rr.resNum = res.resNum AND
	nchi.name3 = res.name3 AND
	res_conf.struct_id = res.struct_id AND
	res_conf.residue_number = res.resNum AND
	res_conf.max_sc_temperature < 30 AND
	nchi.nchi > 0;


SELECT
	rr.res_type,
	count(*) AS count,
	CASE rr.nchi > 0 WHEN 1 THEN AVG(first)*100 ELSE NULL END AS first,
	CASE rr.nchi > 1 WHEN 1 THEN AVG(second)*100 ELSE NULL END AS second,
	CASE rr.nchi > 2 WHEN 1 THEN AVG(third)*100 ELSE NULL END AS third,
	CASE rr.nchi > 3 WHEN 1 THEN AVG(fourth)*100 ELSE NULL END AS fourth
FROM
	chi_recovered AS rr
GROUP BY res_type;"


f <- query_sample_sources(sample_sources, sele)

#reorder the residue types by functional group similarity
f$res_type <- factor(f$res_type,
	levels=c("TRP", "PHE", "GLN", "GLU", "TYR", "SER", "ARG", "HIS", "LEU", "MET", "CYS", "THR", "ASP", "ILE", "VAL", "LYS", "ASN", "PRO"),
	labels=c("TRP", "PHE", "GLN", "GLU", "TYR", "SER", "ARG", "HIS", "LEU", "MET", "CYS", "THR", "ASP", "ILE", "VAL", "LYS", "ASN", "PRO"))

#sort the table first by res_type then by sample_source
f <- f[order(f$res_type, f$sample_source),]

#make recovery values two digit numbers
f[,c(4,5,6,7)] <- round(colwise(as.numeric)(f[,c(4,5,6,7)]),2)


table_id <- "rotamer_recovery_summary"
table_title <- "Recovery of Chi Angle Rotamers after applying MinPack by Residue type, BFactor < 30\n"
save_tables(self,
	f, table_id,
	sample_sources, output_dir, output_formats,
	caption=table_title, caption.placement="top")


})) # end FeaturesAnalysis
