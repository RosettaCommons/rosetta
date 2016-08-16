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
id = "rotamer_recovery_by_secondary_structure",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueFeatures", "RotamerRecoveryFeatures"),
run=function(self, sample_sources, output_dir, output_formats){



sele <-"
DROP TABLE IF EXISTS nchi;

CREATE TEMPORARY TABLE nchi (
	name3 TEXT,
	nchi INTEGER,
	PRIMARY KEY (name3));

INSERT INTO nchi VALUES('ARG', 4);
INSERT INTO nchi VALUES('LYS', 4);
INSERT INTO nchi VALUES('TYR', 3);
INSERT INTO nchi VALUES('MET', 3);
INSERT INTO nchi VALUES('PRO', 3);
INSERT INTO nchi VALUES('GLN', 3);
INSERT INTO nchi VALUES('GLU', 3);
INSERT INTO nchi VALUES('ILE', 2);
INSERT INTO nchi VALUES('ASP', 2);
INSERT INTO nchi VALUES('TRP', 2);
INSERT INTO nchi VALUES('PHE', 2);
INSERT INTO nchi VALUES('CYS', 2);
INSERT INTO nchi VALUES('HIS', 2);
INSERT INTO nchi VALUES('ASN', 2);
INSERT INTO nchi VALUES('THR', 2);
INSERT INTO nchi VALUES('SER', 2);
INSERT INTO nchi VALUES('VAL', 1);
INSERT INTO nchi VALUES('LEU', 1);
INSERT INTO nchi VALUES('ALA', 0);
INSERT INTO nchi VALUES('GLY', 0);

CREATE TEMPORARY TABLE IF NOT EXISTS chi_recovered AS
SELECT
	res.name3 AS res_type,
	nchi.nchi,
	nchi.nchi-rr.divergence > 0 AS first,
	nchi.nchi-rr.divergence > 1 AS second,
	nchi.nchi-rr.divergence > 2 AS third,
	nchi.nchi-rr.divergence > 3 AS fourth,
	dssp_code.label AS dssp,
	CASE bb.sasa_r140 WHEN 0 THEN 0 ELSE 1 END AS buried
FROM
	residues AS res,
	rotamer_recovery AS rr,
	residue_secondary_structure AS ss, dssp_codes AS dssp_code,
	residue_burial AS bb,
	nchi
WHERE
	rr.struct_id = res.struct_id AND rr.resNum = res.resNum AND
	ss.struct_id = res.struct_id AND ss.resNum = res.resNum AND
	bb.struct_id = res.struct_id AND bb.resNum = res.resNum AND
	dssp_code.code = ss.dssp AND
	nchi.name3 = res.name3;

SELECT
	res_type,
	dssp,
	buried,
	AVG(first) AS first * 100,
	COUNT(*) AS count
FROM
	chi_recovered
WHERE
	nchi > 0
GROUP BY
	res_type, dssp;"

f <- query_sample_sources(sample_sources, sele)

cat("Rotamer Recovery by residue type and secondary structure\n")
print(f)


g <- f[order(f$dssp, f$res_type, f$sample_source, f$first),]
g <- cast(g, dssp + res_type ~ sample_source, value="first")
g$diff <- g$top8000_olf_r45890_111114 - g$top8000_r45890_111114
print(ascii(g))

})) # end FeaturesAnalysis
