# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

check_setup()
feature_analyses <- c(feature_analyses, new("FeatureAnalysis",
id = "rotamer_recovery_summary",
filename = "scripts/analysis/statistics/rotamer_recovery/rotamer_recovery_summary.R",
author = "Matthew O'Meara",
brief_description = "",
long_description = "
This feature analysis requires the comparison to be the RRComparerRotBins
",

feature_reporter_dependencies = c("RotamerRecoveryFeatures", "PdbDataFeatures"),
run=function(){

library(ascii)

# Rotamer recovery rate of first chi angle
sele <-"
DROP TABLE IF EXISTS nchi;

CREATE TEMPORARY TABLE nchi (
	name3 TEXT,
	nchi INTEGER,
	PRIMARY KEY (name3));

INSERT INTO nchi VALUES('ARG', 4);
INSERT INTO nchi VALUES('LYS', 4);
INSERT INTO nchi VALUES('MET', 3);
INSERT INTO nchi VALUES('GLN', 3);
INSERT INTO nchi VALUES('GLU', 3);
INSERT INTO nchi VALUES('TYR', 2);
INSERT INTO nchi VALUES('ILE', 2);
INSERT INTO nchi VALUES('ASP', 2);
INSERT INTO nchi VALUES('TRP', 2);
INSERT INTO nchi VALUES('PHE', 2);
INSERT INTO nchi VALUES('HIS', 2);
INSERT INTO nchi VALUES('ASN', 2);
INSERT INTO nchi VALUES('THR', 1);
INSERT INTO nchi VALUES('SER', 1);
INSERT INTO nchi VALUES('PRO', 1);
INSERT INTO nchi VALUES('CYS', 1);
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

f$res_type <- factor(f$res_type,
	levels=c("TRP", "PHE", "GLN", "GLU", "TYR", "SER", "ARG", "HIS", "LEU", "MET", "CYS", "THR", "ASP", "ILE", "VAL", "LYS", "ASN", "PRO"),
	labels=c("TRP", "PHE", "GLN", "GLU", "TYR", "SER", "ARG", "HIS", "LEU", "MET", "CYS", "THR", "ASP", "ILE", "VAL", "LYS", "ASN", "PRO"))

f <- f[order(f$res_type, f$sample_source),]


g <- melt(f, id.vars=c("sample_source", "res_type"), measure.vars=c("first", "second", "third", "fourth"), variable_name = "chi")
#g$value <- as.numeric(levels(g$value))[g$value]
g <- cast(g, chi + sample_source ~ res_type)



counts <- cast(f, sample_source ~ res_type, value="count")
counts <- cbind(data.frame(chi=NA), counts)

different_res_count <- F
for( name in names(counts)[c(-1, -2),]){
	different_res_count <- different_res_count || (length(unique(counts[,name])) != 1)
}
if(different_res_count){
	g <- rbind(counts, g)
} else {
	counts$sample_source <- "Residue Count"
	g <- rbind(counts[1,], g)
}

print(ascii(g, caption="Recovery of Chi Angle Rotamers after applying MinPack by Residue type, BFactor < 30", header=F))



})) # end FeatureAnalysis
