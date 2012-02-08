-- This query will return a table with the average rotamer recovery
-- for all residues with sidechain degrees of freedom and maximum
-- sidechain B-factor less than 30, where the criteria for what counts
-- as recovered is determined by the RRComparer used.

.separator ', '

SELECT
	AVG(rr.recovered) AS avg_recovery
FROM
	residues AS res,
	residue_pdb_confidence AS res_conf,
	rotamer_recovery AS rr
WHERE
	res.name3 != 'ALA' AND res.name3 != 'GLY' AND
	res_conf.struct_id = res.struct_id AND
	res_conf.residue_number = res.resNum AND
	res_conf.max_sc_temperature < 30 AND
	rr.struct_id = res.struct_id AND
	rr.resNum = res.resNum;

