-- To reduce noise in the integration test,
--   1) Only print out a subset of the score values, and none of the score_type_id values
--   2) Remove svn information


.header on
.separator ' '

-- Score Function Weights --
SELECT
	weight.score_function_name,
	score_type.score_type_name,
	weight.weight
FROM
	score_function_weights AS weight,
	score_types AS score_type
WHERE
	weight.batch_id = score_type.batch_id AND
	weight.score_type_id = score_type.score_type_id;

-- Structure Scores --
SELECT
	batch.name AS batch,
	struct.tag,
	score_type.score_type_name AS score_type,
	structure_score.score_value AS raw_score,
	structure_score.score_value * weight.weight AS weighted_score
FROM
	batches AS batch,
	structures AS struct,
	score_types AS score_type,
	structure_scores AS structure_score,
	score_function_weights AS weight
WHERE
	score_type.batch_id = batch.batch_id AND
	struct.batch_id = batch.batch_id AND
	structure_score.struct_id = struct.struct_id AND
	structure_score.score_type_id = score_type.score_type_id AND
	weight.batch_id = score_type.batch_id AND
	weight.score_type_id = score_type.score_type_id;

-- Residue Total Scores --
SELECT
	batch.name AS batch,
	struct.tag,
	residue_score.resNum,
	residue_score.score_value
FROM
	batches AS batch,
	structures AS struct,
	residue_total_scores AS residue_score
WHERE
	struct.batch_id = batch.batch_id AND
	residue_score.struct_id = struct.struct_id;

-- 1-Body Residue Scores --
SELECT
	batch.name AS batch,
	struct.tag,
	residue_score.resNum,
	score_type.score_type_name AS score_type,
	residue_score.score_value AS raw_score,
	residue_score.score_value *  weight.weight AS weighted_score
FROM
	batches AS batch,
	structures AS struct,
	score_types AS score_type,
	residue_scores_1b AS residue_score,
	score_function_weights AS weight
WHERE
	score_type.batch_id = batch.batch_id AND
	struct.batch_id = batch.batch_id AND
	residue_score.struct_id = struct.struct_id AND
	residue_score.score_type_id = score_type.score_type_id AND
	weight.batch_id = score_type.batch_id AND
	weight.score_type_id = score_type.score_type_id;

-- Short Range 2-Body Residue Scores --
SELECT
	batch.name AS batch,
	struct.tag,
	residue_score.resNum1,
	residue_score.resNum2,
	score_type.score_type_name AS score_type,
	residue_score.score_value AS raw_score,
	residue_score.score_value *  weight.weight AS weighted_score
FROM
	batches AS batch,
	structures AS struct,
	score_types AS score_type,
	residue_scores_2b AS residue_score,
	score_function_weights AS weight
WHERE
	score_type.batch_id = batch.batch_id AND
	struct.batch_id = batch.batch_id AND
	residue_score.struct_id = struct.struct_id AND
	residue_score.score_type_id = score_type.score_type_id AND
	weight.batch_id = score_type.batch_id AND
	weight.score_type_id = score_type.score_type_id;

-- Long Range 2-Body Residue Scores --
SELECT
	batch.name AS batch,
	struct.tag,
	residue_score.resNum1,
	residue_score.resNum2,
	score_type.score_type_name AS score_type,
	residue_score.score_value AS raw_score,
	residue_score.score_value *  weight.weight AS weighted_score
FROM
	batches AS batch,
	structures AS struct,
	score_types AS score_type,
	residue_scores_lr_2b AS residue_score,
	score_function_weights AS weight
WHERE
	score_type.batch_id = batch.batch_id AND
	struct.batch_id = batch.batch_id AND
	residue_score.struct_id = struct.struct_id AND
	residue_score.score_type_id = score_type.score_type_id AND
	weight.batch_id = score_type.batch_id AND
	weight.score_type_id = score_type.score_type_id;

UPDATE protocols SET svn_url = "";
UPDATE protocols SET svn_version = "";
--DROP TABLE residue_scores_lr_2b;
--DROP TABLE residue_scores_2b;
--DROP TABLE residue_scores_1b;
--DROP TABLE residue_total_scores;
--DROP TABLE structure_scores;
--DROP TABLE score_function_weights;
--DROP TABLE score_types;


.dump