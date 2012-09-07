.header on
.separator ' '
SELECT
	score_type.score_type_name,
	structure_score.score_value
FROM
	score_types AS score_type,
	structure_scores AS structure_score
WHERE
	structure_score.score_type_id = score_type.score_type_id;