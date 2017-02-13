SELECT batches.name, structure_scores.struct_id, score_types.score_type_name, ROUND(structure_scores.score_value,2) AS 'score'
from structure_scores
INNER JOIN score_types ON score_types.batch_id=structure_scores.batch_id AND score_types.score_type_id=structure_scores.score_type_id
INNER JOIN batches ON structure_scores.batch_id=batches.batch_id
WHERE score_type_name='total_score';

SELECT structures.tag, ROUND(SUM(case when (batches.name LIKE 'bound_mut_%' OR batches.name LIKE 'unbound_wt_%' ) then structure_scores.score_value else - structure_scores.score_value end),2) as ddG
from structure_scores
INNER JOIN score_types ON score_types.batch_id=structure_scores.batch_id AND score_types.score_type_id=structure_scores.score_type_id
INNER JOIN structures ON structures.struct_id=structure_scores.struct_id
INNER JOIN batches ON structure_scores.batch_id=batches.batch_id
WHERE score_type_name='total_score'
GROUP BY structures.tag;
