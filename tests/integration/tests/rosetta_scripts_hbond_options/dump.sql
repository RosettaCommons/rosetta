-- To reduce noise in the integration test,
--   1) Only print out a subset of the score values, and none of the score_type_id values
--   2) Remove svn information


.header on
.separator ' '

UPDATE protocols SET svn_url = "";
UPDATE protocols SET svn_version = "";

.dump