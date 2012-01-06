-- generate table of native structure tags
-- use by:
--   

CREATE TABLE IF NOT EXISTS native_tags (
	struct_id INTEGER,
	native_tag TEXT,
	FOREIGN KEY (struct_id)
		REFERENCES structures (struct_id)
		DEFERRABLE INITIALLY DEFERRED,
	PRIMARY KEY (struct_id));

INSERT INTO native_tags SELECT
	struct_id,
	substr(tag, 5, 4)
FROM
	structures;