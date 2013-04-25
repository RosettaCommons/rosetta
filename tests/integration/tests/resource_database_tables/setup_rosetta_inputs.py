
####
#
# This script creates a table in a database that specifies as a
# resfile the relevant chain in each pdb.
#
# Run this script
#
#    python setup_resfiles.py --data_dir <path_to_pdb_files> --rosetta_inputs_db <rosetta_inputs.db3>
#
#    In the data directory should be pdb files named like ????FH_?.pdb
#    where the last "?" indicates which chain should be used.
#
# To use the resulting rosetta inputs database, use the task operation
#
#    <ReadResfileFromDB
#        name=relevant_chains
#        db=<rosetta_inputs.db3>
#        table=resfiles/>
#
##

import sqlite3

rosetta_inputs_db = "rosetta_inputs.db3"
conn = sqlite3.connect(rosetta_inputs_db)

conn.execute("""
CREATE TABLE starting_structures (
	pdb_code TEXT,
	pdb_file TEXT,
	PRIMARY KEY ( pdb_code ) );
""")

f = open("../../tests/make_symmdef_file/inputs/1xu1FH_D.pdb")
conn.execute("INSERT INTO starting_structures VALUES ( ?, ? );", [ "1xu1FH_D", f.read() ])
f.close()


conn.execute("""
CREATE TABLE symmetry_definitions (
	pdb_code TEXT,
	symmdef_file TEXT,
	PRIMARY KEY ( pdb_code ) );
""")

f = open("../../tests/repack_with_elec_dens/inputs/1xu1FH_D.symm")
conn.execute("INSERT INTO symmetry_definitions VALUES ( ?, ? );", [ "1xu1FH_D", f.read() ])
f.close()

conn.executescript("""
CREATE TABLE resource_definitions (
	resource_tag TEXT,
	locator_tag TEXT,
	locator_id TEXT,
	loader_type TEXT,
	options_tag TEXT);

INSERT INTO resource_definitions VALUES (
	'1xu1FH_D_starting_struct',
	'starting_structures',
	'1xu1FH_D',
	'PoseFromPDB',
	'pose_options');

INSERT INTO resource_definitions VALUES (
	'1xu1FH_D_symmdef',
	'symmetry_definitions',
	'1xu1FH_D',
	'SymmData',
	NULL);
""")




conn.executescript("""
CREATE TABLE jobs (
	job_name TEXT,
	data_type TEXT,
	key TEXT,
	value TEXT);

INSERT INTO jobs VALUES (
	'1xu1FH_D',
	'Resource',
	'startstruct',
	'1xu1FH_D_starting_struct');

INSERT INTO jobs VALUES (
	'1xu1FH_D',
	'Resource',
	'symmdef',
	'1xu1FH_D_symmdef');

INSERT INTO jobs VALUES (
	'1xu1FH_D',
	'Resource',
	'scores_db',
	'scores_db');
""")

conn.commit()
conn.close()

