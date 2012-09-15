
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

conn.commit()
conn.close()

