####
#
# This script creates a table in a database that specifies the symmetry
# definitions for each structure
# 
# To use, use the task operation <ReadResfileFromDB db_name=<resfiles_db>/>
#
##

import glob
import sqlite3
import subprocess

data_dir = "top8000_chains_eds_70"
rosetta_inputs_db = "rosetta_inputs.db3"
rosetta = "/home/momeara/rosetta/mini"

conn = sqlite3.connect(rosetta_inputs_db)
conn.execute("CREATE TABLE IF NOT EXISTS symmetry_definitions (tag TEXT, symmetry_definition TEXT, PRIMARY KEY(tag));")

structure_list = open("symm_todo")
for fname in structure_list.read().split("\n"):
    p = subprocess.Popen(
        args=["perl", rosetta + "/src/apps/public/symmetry/make_symmdef_file.pl",
              "-m", "CRYST", "-p", fname, "-s", "R 3 :H", "-r", "12"],
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE)
    (symmdef, err) = p.communicate()
    print "fname:", fname
    print err
    conn.execute("DELETE FROM symmetry_definitions WHERE tag = ?;", (fname,))
    conn.execute("INSERT INTO symmetry_definitions VALUES (?,?);",
                 (fname, symmdef))

conn.commit()
conn.close()
