#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/tests/Test_Module_Structure.py
## @brief  PyRosetta Toolkit Test script
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Python Imports
import os
import sys
rel = "../"
sys.path.append(os.path.abspath(rel))

from modules.SQLPDB import SQLPDB
from modules.SQLPDB import PDB_database

print "Testing SQLPDB"

outdir = os.path.dirname(os.path.abspath(__file__))+"/outputs"
pdb_path = os.path.dirname(os.path.abspath(__file__))+"/data/2j88.pdb"

if os.path.exists(outdir+"/test_db.db"):
    os.remove(outdir+"/test_db.db")
    
sql = SQLPDB("2J88", "scfv_antibody_renumbered", 1, False, outdir+"/test_db.db")

print "Reading 2j88 into the database"
sql.read_pdb_into_database_flat(pdb_path, "L")

print "Fetching 2j88 and reading into database"
sql.set_structID(2)
sql.set_modelID("antibody_fetched")
sql.fetch_and_read_pdb_into_database("2j88")

print "Testing PDBDatabase"
pdb_database = PDB_database(outdir+"/test_db.db")
pdb_database.query_strucID("pdb", 1)
pdb_database.save_cur_as_pdb(outdir+"/2j88_L_db_test_out.pdb")

pdb_database._reset_cursor()
pdb_database.query_all()
pdb_database.save_cur_as_pdb(outdir+"/all_db_test_out.pdb")
