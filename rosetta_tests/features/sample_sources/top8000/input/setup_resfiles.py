
####
#
# This script sets up the  

# resfile the relevant chain in each pdb.
#
# Run this script
#
#    python setup_resfiles.py --data_dir <path_to_pdb_files>
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

import glob, os, sqlite3, sys
from optparse import OptionParser

def setup_resfiles(data_dir, rosetta_inputs_db):

    if not os.path.exists(data_dir): 
        print "ERROR: The supplied data directory '%s' does not exist." % data_dir
        exit(1)

    conn = sqlite3.connect(rosetta_inputs_db)
    conn.execute("PRAGMA cache_size = 10000;")
    conn.execute("DROP TABLE IF EXISTS resfiles;")
    conn.execute("CREATE TABLE resfiles (tag TEXT, resfile TEXT, PRIMARY KEY(tag));")
    
    n_pdbs = 0
    for fname in glob.glob(data_dir+"/????FH_?.pdb"):
        n_pdbs += 1
        chain = fname.split("/")[-1].split(".")[0].split("_")[1]
        if chain == " ": chain = "_"
        tag = fname.split("/")[-1]
        resfile = "NATRO\nSTART\n* " + chain + " NATAA"
        conn.execute("INSERT INTO resfiles VALUES (?,?);", (tag, resfile))

    conn.commit()
    conn.close()

    print "Generated resfiles for %s pdb files." % n_pdbs


def res_task(res_num, res_atoms, chain):
    #does it have enough main chain atoms?
    if (" N  " in res_atoms and " CA " in res_atoms and " C  " in res_atoms) or \
            (" CA " in res_atoms and " C  " in res_atoms and " O  " in res_atoms):
        return "%s %s NATAA\n" % (res_num, chain)
    else:
        return ""


def main(argv):

    parser = OptionParser(usage="Usage: %prog [OPTIONS]")


    parser.add_option("--data_dir",
      help="A directory containing pdb files of the form ???FH_?.pdb where the final ? indicates the chain letter of the chain that should be used.")

    parser.add_option("--rosetta_inputs_db",
      default="rosetta_inputs.db3",
      help="The path to a sqlite3 database where the resfiles should be stored.")

    (options, args) = parser.parse_args(args=argv)
    
    if len(argv) == 0:
        parser.print_usage()

    setup_resfiles(options.data_dir, options.rosetta_inputs_db)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

