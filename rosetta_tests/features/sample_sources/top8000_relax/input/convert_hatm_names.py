
####
#
# Run this script
#
#    python setup_resfiles.py --data_dir <path_to_pdb_files>  --output_dir
#
# For all pdbs in the data_dir, convert the names of the hydrogens
# from the PDB3.2 format to Rosetta's format
#
###

#this messes up prolines at the N-termini, so handle those with the
#other proline atoms below
termini = {
  " H1 ": "1H  ",  
  " H2 ": "2H  ", 
  " H3 ": "3H  ", 
}


rename = {
# PDB3.2       ->  rosetta  
  ("ALA", " HB1"): "1HB ", 
  ("ALA", " HB2"): "2HB ", 
  ("ALA", " HB3"): "3HB ", 
  ("ARG", " HB2"): "1HB ", 
  ("ARG", " HB3"): "2HB ", 
  ("ARG", " HG2"): "1HG ", 
  ("ARG", " HG3"): "2HG ", 
  ("ARG", " HD2"): "1HD ", 
  ("ARG", " HD3"): "2HD ", 
  ("ARG", "HH11"): "1HH1", 
  ("ARG", "HH12"): "2HH1", 
  ("ARG", "HH21"): "1HH2", 
  ("ARG", "HH22"): "2HH2", 
  ("ASN", " HB2"): "1HB ", 
  ("ASN", " HB3"): "2HB ", 
  ("ASN", "HD21"): "1HD2", 
  ("ASN", "HD22"): "2HD2", 
  ("ASP", " HB2"): "1HB ", 
  ("ASP", " HB3"): "2HB ", 
  ("CYS", " HB2"): "1HB ", 
  ("CYS", " HB3"): "2HB ", 
  ("GLN", " HB2"): "1HB ", 
  ("GLN", " HB3"): "2HB ", 
  ("GLU", " HB2"): "1HB ", 
  ("GLU", " HB3"): "2HB ", 
  ("GLU", " HG2"): "1HG ", 
  ("GLU", " HG3"): "2HG ", 
  ("GLN", " HG2"): "1HG ", 
  ("GLN", " HG2"): "1HG ", 
  ("GLN", " HG3"): "2HG ", 
  ("GLN", " HG3"): "2HG ", 
  ("GLN", " HG2"): "1HG ", 
  ("GLN", " HG3"): "2HG ", 
  ("GLN", "HE21"): "1HE2",
  ("GLN", "HE22"): "2HE2",
  ("GLY", " HA2"): "1HA ", 
  ("GLY", " HA3"): "2HA ",
  ("HIS", " HB2"): "1HB ", 
  ("HIS", " HB3"): "2HB ", 
  ("ILE", "HG12"): "1HG1", 
  ("ILE", "HG13"): "2HG1", 
  ("ILE", "HG21"): "1HG2", 
  ("ILE", "HG22"): "2HG2", 
  ("ILE", "HG23"): "3HG2", 
  ("ILE", "HD11"): "1HD1", 
  ("ILE", "HD12"): "2HD1", 
  ("ILE", "HD13"): "3HD1",
  ("LEU", " HB2"): "1HB ", 
  ("LEU", " HB3"): "2HB ",
  ("LEU", "HD11"): "1HD1", 
  ("LEU", "HD12"): "2HD1", 
  ("LEU", "HD13"): "3HD1",  
  ("LEU", "HD21"): "1HD2", 
  ("LEU", "HD22"): "2HD2", 
  ("LEU", "HD23"): "3HD2", 
  ("LYS", " HB2"): "1HB ", 
  ("LYS", " HB3"): "2HB ", 
  ("LYS", " HG2"): "1HG ", 
  ("LYS", " HG3"): "2HG ", 
  ("LYS", " HD2"): "1HD ", 
  ("LYS", " HD3"): "2HD ", 
  ("LYS", " HE2"): "1HE ", 
  ("LYS", " HE3"): "2HE ", 
  ("LYS", " HZ1"): "1HZ ", 
  ("LYS", " HZ2"): "2HZ ", 
  ("LYS", " HZ3"): "3HZ ", 
  ("MET", " HB2"): "1HB ", 
  ("MET", " HB3"): "2HB ", 
  ("MET", " HG2"): "1HG ", 
  ("MET", " HG3"): "2HG ", 
  ("MET", " HE1"): "1HE ", 
  ("MET", " HE2"): "2HE ", 
  ("MET", " HE3"): "3HE ",   
  ("PHE", " HB2"): "1HB ", 
  ("PHE", " HB3"): "2HB ", 
  ("PRO", " H2 "): "1H  ", # this is a special case for termini
  ("PRO", " H3 "): "2H  ", # this is a special case for termini
  ("PRO", " HB2"): "1HB ", 
  ("PRO", " HB3"): "2HB ", 
  ("PRO", " HG2"): "1HG ", 
  ("PRO", " HG3"): "2HG ", 
  ("PRO", " HD2"): "1HD ", 
  ("PRO", " HD3"): "2HD ", 
  ("SER", " HB2"): "1HB ", 
  ("SER", " HB3"): "2HB ", 
  ("THR", "HG21"): "1HG2", 
  ("THR", "HG22"): "2HG2", 
  ("THR", "HG23"): "3HG2", 
  ("TRP", " HB2"): "1HB ", 
  ("TRP", " HB3"): "2HB ", 
  ("TYR", " HB2"): "1HB ", 
  ("TYR", " HB3"): "2HB ", 
  ("TYR", "HG21"): "1HG2", 
  ("TYR", "HG22"): "2HG2", 
  ("TYR", "HG23"): "3HG2", 
  ("VAL", "HG11"): "1HG1", 
  ("VAL", "HG12"): "2HG1", 
  ("VAL", "HG13"): "3HG1", 
  ("VAL", "HG21"): "1HG2", 
  ("VAL", "HG22"): "2HG2", 
  ("VAL", "HG23"): "3HG2", 
}
  


import glob, os, sqlite3, sys
from optparse import OptionParser

def run(data_dir, output_dir):

    if not os.path.exists(data_dir): 
        print "ERROR: The supplied data directory '%s' does not exist." % data_dir
        exit(1)

    n_pdbs = 0
    for fname in glob.glob(data_dir+"/*.pdb"):
        n_pdbs += 1
        f = open(fname)
        fout = open(output_dir +"/" + fname.split("/")[-1], 'w')
        for line in f.read().split("\n"):
            res_name = line[17:20]
            atom_name = line[12:16]
            rosetta_name = atom_name # use this by default
            if (res_name, atom_name) in rename:
                rosetta_name = rename[(res_name, atom_name)]
            else:
                if atom_name in termini:
                    rosetta_name = termini[atom_name]

            
            line = line[:12] + rosetta_name + line[16:]
            fout.write("%s\n" % line)

def main(argv):

    parser = OptionParser(usage="Usage: %prog [OPTIONS]")


    parser.add_option("--data_dir")
    parser.add_option("--output_dir")

    (options, args) = parser.parse_args(args=argv)
    
    if len(argv) == 0:
        parser.print_usage()

    run(options.data_dir, options.output_dir)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

