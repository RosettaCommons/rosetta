#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

import os, sys, glob, re
import subprocess

sys.path.insert(0, "/".join(os.getcwd().split("/")[:-1]))

from submit import BaseSampleSource

class SampleSource(BaseSampleSource):


    def __init__(self, argv):
        BaseSampleSource.__init__(self)

        self.sample_source_description_default_value = \
            "Top 4400 chains in PDB from the Richardson Lab with the O'Meara Leaver-Fay corrections and Reduce placed hydrogens"

        self.initialize_options_parser()
        self.parse_options(argv)
        self.setup()

#    def setup_input_data(self):
#
#        print "Checkout %(sample_source_id)s dataset..." % self.mvars
#        try:
#	    p = subprocess.Popen(
#	        args=['svn', 'checkout', "--non-interactive", "--no-auth-cache",
#	              'https://svn.rosettacommons.org/source/trunk/mini.data/tests/scientific/cluster/features/sample_sources/top4400/input',
#	              'input'],
#	        stderr=subprocess.PIPE)
#	    err = p.communicate()[1]
#	    if err != '': print err
#        except OSError, e:
#            print >>sys.stderr, "Execution Failed:", e
#            
#
#        print "unzipping top4400pdbs.tar.gz..."
#        try:
#            p = subprocess.Popen(['tar', '-xzf', 'top4400pdbs.tar.gz'],
#                  cwd='input', stderr=subprocess.PIPE)
#            err = p.communicate()[1]
#            if err != '': print err
#        except OSError, e:
#            print >>sys.stderr, "Execution Failed:", e
#
#        print "restricting to the relevant chain in each pdb..."
#        for pdb_fname in glob.glob("input/top4400pdbs/*pdb"):
#            chain = re.match(r".*FH([A-Z,0-9])\.pdb", pdb_fname).group(1)
#            f = open(pdb_fname, "r")
#            relevant_atoms = re.findall(r"ATOM.{17}%s.*\n" % chain, f.read())
#            f.close()
#            f = open(pdb_fname, "w")
#            f.write("".join(relevant_atoms))
#            f.close()
#   
#        print "generating all_pdbs.list..."
#        try:
#            f = open("input/all_pdbs.list", 'w')
#            p = subprocess.Popen(["find", ".", "-name", "*pdb"],
#                  stdout=f, stderr=subprocess.PIPE, cwd=os.getcwd() + "/input/top4400pdbs")
#                 stdout=f, stderr=subprocess.PIPE)
#            err = p.communicate()[1]
#            if err != '': print err
#            f.close()
#        except OSError, e:
#            print >>sys.stderr, "Execution Failed:", e
#
#        print "Done setting up input files"


def main(argv):

    ss = SampleSource(argv)
    ss.submit()


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))


