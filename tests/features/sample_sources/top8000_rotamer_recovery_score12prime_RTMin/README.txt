# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file test/scientific/tests/rotamer_recovery/README.txt


____________________________________________________

HOW TO RUN
____________________________________________________

0) Check Prerequisites

   a) Rosetta
   b) Input Structures ( see next )
   c) R  version >= 1.11
      c.1) with packages 'plyr', 'ggplot2', 'reshape', 'RSQLite'  (see ?install.packages in R)

1) Download input structures

   cd %(minidir)s/tests/rotamer_recovery
   svn checkout https://svn.rosettacommons.org/source/trunk/mini.data/tests/scientific/tests/rotamer_recovery/inputs inputs


2) Compile Rosetta

   cd mini_base_dir
   ./scons.py bin mode=release -j<num_processors>

3) Run test

   cd mini_base_dir/test/scientific
   ./scientific.py rotamer_recovery -d mini_database_dir

4) Investigate results

   cd mini_base_dir/test/scientific/statistics/rotamer_recovery
   -> outputs

____________________________________________________

ABOUT
____________________________________________________


 Rotamer Recovery Scientific Benchmark

The rotamer recovery scientific benchmark addresses the question,
"Given the sequence identity and backbone conformation of a
experimentally characterized protein structure, how
accurately can Rosetta predict the conformation of the sidechains?"
Having an accurate answer will:

 a) Inform researchers on the trustworthiness of Rosetta as molecular modeler
 b) Guide tuning Rosetta to improve structural predictions

The input structures is a selection of the Richardson's Top5200
dataset.  The Top5200 set was constructed as follows.  Daniel Keedy
and others the Richardson Lab clustered all structures in the PDB on
April 5th 2007 into 70% sequence homology clusters.  Each structure
with resolution 2.2A or better and was not filtered by hand was run
through MolProbity.  Then from each group the structure with the best
average resolution and MolProbity score was selected provided it had
resolution at least 2.0A.  All structures from the Top5200 having
between 50 and 200 residues and resolution less than 1.2A were
selected for this benchmark. This leaves 152 structures with 17463
residues.

