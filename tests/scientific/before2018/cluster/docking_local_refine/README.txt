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

1) Download input structures

   cd %(minidir)s/cluster/docking
   svn checkout https://svn.rosettacommons.org/source/trunk/mini.data/tests/scientific/cluster/docking/ inputs
   cd inputs


2) Compile Rosetta

   cd mini_base_dir
   ./scons.py bin mode=release -j<num_processors>

3) Run test

   cd mini_base_dir/test/scientific
   ./scientific.py docking_local_refine -d mini_database_dir

4) Investigate results

   cd mini_base_dir/test/scientific/statistics/docking_local_refine
   -> outputs

____________________________________________________

ABOUT
____________________________________________________


Docking Local Refine Scientific Benchmark

Docking is a an approach to predicting how molecules bind
together. Given the fixed chemical sequences and initial conformations
of the individual partners as input, rigid body alignment, sidechain
and backbone degrees of freedom are sampled to search for bound
conformations that have favorable over all energy.

As a scientific benchmark, accurately predicting experimentally
validated binding interactions is a stringent test of the energy
function and conformational sampling. Depending upon the level of
detail of the input data and how successful predictions are evaluated,
the docking scientific benchmark can emphasize the energy function or
the conformational sampling.

The docking local refine scientific benchmark emphasizes evaluation of
the energy function by providing the experimentally validated bound
conformation of both interaction partners as input and limiting the
allowed confomational sampling. A physically realisitic score function
should recovery the native conformation a substantial fraction of the
time. Since sampling is not required to recover the native
conformation, any significant deviation from the native conformation
indicates that configuration molecular interactions in the native
conformation is inappropriately disfavored relative to alternative
binding conformations.

For documentation about how to run the docking_protocol application see
http://graylab.jhu.edu/Rosetta.Developer.Documentation/all_else/d0/de4/docking_protocol.html



