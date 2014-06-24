#!/usr/bin/env bash
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

find ./src/{apps,core,devel,numeric,protocols,utility} ./test -name '*.cc' -o -name '*.hh' -print > tagged_files.txt
# The -L flag requires Exuberant Ctags
ctags -L tagged_files.txt
cscope -bqk -i tagged_files.txt
