// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file external/calibur/pdb_util.hh
/// @author SC Li & YK Ng (kalngyk@gmail.com)
/// @author Andy Watkins (amw579@stanford.edu)

#ifndef external_calibur_pdb_util_HH
#define external_calibur_pdb_util_HH

#define LONGEST_CHAIN 4000

#include <iostream>
#include <fstream>

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>

#include <map>
#include <vector>

#include <protocols/cluster/calibur/SimPDB.fwd.hh>

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace cluster {
namespace calibur {

typedef std::vector< std::string > StringVec;
typedef utility::pointer::shared_ptr< StringVec > StringVecOP;

enum INPUT_FILE_TYPE { UNKNOWN=-1, SILENT_FILE, PDB_LIST };
INPUT_FILE_TYPE filetype( std::string const & filename );
unsigned int num_lines_in_file( std::string const & filename );

//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

int toInt( std::string const & aString );
double toFloat( std::string const & aString );
void center_residues( std::vector<double> & calpha_vector, int num_residue_ );

INPUT_FILE_TYPE
filetype( std::string const & filename );

unsigned int num_lines_in_file(std::string const & filename);

}
}
}

#endif
