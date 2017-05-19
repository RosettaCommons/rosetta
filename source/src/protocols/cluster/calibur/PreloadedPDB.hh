// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file external/calibur/PreloadedPDB.hh
/// @author SC Li & YK Ng (kalngyk@gmail.com)
/// @author Andy Watkins (amw579@stanford.edu)

#ifndef external_calibur_PreloadedPDB_HH
#define external_calibur_PreloadedPDB_HH

#define LONGEST_CHAIN 4000

#include <iostream>
#include <fstream>

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>

#include <map>
#include <vector>

#include <protocols/cluster/calibur/SimPDB.hh>
#include <protocols/cluster/calibur/pdb_util.hh>

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace cluster {
namespace calibur {

class PreloadedPDB
{
public:

	static unsigned int ADVISED_THRESHOLD;

private:
	std::string silentfilename_ = "";
	std::string pdblistfilename_ = "";

public:
	int num_residue_;
	int mNumDecoy;
	StringVecOP names_;
	std::map<std::string, SimPDBOP > filename2PDB; // fix this if it is deemed too slow


public:
	PreloadedPDB();
	~PreloadedPDB() {};
	void loadSilentFile( std::string const & silentfilename );
	void loadPDBFromList( std::string const & pdblistfilename );

	SimPDBOP getSimPDB( std::string const & pdbfilename ); // unused. for internal testing
};


}
}
}

#endif

