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

#ifndef INCLUDED_external_calibur_PreloadedPDB_HH
#define INCLUDED_external_calibur_PreloadedPDB_HH

#define LONGEST_CHAIN 4000



#include <map>
#include <string>

#include <protocols/cluster/calibur/SimPDB.fwd.hh>
#include <protocols/cluster/calibur/pdb_util.hh>



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

