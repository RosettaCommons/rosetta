// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file external/calibur/SimPDB.hh
/// @author SC Li & YK Ng (kalngyk@gmail.com)
/// @author Andy Watkins (amw579@stanford.edu)

#ifndef external_calibur_SimPDB_HH
#define external_calibur_SimPDB_HH

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
#include <protocols/cluster/calibur/PreloadedPDB.fwd.hh>
#include <protocols/cluster/calibur/pdb_util.hh>

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace cluster {
namespace calibur {


class SimPDB
{
public:
	/**
	* These parameters control how PDB files are to be loaded.
	* They do not apply to silent file, which are assumed to be
	* pre-processed.
	*/
	static int s_residue;
	static int e_residue;
	static char * chains;

	/**
	* This feature allows the preloading of SimPDB objects.
	* SimPDB will then be obtained from preloadedPDB instead of from disk.
	*/
	static PreloadedPDB * preloadedPDB;
	static bool preloadPDB;

public:
	std::string protein_file_name_;
	int num_residue_;
	//double mSquaredSum;
	std::vector<double> calpha_vector_;
	void read();

public:
	SimPDB( std::string const & aProteinFileName );
	SimPDB( std::string const & aProteinFileName, int const len );
	~SimPDB();

	// Special constructors used only by PreloadedPDB. Don't touch.
	SimPDB();
	SimPDB(int len);
};


}
}
}

#endif
