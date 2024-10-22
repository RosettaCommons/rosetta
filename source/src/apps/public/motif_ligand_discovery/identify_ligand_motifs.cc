// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief

// libRosetta headers
//#include <core/options/option.hh>

#include <sstream>
#include <string>

#include <devel/init.hh>
#include <utility/excn/Exceptions.hh>

#include <protocols/motifs/IdentifyLigandMotifs.fwd.hh>
#include <protocols/motifs/IdentifyLigandMotifs.hh>

#include <protocols/motifs/MotifLibrary.hh>
#include <protocols/motifs/motif_utils.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <utility/io/ozstream.hh>
#include <utility/file/FileName.hh>

// Time profiling header
#include <time.h>

using namespace core;
//using namespace ObjexxFCL;
using namespace pose;
//using namespace chemical;
using namespace scoring;
//using namespace optimization;



int
main( int argc, char * argv [] )
{
	try {
		devel::init( argc, argv );



		//create and use IdentifyLigandMotifs
		IdentifyLigandMotifs ilm;
		//create MotifLibrary to hold motif library that ILM will generate

		ilm.process_file_list();
		ilm.write_motifs_to_disk();
	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}
