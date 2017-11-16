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


#include <core/chemical/ChemicalManager.hh>
//#include <protocols/moves/ResidueMover.hh>

#include <protocols/viewer/viewers.hh>
#include <core/types.hh>

#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>

#include <core/scoring/etable/Etable.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/pose/Pose.hh>

#include <basic/options/option.hh>
#include <basic/options/util.hh>

#include <devel/init.hh>

#include <core/io/pdb/pdb_writer.hh>

// Utility Headers
#include <utility/vector1.hh>

// Numeric Headers
#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <ctime>

//silly using/typedef


#include <core/import_pose/import_pose.hh>
#include <basic/Tracer.hh>
using basic::Error;
using basic::Warning;
using namespace core;

using utility::vector1;



typedef std::map< std::string, std::map< std::string, numeric::xyzVector< Real > > > Coords;

typedef vector1< std::string > Strings;


///////////////////////////////////////////////////////////////////////////////
void
score_pdb( std::string pdb )
{
	using namespace core;
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace pose;

	ScoreFunctionOP scorefxn( get_score_function_legacy( PRE_TALARIS_2013_STANDARD_WTS ) );

	Pose pose;
	core::import_pose::pose_from_file( pose, pdb , core::import_pose::PDB_file);

	/*Energy score_orig = */ (*scorefxn)( pose );
	EnergyMap emap_orig = pose.energies().total_energies();
	// emap_orig *= scorefxn->weights();

	std::cout << emap_orig << std::endl;

}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {

		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		devel::init(argc, argv);

		utility::vector1< std::string > pdbnames( basic::options::start_files() );
		for ( Size ii = 1; ii <= pdbnames.size(); ++ii ) {
			score_pdb( pdbnames[ ii ]);
		}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

