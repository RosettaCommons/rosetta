// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file
/// @brief

// libRosetta headers
//#include <basic/options/option.hh>

#include <core/types.hh>
#include <core/id/AtomID.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <protocols/simple_moves/MinMover.hh>

#include <devel/init.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/basic.hh>
#include <basic/database/open.hh>

#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>

#include <numeric/xyzVector.io.hh>
#include <numeric/xyzMatrix.io.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <map>
#include <cstdlib>

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>

//silly using/typedef

using namespace core;
using namespace core::id;
using namespace protocols::simple_moves;
using namespace pose;
using namespace chemical;
using namespace utility;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void
calculate_sasa()
{
	using namespace pose;
	using namespace conformation;
	using namespace chemical;
	using namespace core::import_pose;
	using namespace core::scoring;


	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	Pose pose;
	core::import_pose::pose_from_file( pose, option[ in::file::s ]().vector().front() );

	Pose start_pose = pose;

	std::string weights( "sasa_only" );
	ScoreFunctionOP score_fxn( ScoreFunctionFactory::create_score_function( weights ) );

	Real sasa( (*score_fxn)( pose ) );
	std::cout << "Solvent Accessible Surface Area is: " << sasa << std::endl;

	return;

}

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {
		//using namespace core;
		devel::init( argc, argv );

		calculate_sasa();
	} catch (utility::excn::Exception const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
