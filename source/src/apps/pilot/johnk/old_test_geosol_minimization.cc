// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @author jk

// Project Headers
#include <devel/init.hh>
#include <core/types.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pose/Pose.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Atom.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/EtableOptions.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/hbonds/constants.hh>
#include <basic/Tracer.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>

#include <numeric/constants.hh>

//#include <core/scoring/ScoreFunction.hh>
//#include <core/scoring/ScoreFunctionFactory.hh>

#include <basic/options/option_macros.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/excn/Exceptions.hh>

// C++ Headers
#include <cmath>
#include <iostream>
#include <iomanip>
#include <map>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace core;
using namespace core::scoring;
using namespace core::scoring::hbonds;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::optimization;

static basic::Tracer TR( "apps.pilot.johnk_test_geosol_minimization.main" );

/// General testing code
int
main( int argc, char * argv [] )
{
	try {
		devel::init(argc, argv);

		TR << "jk testing derivatives for geometric solvation" << std::endl;

		// scoring function
		scoring::ScoreFunctionOP scorefxn( get_score_function() );

		// scorefxn->reset();
		// scorefxn->set_weight( core::scoring::fa_sol, 0.65 );
		scorefxn->set_weight( core::scoring::occ_sol_pw, 0.65 );

		// Build a peptide with multiple aa's, to use as a template for measuring bond lengths and angles...
		pose::Pose pose;
		std::string const input_pdb_name ( basic::options::start_file() );
		core::import_pose::pose_from_file( pose, input_pdb_name , core::import_pose::PDB_file);

		// setting degrees of freedom which can move during minimization - everything
		kinematics::MoveMap mm_all;
		mm_all.set_chi( true );
		mm_all.set_bb( true );
		mm_all.set_jump( true );

		// minimize protein with deriv check on
		TR << "Starting minimization...." << std::endl;
		AtomTreeMinimizer minimizer;
		MinimizerOptions min_options( "lbfgs_armijo_nonmonotone", 0.0001, true, true, false );

		minimizer.run( pose, mm_all, *scorefxn, min_options );
		(*scorefxn)(pose);

		TR << "jk finished testing derivatives" << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

