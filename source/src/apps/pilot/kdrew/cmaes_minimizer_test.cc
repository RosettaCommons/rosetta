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
//#include <basic/options/option.hh>

// Utility headers
#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <core/id/DOF_ID.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/optimization/MinimizerMap.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/DOF_Node.hh>

#include <protocols/minimization_packing/MinMover.hh>


/*
#include <core/conformation/symmetry/SymmData.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>
*/


#include <stdio.h>
#include <stdlib.h> /* free() */
#include <stddef.h> /* NULL */
#include <algorithm>

using namespace core;
using namespace kinematics;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace pose;

//using namespace core::conformation::symmetry;

static basic::Tracer TR( "CMAES_MINIMIZER_TEST" );


int main(int argc, char *argv[]) {

	try {

		devel::init(argc, argv);
		core::pose::Pose pose;
		core::import_pose::pose_from_file(pose,basic::options::option[basic::options::OptionKeys::in::file::s]()[1]);

		core::scoring::ScoreFunctionOP score_fxn_ = scoring::get_score_function();
		TR << (*score_fxn_)(pose) << std::endl;

		//kdrew: change from "new" to "utility::pointer::make_shared"
		//core::kinematics::MoveMapOP move_map( new core::kinematics::MoveMap() );
		core::kinematics::MoveMapOP move_map( utility::pointer::make_shared< core::kinematics::MoveMap >() );

		move_map->set_bb( true );
		move_map->set_chi( true );
		move_map->set_jump( 1, true );

		move_map->show();

		core::optimization::MinimizerMap min_map;
		min_map.setup( pose, *move_map );

		TR << "pose: " << pose<< std::endl;

		Size dof_count = 0;
		TR <<  "dof"<< "  " << "rsd"<< "   "<< "AtomNo" << " " << "Type" << std::endl;
		for ( std::list< core::optimization::DOF_NodeOP >::const_iterator it = min_map.begin(), ie = min_map.end(); it != ie; ++it ) {
			//std::cout <<  "dof " << dof_count << "  " << (**it).rsd() << "   "<< (**it).atom_id() << std::endl;
			TR <<  dof_count << "  " << (**it).rsd() << "   "<< (**it).atomno() << " " << (**it).type() << std::endl;
			dof_count++;
		}

		//kdrew: extra options may not mean anything
		protocols::minimization_packing::MinMoverOP minmover( utility::pointer::make_shared< protocols::minimization_packing::MinMover >( move_map, score_fxn_, "cmaes", 0.001, true ) );
		minmover->apply(pose);

		std::stringstream minpdbname;
		minpdbname << "min_pose.pdb";
		pose.dump_scored_pdb( minpdbname.str(), *score_fxn_ );


	} catch ( utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;

}
