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

#include <iostream>
#include <iomanip>

#include <devel/init.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <basic/options/util.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <basic/options/option_macros.hh>

#include <core/scoring/Energies.hh>
#include <protocols/moves/RigidBodyMover.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
#include <numeric/xyzVector.hh>
#include <basic/Tracer.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace core;
using namespace core::scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;

/*
OPT_KEY( Integer, phosphotyr_num )
OPT_KEY( String, phosphotyr_chain )
OPT_KEY( Real, match_distance_cutoff )
OPT_KEY( Real, phosphate_force_constant )
*/
OPT_KEY( Boolean, do_minimization )

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.phosphoArg_minimize.main" );


int
main( int argc, char * argv [] )
{

	/*
	NEW_OPT( phosphotyr_num, "which residue the starting pY is", 0 );
	NEW_OPT( phosphotyr_chain, "which chain is the starting pY is on", "P" );
	NEW_OPT( match_distance_cutoff, "required distance from Arg N to pTyr O", 1.0 );
	NEW_OPT( phosphate_force_constant, "force constant for the coordinate constraint used in minimization", 25.0 );
	*/
	NEW_OPT( do_minimization, "whether or not to do minimzation step", true );

	devel::init(argc, argv);

	TR << "Starting phospho-Arg calculations" << std::endl;

	TR << "Reading pose" << std::endl;

	pose::Pose input_pose;
	std::string const input_pdb_name( basic::options::start_file() );
	core::import_pose::pose_from_file( input_pose, input_pdb_name , core::import_pose::PDB_file);

	Size const totres = input_pose.size();
	pose::Pose pose = input_pose;

	// Setup for scoring/repacking
	scoring::ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function(PRE_TALARIS_2013_STANDARD_WTS, SCORE12_PATCH) );
	scorefxn->set_weight( core::scoring::fa_dun, 0. );
	scorefxn->set_weight( core::scoring::fa_intra_rep, 0. );
	(*scorefxn)(pose);

	if ( option[ do_minimization ] ) {

		// carry out minimzation of (only) sidechains
		kinematics::MoveMap mm_sc;
		mm_sc.set_chi( true );
		mm_sc.set_bb( false );
		mm_sc.set_jump( false );

		TR << "Running first minimization" << std::endl;

		core::optimization::AtomTreeMinimizer minimizer;
		core::optimization::MinimizerOptions min_options( "lbfgs_armijo_nonmonotone", 0.00001, true, false );
		minimizer.run( pose, mm_sc, *scorefxn, min_options );

		// JK DEBUG
		pose.dump_scored_pdb( basic::options::start_file()+".min1.pdb", *scorefxn );

		TR << "Running second minimization" << std::endl;
		kinematics::MoveMap mm_all;
		mm_all.set_chi( true );
		mm_all.set_bb( true );
		mm_all.set_jump( true );
		minimizer.run( pose, mm_all, *scorefxn, min_options );

		// JK DEBUG
		TR << "Printing minimized structure" << std::endl;
		pose.dump_scored_pdb( basic::options::start_file()+".min2.pdb", *scorefxn );

	}

	core::Real const bound_score = pose.energies().total_energies()[ total_score ];
	(*scorefxn)(pose);

	// rb_jump = rigid_backbone - changes the distance between the backbones of the molecules
	core::Real const unbound_dist = 40.; // create unbound distance
	Size const rb_jump = 1;
	protocols::moves::RigidBodyTransMover trans_mover( pose, rb_jump );
	trans_mover.trans_axis( trans_mover.trans_axis() );
	trans_mover.step_size(unbound_dist);
	trans_mover.apply( pose );

	(*scorefxn)(pose);
	core::Real const unbound_score = pose.energies().total_energies()[ total_score ];
	//	pose.dump_scored_pdb( basic::options::start_file()+".unbound.pdb", *scorefxn );

	TR << "Completed, bound - free score diff is " << bound_score - unbound_score << std::endl;

	return 0;
}


