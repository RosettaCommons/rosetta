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
#include <core/io/pdb/pdb_writer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <core/scoring/Energies.hh>

#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>


// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace core;
using namespace core::scoring;
using namespace basic::options;
using namespace core::optimization;

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.whitney_prepare_native.main" );


/// General testing code
int
main( int argc, char * argv [] )
{

	try {


	devel::init(argc, argv);

	TR << "Starting minimization and repacking" << std::endl;

	// scoring function
	scoring::ScoreFunctionOP scorefxn( get_score_function() );
	scoring::ScoreFunctionOP repack_scorefxn( get_score_function() );
	//	scoring::ScoreFunctionOP repack_scorefxn( ScoreFunctionFactory::create_score_function(SOFT_REP_WTS) );
	//	scoring::ScoreFunctionOP repack_scorefxn( ScoreFunctionFactory::create_score_function(SOFT_REP_DESIGN_WTS) );

	repack_scorefxn->set_weight( core::scoring::fa_dun, 0.1 ); // this works for BAFF (1oqe)
	//	repack_scorefxn->set_weight( core::scoring::fa_dun, 0.2 ); // this does NOT work for BAFF (1oqe)

	// create pose for native pose from pdb
	pose::Pose native_pose;
	std::string const input_pdb_name ( basic::options::start_file() );
	core::import_pose::pose_from_file( native_pose, input_pdb_name , core::import_pose::PDB_file);

	(*scorefxn)(native_pose);
	native_pose.dump_scored_pdb( "init_pose.pdb", *scorefxn );
	TR << "Initial score: " << native_pose.energies().total_energies()[ total_score ] << std::endl;

	/*
	// rearrange fold-tree for protein-protein interface case, put jumps connecting the centers of mass
	{ // fold-tree scope

		using namespace kinematics;

		// JK NOTE: EXPECTS TWO JUMPS (IE. THREE CHAINS) WITH THE FIRST ONE THE LIGAND (IE. RECEPTORS MOVE RELATIVE TO THIS). THIS SHOULD BE GENERALIZED....

		FoldTree f( native_pose.fold_tree() );

		if ( f.num_jump() != 2 ) {
			TR << "BAFF is expected to have 2 jumps, instead found " << f.num_jump() << std::endl;
			exit(1);
		}

		FArray1D_int cuts( 2 );
		char first_chain = native_pose.pdb_info()->chain( 1 );
		char second_chain = '_';
		for ( Size i=2; i<= native_pose.total_residue(); ++i ) {
			char curr_chain = native_pose.pdb_info()->chain( i );
			if ( curr_chain != first_chain ) {
				cuts(1) = i-1;
				second_chain = curr_chain;
				break;
			}
		}
		for ( Size i=2; i<= native_pose.total_residue(); ++i ) {
			char curr_chain = native_pose.pdb_info()->chain( i );
			if ( ( curr_chain != first_chain ) && ( curr_chain != second_chain ) ) {
				cuts(2) = i-1;
				break;
			}
		}
		//		TR << "cuts are after residues " << cuts(1) << " and " << cuts(2) << std::endl;

		FArray2D_int jump_points(2,2);
		jump_points(1,1) = core::pose::residue_center_of_mass( native_pose, 1, cuts(1) );
		jump_points(2,1) = core::pose::residue_center_of_mass( native_pose, cuts(1)+1, cuts(2) );
		jump_points(1,2) = jump_points(1,1);
		jump_points(2,2) = core::pose::residue_center_of_mass( native_pose, cuts(2)+1, native_pose.total_residue() );
		//		TR << "jump1 is " << jump_points(1,1) << " and " << jump_points(2,1) << std::endl;
		//		TR << "jump2 is " << jump_points(1,2) << " and " << jump_points(2,2) << std::endl;

		bool successful_tree = f.tree_from_jumps_and_cuts( native_pose.total_residue(), 2, jump_points, cuts, 1, true );

		if ( ! successful_tree ) {
			TR << "tree_from_jumps_and_cuts was NOT successful" << std::endl;
			exit(1);
		}

		f.check_fold_tree();
		native_pose.fold_tree( f );

	} // fold-tree scope
*/

	// setting degrees of freedom which can move during minimization - sidechains only
	kinematics::MoveMap mm_sc;
	mm_sc.set_chi( true );
	mm_sc.set_bb( false );
	mm_sc.set_jump( false );

	// setting degrees of freedom which can move during minimization - everything
	kinematics::MoveMap mm_all;
	mm_all.set_chi( true );
	mm_all.set_bb( true );
	mm_all.set_jump( true );

	// minimize protein
	TR << "Starting minimization...." << std::endl;
	AtomTreeMinimizer minimizer;
	MinimizerOptions min_options( "lbfgs_armijo_nonmonotone", 0.00001, true, false );

	//	minimizer.run( native_pose, mm_sc, *scorefxn, min_options );
	//	minimizer.run( native_pose, mm_sc, *scorefxn, min_options );
	//	minimizer.run( native_pose, mm_sc, *scorefxn, min_options );
	//	minimizer.run( native_pose, mm_sc, *scorefxn, min_options );
	//	(*scorefxn)(native_pose);

	minimizer.run( native_pose, mm_all, *scorefxn, min_options );
	minimizer.run( native_pose, mm_all, *scorefxn, min_options );
	minimizer.run( native_pose, mm_all, *scorefxn, min_options );
	minimizer.run( native_pose, mm_all, *scorefxn, min_options );
	(*scorefxn)(native_pose);

	// calculate score and print it out
	(*scorefxn)(native_pose);
	TR << "Post minimization 1 score: " << native_pose.energies().total_energies()[ total_score ] << std::endl;

	// Setup packer task for repacking
  pack::task::PackerTaskOP base_packer_task( pack::task::TaskFactory::create_packer_task( native_pose ));
	base_packer_task->set_bump_check( false );
	base_packer_task->initialize_from_command_line();
	base_packer_task->or_include_current( true ); // jk absolutely critical for BAFF case, Tyr65 is unusual

	for ( Size ii = 1; ii <= native_pose.total_residue(); ++ii ) {
		base_packer_task->nonconst_residue_task(ii).restrict_to_repacking();
	}

	// First repack
	pack::pack_rotamers( native_pose, *repack_scorefxn, base_packer_task );

	// Report Scores
	(*scorefxn)(native_pose);
	native_pose.dump_scored_pdb( "repacked_once.pdb", *scorefxn );
	TR << "Score after repacking once: " << native_pose.energies().total_energies()[ total_score ] << std::endl << std::endl;

	// iterate over minimizing and repacking
	for ( Size iter = 1; iter <= 5; ++iter ) {
		minimizer.run( native_pose, mm_all, *scorefxn, min_options );
		minimizer.run( native_pose, mm_all, *scorefxn, min_options );
		minimizer.run( native_pose, mm_all, *scorefxn, min_options );
		(*scorefxn)(native_pose);
		TR << "Current score after minimizing: " << native_pose.energies().total_energies()[ total_score ] << std::endl << std::endl;
		pack::pack_rotamers( native_pose, *repack_scorefxn, base_packer_task );
		(*scorefxn)(native_pose);
		TR << "Current score after repacking: " << native_pose.energies().total_energies()[ total_score ] << std::endl << std::endl;
	}

	// final minimization
	minimizer.run( native_pose, mm_all, *scorefxn, min_options );
	minimizer.run( native_pose, mm_all, *scorefxn, min_options );
	minimizer.run( native_pose, mm_all, *scorefxn, min_options );
	minimizer.run( native_pose, mm_all, *scorefxn, min_options );
	(*scorefxn)(native_pose);
	native_pose.dump_scored_pdb( "minimized_final.pdb", *scorefxn );
	TR << "Final score: " << native_pose.energies().total_energies()[ total_score ] << std::endl << std::endl;

	TR << "Successfully finished minimizing input." << std::endl;

	return 0;


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}


