// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @brief
/// @author jk

#include <iostream>
#include <iomanip>

// Protocol Headers
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/backrub/BackrubMover.hh>
#include <protocols/moves/MonteCarlo.hh>

// Core Headers
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/pose/PDBInfo.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <basic/options/util.hh>
#include <basic/options/after_opts.hh>
#include <core/id/AtomID_Map.hh>
#include <protoc/PocketConstraint.hh>
#include <core/scoring/func/ScalarWeightedFunc.hh>

#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <basic/options/option_macros.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/excn/Exceptions.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>



using namespace core;
using namespace basic::options;
using namespace core::scoring;
using namespace core::optimization;
using namespace basic::options::OptionKeys;

OPT_KEY( Integer, central_relax_pdb_num )
OPT_KEY( Real, pocket_constraint_weight )

/// General testing code
int
main( int argc, char * argv [] )
{
    try {
	NEW_OPT( central_relax_pdb_num, "which residue to carry out backrub around", 1 );
	NEW_OPT( pocket_constraint_weight, "weight to use for the pocket constraint", 1. );

	devel::init(argc, argv);

	std::cout << "Starting to relax backbone around central residue" << std::endl;

	core::Real const mc_kT = 0.6;
	std::cout << "Using kT = " << mc_kT << " for backrub MC...." << std::endl;

	std::string const output_tag = option[ OptionKeys::out::output_tag ]();

	pose::Pose input_pose;

	//read in pdb file from command line
	std::string const input_pdb_name ( basic::options::start_file() );
	core::import_pose::pose_from_pdb( input_pose, input_pdb_name );

	// This is the residue we'll backrub around!!
	int const central_relax_pdb_number = option[ central_relax_pdb_num ];
	core::Size central_relax_res = 0;
	for ( int j = 1, resnum = input_pose.total_residue(); j <= resnum; ++j ) {
		if ( input_pose.pdb_info()->number(j) == central_relax_pdb_number ) {
			central_relax_res = j;
		}
	}
	if ( central_relax_res == 0 ) {
		std::cerr << "ERROR!! Could not find residue to backrub around" << std::endl;
		exit(1);
	}

	scoring::ScoreFunctionOP scorefxn( getScoreFunction() );
	scorefxn->set_weight( core::scoring::pocket_constraint, option[ pocket_constraint_weight ] );

	// report starting score
	(*scorefxn)(input_pose);
  core::Real const starting_total_score = input_pose.energies().total_energies()[ total_score ];
	std::cout << "Total score at start without constraint is: " << starting_total_score << std::endl;
  core::Real const starting_pocket_score = input_pose.energies().total_energies()[ pocket_constraint ];
	std::cout << "Constraint score (unweighted) at start without constraint is: " << starting_pocket_score << std::endl;

	// apply pocket constraint
	input_pose.add_constraint( new core::scoring::constraints::PocketConstraint( input_pose, central_relax_res ) );

	// rescore, report new score
	(*scorefxn)(input_pose);
  core::Real const constraint_total_score = input_pose.energies().total_energies()[ total_score ];
	std::cout << "Total score at start with constraint is: " << constraint_total_score << std::endl;
  core::Real const constraint_pocket_score = input_pose.energies().total_energies()[ pocket_constraint ];
	std::cout << "Constraint score (unweighted) at start with constraint is: " << constraint_pocket_score << std::endl;

	// JK DEBUG - DON'T DO THE BACKRUB
	std::cerr << "JK DEBUG - JUMPOUT" << std::endl;
	exit(1);

	// create and open file for printing scores
	utility::io::ozstream ddg_outstream;
	ddg_outstream.open( "backrub."+output_tag+".out", std::ios::out );

	// create string for header for file
	std::ostringstream header_string_stream;
	header_string_stream << std::setw(9) << "Trajectory";
	header_string_stream << std::setw(15) << "total_score";
	header_string_stream << std::setw(15) << "pocket_score";
	ddg_outstream << header_string_stream.str() << std::endl;

	std::ostringstream start_data_string_stream;
	start_data_string_stream << std::setw(9) << "0";
	start_data_string_stream << std::setw(15) << constraint_total_score;
	start_data_string_stream << std::setw(15) << constraint_pocket_score;
	ddg_outstream << start_data_string_stream.str() << std::endl;

	utility::vector1 <bool>	allow_moving( input_pose.total_residue(), false );
	allow_moving.at(central_relax_res) = true;

	// find the neighbors for the central residue
	core::scoring::TenANeighborGraph const & graph = input_pose.energies().tenA_neighbor_graph();
	for ( core::graph::Graph::EdgeListConstIter
					iter = graph.get_node( central_relax_res )->const_edge_list_begin(),
					iter_end = graph.get_node( central_relax_res )->const_edge_list_end();
				iter != iter_end; ++iter ) {
		Size const neighbor_id( (*iter)->get_other_ind( central_relax_res ) );
		allow_moving.at(neighbor_id) = true;
	}


	//  Run this many trajectories
	for ( int traj = 1, ntraj = 10; traj <= ntraj; ++traj ) {

		std::cout << "Starting trajectory: " << traj << std::endl;

		// set up BackrubMover and read from the database
		protocols::backrub::BackrubMover backrubmover;
		backrubmover.branchopt().read_database();

		// Create Owning Pointer for pose
		pose::PoseOP relax_poseOP ( new pose::Pose() );
		*relax_poseOP = input_pose;

		// Setup MC
		protocols::moves::MonteCarlo mc(*relax_poseOP, *scorefxn, mc_kT );
		mc.reset(*relax_poseOP);

		//clear segments and set the input pose
		backrubmover.clear_segments();
		backrubmover.set_input_pose( relax_poseOP );

		(*scorefxn)(*relax_poseOP);

		// setup segments to move
		for ( int j = 2, resnum = input_pose.total_residue(); j < resnum; ++j ) {
			if ( ! allow_moving.at(j) ) continue;
				//			if ( j == 1 ) continue;
				//			if ( j == input_pose.total_residue() ) continue;
			if ( input_pose.pdb_info()->chain( j ) != input_pose.pdb_info()->chain( j + 1 ) ) continue;
			if ( input_pose.pdb_info()->chain( j ) != input_pose.pdb_info()->chain( j - 1 ) ) continue;
			// add current 3 residue segment to the backbone mover unless it is on another chain or isn't there
			// check to see if all 3 resiudes are on the same chain
			id::AtomID start_atom_id = id::AtomID( input_pose.residue( j - 1 ).atom_index("CA") , j - 1  );
			id::AtomID end_atom_id = id::AtomID( input_pose.residue( j + 1 ).atom_index("CA"), j + 1 );
			// report segment about to be added
			chemical::AA  start_res_neighbor( input_pose.residue( j - 1).aa() );
			chemical::AA  end_res_neighbor( input_pose.residue( j + 1).aa() );
			std::ostringstream segment1;
			segment1 << oneletter_code_from_aa( start_res_neighbor ) << input_pose.pdb_info()->number( j ) << oneletter_code_from_aa( end_res_neighbor );
			std::cout << "Current segment being added: " << segment1.str() << std::endl;
			// add current segment to the mover
			backrubmover.add_segment( start_atom_id, end_atom_id, 0 );
		}

		std::cout << "Number of segments added: " << backrubmover.num_segments() << std::endl;

		std::cout << "Score before backrub: " << relax_poseOP->energies().total_energies()[ total_score ] << std::endl;

		(*scorefxn)(*relax_poseOP);
		// optimize branch angles and idealize side chains
		backrubmover.optimize_branch_angles( *relax_poseOP );
		(*scorefxn)(*relax_poseOP);
		mc.reset(*relax_poseOP);

		std::string move_type = backrubmover.type();
		for ( int trial = 1, ntrials = 1000; trial <= ntrials; ++trial ) {
			backrubmover.apply(*relax_poseOP);
			mc.boltzmann(*relax_poseOP, move_type);
		}
		*relax_poseOP = mc.lowest_score_pose();
		(*scorefxn)(*relax_poseOP);

		// calculate score for output
		core::Real const final_total_score = relax_poseOP->energies().total_energies()[ total_score ];
		std::cout << "Total score after backrub: " << final_total_score << std::endl;
		core::Real const final_pocket_score = relax_poseOP->energies().total_energies()[ pocket_constraint ];
		std::cout << "Constraint score (unweighted) after backrub: " << final_pocket_score << std::endl;

		// Reporting Data
		std::ostringstream data_string_stream;
		data_string_stream << std::setw(9) << output_tag << "_" << traj;
		data_string_stream << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::setw(15) << final_total_score;
		data_string_stream << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::setw(15) << final_pocket_score;

		// Print to file
		ddg_outstream << data_string_stream.str() << std::endl;

		// writes out the current structure
		std::ostringstream outPDB_name;
		outPDB_name << "traj." << output_tag << "_" << traj << ".pdb";
		relax_poseOP->dump_scored_pdb( outPDB_name.str(), *scorefxn );

	}

	// close file
	ddg_outstream.close();
	ddg_outstream.clear();

	std::cout << "Successfully finished relaxing backbone around central residue" << std::endl;

    } catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
    }
    return 0;
}



