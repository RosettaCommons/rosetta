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

// Protocol Headers
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/backrub/BackrubMover.hh>
#include <protocols/moves/MonteCarlo.hh>

// Core Headers
#include <core/conformation/Residue.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
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
#include <basic/Tracer.hh>
#include <core/scoring/Energies.hh>

#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <basic/options/option_macros.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace core;
using namespace basic::options;
using namespace core::scoring;
using namespace core::optimization;
using namespace basic::options::OptionKeys;

OPT_KEY( Integer, central_relax_pdb_num )
OPT_KEY( String, central_relax_chain )

static basic::Tracer TR( "apps.pilot.whitney_local_backrub.main" );

/// General testing code
int
main( int argc, char * argv [] )
{

	try {


		NEW_OPT( central_relax_pdb_num, "which residue to carry out backrub around", 1 );
		NEW_OPT( central_relax_chain, "which chain the residue is on", "A" );

		devel::init(argc, argv);

		std::string const tmp_chain = option[ central_relax_chain ];
		if ( tmp_chain.length() != 1 ) {
			std::cerr << "ERROR!! Chain ID should be one character" << std::endl;
			exit(1);
		}
		char const central_relax_pdb_chain = tmp_chain[0];
		int const central_relax_pdb_number = option[ central_relax_pdb_num ];

		TR << "Starting to relax backbone around central residue" << std::endl;

		core::Real const mc_kT = 0.6;
		TR << "Using kT = " << mc_kT << " for backrub MC...." << std::endl;

		std::string const output_tag = option[ OptionKeys::out::output_tag ]();

		// create pose for wild type native bound and unbound
		pose::Pose input_pose, unbound_input_pose;

		//read in pdb file from command line
		std::string const input_pdb_name ( basic::options::start_file() );
		core::import_pose::pose_from_file( input_pose, input_pdb_name , core::import_pose::PDB_file);

		// This is the residue we'll backrub around!!
		core::Size central_relax_res = 0;
		for ( int j = 1, resnum = input_pose.size(); j <= resnum; ++j ) {
			if ( ( input_pose.pdb_info()->chain(j) == central_relax_pdb_chain ) &&
					( input_pose.pdb_info()->number(j) == central_relax_pdb_number ) ) {
				central_relax_res = j;
			}
		}
		if ( central_relax_res == 0 ) {
			std::cerr << "ERROR!! Could not find residue/chain" << std::endl;
			exit(1);
		}

		scoring::ScoreFunctionOP scorefxn( get_score_function() );
		scorefxn->set_weight( core::scoring::fa_dun, 0.1 );
		(*scorefxn)(input_pose);

		utility::vector1 <bool> allow_moving( input_pose.size(), false );
		allow_moving.at(central_relax_res) = true;

		// find the neighbors for the central residue
		core::scoring::TenANeighborGraph const & graph = input_pose.energies().tenA_neighbor_graph();
		for ( utility::graph::Graph::EdgeListConstIter
				iter = graph.get_node( central_relax_res )->const_edge_list_begin(),
				iter_end = graph.get_node( central_relax_res )->const_edge_list_end();
				iter != iter_end; ++iter ) {
			Size const neighbor_id( (*iter)->get_other_ind( central_relax_res ) );
			allow_moving.at(neighbor_id) = true;
		}

		// Setup pose for wild type unbound complex
		// create an unbound pose
		unbound_input_pose = input_pose;

		core::Real const unbound_dist = 40.; // create unbound distance

		// rb_jump = rigid_backbone - changes the distance between the backbones of the molecules
		Size const rb_jump = 1;
		protocols::rigid::RigidBodyTransMover trans_mover( unbound_input_pose, rb_jump );
		trans_mover.trans_axis( trans_mover.trans_axis() );
		trans_mover.step_size(unbound_dist);
		trans_mover.apply( unbound_input_pose );

		// repacking unbound wt
		pack::task::PackerTaskOP base_packer_task_unbound( pack::task::TaskFactory::create_packer_task( unbound_input_pose) );
		base_packer_task_unbound->set_bump_check( false );
		base_packer_task_unbound->initialize_from_command_line();
		base_packer_task_unbound->or_include_current( true );
		for ( int ii = 1, nres = unbound_input_pose.size(); ii < nres; ++ii ) {
			base_packer_task_unbound->nonconst_residue_task( ii ).restrict_to_repacking();
		}
		base_packer_task_unbound->restrict_to_residues( allow_moving );
		pack::pack_rotamers( unbound_input_pose, *scorefxn, base_packer_task_unbound );

		//write out pdb file for unbound native
		unbound_input_pose.dump_scored_pdb( "unbound_native.pdb", *scorefxn );

		// calculate score for wild type native bound and unbound
		(*scorefxn)(input_pose);
		(*scorefxn)(unbound_input_pose);

		// calculating wt bound and unbound for native
		core::Real const starting_bound_score = input_pose.energies().total_energies()[ total_score ];
		core::Real const starting_unbound_score = unbound_input_pose.energies().total_energies()[ total_score ];

		TR << "Starting bound score for native is: " << starting_bound_score << std::endl;
		TR << "Starting unbound score for native is: " << starting_unbound_score << std::endl;

		// create delta G for starting structure
		core::Real const dg_start = starting_bound_score - starting_unbound_score;

		TR << "Delta G score for start: " << dg_start << std::endl;

		// create and open file for printing scores
		utility::io::ozstream ddg_outstream;
		ddg_outstream.open( "backrub."+output_tag+".out", std::ios::out );

		// create string for header for file
		std::ostringstream header_string_stream;
		header_string_stream << std::setw(9) << "Trajectory";
		header_string_stream << std::setw(15) << "total_score";
		header_string_stream << std::setw(15) << "delta_score";
		ddg_outstream << header_string_stream.str() << std::endl;

		std::ostringstream start_data_string_stream;
		start_data_string_stream << std::setw(9) << "0";
		start_data_string_stream << std::setw(15) << starting_bound_score;
		start_data_string_stream << std::setw(15) << dg_start;
		ddg_outstream << start_data_string_stream.str() << std::endl;


		//  Run this many trajectories
		int const ntraj = option[ OptionKeys::out::nstruct ]();
		for ( int traj = 1; traj <= ntraj; ++traj ) {

			TR << "Starting trajectory: " << traj << std::endl;

			// set up BackrubMover and read from the database
			protocols::backrub::BackrubMover backrubmover;
			backrubmover.branchopt().read_database();

			// setup minimizer
			//  AtomTreeMinimizer minimizer;
			//  MinimizerOptions min_options( "lbfgs_armijo_nonmonotone", 0.00001, true, false );
			//  kinematics::MoveMap min_mm;

			// setting degrees of freedom which can move during minimization (backbone and sidechains, not jump)
			//  min_mm.set_jump( false );
			//  for ( int j = 2, resnum = input_pose.size(); j < resnum; ++j ) {
			//   if ( ! allow_moving.at(j) ) continue;
			//   if ( input_pose.pdb_info()->chain( j ) != input_pose.pdb_info()->chain( j + 1 ) ) continue;
			//   if ( input_pose.pdb_info()->chain( j ) != input_pose.pdb_info()->chain( j - 1 ) ) continue;
			//   min_mm.set_bb( j, true );
			//   min_mm.set_chi( j, true );
			//  }

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
			for ( int j = 2, resnum = input_pose.size(); j < resnum; ++j ) {
				if ( ! allow_moving.at(j) ) continue;
				//   if ( j == 1 ) continue;
				//   if ( j == input_pose.size() ) continue;
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
				segment1 << oneletter_code_from_aa( start_res_neighbor ) <<
					input_pose.pdb_info()->number( j ) << oneletter_code_from_aa( end_res_neighbor );
				TR << "Current segment being added: " << segment1.str() << std::endl;
				// add current segment to the mover
				backrubmover.add_segment( start_atom_id, end_atom_id, 0 );
			}

			TR << "Number of segments added: " << backrubmover.num_segments() << std::endl;

			TR << "Score before backrub: " << relax_poseOP->energies().total_energies()[ total_score ] << std::endl;

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
			TR << "Score after backrub: " << relax_poseOP->energies().total_energies()[ total_score ] << std::endl;

			// final repack
			pack::task::PackerTaskOP final_repack_task( pack::task::TaskFactory::create_packer_task( *relax_poseOP ) );
			final_repack_task->set_bump_check( false );
			final_repack_task->initialize_from_command_line();
			final_repack_task->or_include_current( true );
			for ( int ii = 1, nres = input_pose.size(); ii < nres; ++ii ) {
				final_repack_task->nonconst_residue_task( ii ).restrict_to_repacking();
			}
			final_repack_task->restrict_to_residues( allow_moving );
			pack::pack_rotamers( *relax_poseOP, *scorefxn, final_repack_task );
			(*scorefxn)(*relax_poseOP);
			TR << "Score after repack: " << relax_poseOP->energies().total_energies()[ total_score ] << std::endl;

			// Gradient minimization
			//  TR << "Starting minimization...." << std::endl;
			//  minimizer.run( *relax_poseOP, min_mm, *scorefxn, min_options );
			//  (*scorefxn)(*relax_poseOP);
			//  TR << "Score after minimize: " << relax_poseOP->energies().total_energies()[ total_score ] << std::endl;

			// Calculating delta score
			core::Real const bound_optimized_bb_score = relax_poseOP->energies().total_energies()[ total_score ];
			core::Real const delta_score = bound_optimized_bb_score - starting_unbound_score;

			// Printing Results to screen
			TR << "Bound score after backbone optimization " << bound_optimized_bb_score << std::endl;
			TR << "Delta score: " << delta_score << std::endl;

			// Reporting Data
			std::ostringstream data_string_stream;
			data_string_stream << std::setw(9) << output_tag << "_" << traj;
			data_string_stream << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::setw(15) << bound_optimized_bb_score;
			data_string_stream << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::setw(15) << delta_score;

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

		TR << "Successfully finished relaxing backbone around central residue" << std::endl;

		return 0;

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}


