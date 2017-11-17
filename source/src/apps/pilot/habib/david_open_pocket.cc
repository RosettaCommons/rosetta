// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @author jk+dj

#include <iostream>
#include <iomanip>

// Protocol Headers
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/backrub/BackrubMover.hh>
#include <protocols/backrub/BackrubMover.fwd.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/moves/TrialMover.hh>


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
#include <protocols/pockets/PocketConstraint.hh>
#include <core/scoring/func/ScalarWeightedFunc.hh>
#include <core/scoring/Energies.hh>

#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <basic/options/option_macros.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/pocket_grid.OptionKeys.gen.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <utility/excn/Exceptions.hh>


using namespace core;
using namespace basic::options;
using namespace core::scoring;
using namespace core::optimization;
using namespace basic::options::OptionKeys;

OPT_KEY( Real, pocket_kT )
OPT_KEY( Boolean, pocket_SA)

/// General testing code
int
main( int argc, char * argv [] )
{
	try {

		using namespace protocols::moves;

		NEW_OPT( pocket_kT, "kT to be used for backrub MC", 1.2 );
		NEW_OPT( pocket_SA, "Use simulated annealing during backrub", false );

		devel::init(argc, argv);

		std::cout << "Starting to relax backbone around central residue" << std::endl;

		bool const useSA=option [ pocket_SA ];

		core::Real const mc_kT = option [ pocket_kT ];
		std::cout << "Using kT = " << mc_kT << " for backrub MC...." << std::endl;

		std::string const output_tag = option[ OptionKeys::out::output_tag ]();

		pose::Pose input_pose;

		//read in pdb file from command line
		std::string const input_pdb_name ( basic::options::start_file() );
		core::import_pose::pose_from_file( input_pose, input_pdb_name , core::import_pose::PDB_file);

		scoring::ScoreFunctionOP scorefxn( get_score_function() );
		//scorefxn->set_weight( core::scoring::pocket_constraint, option[ pocket_constraint_weight ] );

		// report starting score
		(*scorefxn)(input_pose);
		core::Real const starting_total_score = input_pose.energies().total_energies()[ total_score ];
		std::cout << "Total score at start without constraint is: " << starting_total_score << std::endl;
		core::Real const starting_pocket_score = input_pose.energies().total_energies()[ pocket_constraint ];
		std::cout << "Constraint score (unweighted) at start without constraint is: " << starting_pocket_score << std::endl;

		// apply pocket constraint
		protocols::pockets::PocketConstraintOP pcons( new protocols::pockets::PocketConstraint(input_pose) );
		//input_pose.add_constraint( pcons );

		// rescore, report new score
		(*scorefxn)(input_pose);
		core::Real const constraint_total_score = input_pose.energies().total_energies()[ total_score ];
		std::cout << "Total score at start with constraint is: " << constraint_total_score << std::endl;
		core::Real const constraint_pocket_score = input_pose.energies().total_energies()[ pocket_constraint ];
		std::cout << "Constraint score (unweighted) at start with constraint is: " << constraint_pocket_score << std::endl;

		// JK DEBUG - DON'T DO THE BACKRUB
		//std::cerr << "JK DEBUG - JUMPOUT" << std::endl;
		//exit(1);

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

		utility::vector1 <bool> allow_moving( input_pose.size(), false );
		allow_moving.at(pcons->target_res()) = true;

		// find the neighbors for the central residue
		core::scoring::TenANeighborGraph const & graph = input_pose.energies().tenA_neighbor_graph();
		for ( utility::graph::Graph::EdgeListConstIter
				iter = graph.get_node( pcons->target_res() )->const_edge_list_begin(),
				iter_end = graph.get_node( pcons->target_res() )->const_edge_list_end();
				iter != iter_end; ++iter ) {
			Size const neighbor_id( (*iter)->get_other_ind( pcons->target_res() ) );
			allow_moving.at(neighbor_id) = true;
		}


		//  Run this many trajectories
		int const ntraj = option[ OptionKeys::out::nstruct ]();
		for ( int traj = 1; traj <= ntraj; ++traj ) {

			std::cout << "Starting trajectory: " << traj << std::endl;

			// set up BackrubMover and read from the database
			protocols::backrub::BackrubMoverOP backrubmover( new protocols::backrub::BackrubMover() );
			backrubmover->branchopt().read_database();

			// Create Owning Pointer for pose
			pose::PoseOP relax_poseOP( new pose::Pose() );
			*relax_poseOP = input_pose;

			// Setup MC
			protocols::moves::MonteCarloOP mc( new protocols::moves::MonteCarlo(*relax_poseOP, *scorefxn, mc_kT ) );
			mc->reset(*relax_poseOP);

			//clear segments and set the input pose
			backrubmover->clear_segments();
			backrubmover->set_input_pose( relax_poseOP );

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
				segment1 << oneletter_code_from_aa( start_res_neighbor ) << input_pose.pdb_info()->number( j ) << oneletter_code_from_aa( end_res_neighbor );
				std::cout << "Current segment being added: " << segment1.str() << std::endl;
				// add current segment to the mover
				backrubmover->add_segment( start_atom_id, end_atom_id, 0 );
			}

			std::cout << "Number of segments added: " << backrubmover->num_segments() << std::endl;

			std::cout << "Score before backrub: " << relax_poseOP->energies().total_energies()[ total_score ] << std::endl;

			(*scorefxn)(*relax_poseOP);
			// optimize branch angles and idealize side chains
			backrubmover->optimize_branch_angles( *relax_poseOP );
			(*scorefxn)(*relax_poseOP);
			mc->reset(*relax_poseOP);

			std::string move_type = backrubmover->type();

			protocols::moves::TrialMoverOP backrub_trial( new TrialMover (backrubmover, mc) );
			core::pack::task::PackerTaskOP repack_task( pack::task::TaskFactory::create_packer_task( *relax_poseOP ) );
			repack_task->set_bump_check( false );
			repack_task->initialize_from_command_line();
			repack_task->or_include_current( true );
			for ( int ii = 1, nres = input_pose.size(); ii < nres; ++ii ) {
				repack_task->nonconst_residue_task( ii ).restrict_to_repacking();
			}
			repack_task->restrict_to_residues( allow_moving );

			protocols::simple_moves::RotamerTrialsMoverOP pack_rottrial( new protocols::simple_moves::RotamerTrialsMover( scorefxn, *repack_task) );
			protocols::moves::TrialMoverOP pack_rottrial_trial( new TrialMover (pack_rottrial, mc) );

			protocols::moves::SequenceMoverOP repack_step( new SequenceMover );
			repack_step->add_mover( pack_rottrial_trial );


			protocols::moves::CycleMoverOP backrub_cycle( new CycleMover );
			int repack_period_=1;
			for ( Size i=0; (int)i<repack_period_; ++i ) backrub_cycle->add_mover( backrub_trial );
			backrub_cycle->add_mover( repack_step );

			int const ntrials = option[ OptionKeys::pocket_grid::pocket_ntrials ];
			//for ( int trial = 1, ntrials = 1000; trial <= ntrials; ++trial ) {
			for ( int trial = 1; trial <= ntrials; ++trial ) {
				if ( useSA ) mc->set_temperature(mc_kT+500-(trial/ntrials)*500);
				backrub_cycle->apply(*relax_poseOP);
				mc->boltzmann(*relax_poseOP, move_type);
			}
			*relax_poseOP = mc->lowest_score_pose();
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

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}

