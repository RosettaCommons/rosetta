// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief
/// @author jk

#include <iostream>
#include <iomanip>

#include <devel/init.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID_Map.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/TwelveANeighborGraph.hh>
#include <core/scoring/Energies.hh>

#include <protocols/moves/RigidBodyMover.hh>
#include <protocols/moves/BackrubMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/KinematicMover.hh>
#include <protocols/loops/KinematicWrapper.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>

#include <basic/Tracer.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/after_opts.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>


using namespace core;
using namespace basic::options;
using namespace core::scoring;
using namespace core::optimization;
using namespace basic::options::OptionKeys;

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.chet_trp_to_gly.main" );


void do_minimization(
										 core::pose::Pose & pose,
										 utility::vector1 <bool> & allow_moving,
										 pack::task::PackerTaskCOP repack_task,
										 scoring::ScoreFunctionCOP scorefxn,
										 scoring::ScoreFunctionCOP repack_scorefxn
) {

		// setting degrees of freedom which can move during minimization
	kinematics::MoveMap mm;
	mm.set_jump( false );
	for ( core::Size j = 1; j <= pose.total_residue(); ++j ) {
		mm.set_chi( j, allow_moving.at(j) );
		mm.set_bb( j, allow_moving.at(j) );
	}

	// minimize protein
	TR << "Starting minimization...." << std::endl;
	AtomTreeMinimizer minimizer;
	MinimizerOptions min_options( "dfpmin", 0.00001, true, false );

	minimizer.run( pose, mm, *scorefxn, min_options );
	(*scorefxn)(pose);

	// calculate score and print it out
	(*scorefxn)(pose);
	TR << "Post minimization 1 score: " << pose.energies().total_energies()[ total_score ] << std::endl;

	// First repack
	pack::pack_rotamers( pose, *repack_scorefxn, repack_task );

	// Report Scores
	(*scorefxn)(pose);
	//	pose.dump_scored_pdb( "repacked_once.pdb", *scorefxn );
	TR << "Score after repacking once: " << pose.energies().total_energies()[ total_score ] << std::endl << std::endl;

	// iterate over minimizing and repacking
	Size const num_iter = 1;
	for ( Size iter = 1; iter <= num_iter; ++iter ) {
		minimizer.run( pose, mm, *scorefxn, min_options );
		(*scorefxn)(pose);
		TR << "Current score after minimizing: " << pose.energies().total_energies()[ total_score ] << std::endl << std::endl;
		pack::pack_rotamers( pose, *repack_scorefxn, repack_task );
		(*scorefxn)(pose);
		TR << "Current score after repacking: " << pose.energies().total_energies()[ total_score ] << std::endl << std::endl;
	}

	// final minimization
	minimizer.run( pose, mm, *scorefxn, min_options );
	(*scorefxn)(pose);

	return;
}


void do_KIC(
										 core::pose::Pose & pose,
										 utility::vector1 <bool> & allow_moving,
										 pack::task::PackerTaskCOP repack_task,
										 scoring::ScoreFunctionCOP scorefxn,
										 scoring::ScoreFunctionCOP repack_scorefxn
) {

	(*scorefxn)(pose);
	//	pose.dump_scored_pdb( "init_pose.pdb", *scorefxn );
	TR << "Initial score: " << pose.energies().total_energies()[ total_score ] << std::endl;

	// figure out loop regions from allow_moving, build a KIC mover for each
	protocols::moves::RandomMoverOP loop_move_set( new protocols::moves::RandomMover() );
	for ( core::Size j = 2; j <= pose.total_residue(); ++j ) {
		if ( allow_moving.at(j) ) {
			core::Size loop_begin = j;
			core::Size loop_end = j;
			for ( core::Size k = j; k <= pose.total_residue(); ++k ) {
				if ( ! allow_moving.at(k) ) {
					break;
				}
				loop_end = k;
			}

			if ( ( loop_end - loop_begin + 1 ) >= 3 ) {
				protocols::moves::KinematicMoverOP kin_mover( new protocols::moves::KinematicMover() );
				kin_mover->set_temperature( 0.8 );
				kin_mover->set_vary_bondangles( true );
				kin_mover->set_sample_nonpivot_torsions( true );
				kin_mover->set_rama_check( true );
				protocols::loops::KinematicWrapperOP kin_wrapper( new protocols::loops::KinematicWrapper(kin_mover, loop_begin, loop_end));
				loop_move_set->add_mover(kin_wrapper, 5);
			}

			j=loop_end;
		}
	}


	Size const num_iter = 5;
	for ( Size iter = 1; iter <= num_iter; ++iter ) {

		// Setup MC, run KIC
		Real const mc_kT = 10. * ( num_iter + 1 - iter );
		protocols::moves::MonteCarlo mc(pose, *scorefxn, mc_kT );
		mc.reset(pose);
		(*scorefxn)(pose);
		// run the backrub
		std::string move_type = loop_move_set->type();
		for ( int trial = 1, ntrials = 200; trial <= ntrials; ++trial ) {
			loop_move_set->apply(pose);
			mc.boltzmann(pose, move_type);
		}

		//		if ( iter == num_iter ) pose = mc.lowest_score_pose(); // if it's the last iteration, keep the best we've seen

		(*scorefxn)(pose);
		TR << "Current score after kinematic loop closure: " << pose.energies().total_energies()[ total_score ] << std::endl << std::endl;
		// repack
		pack::pack_rotamers( pose, *repack_scorefxn, repack_task );
		(*scorefxn)(pose);
		TR << "Current score after repacking: " << pose.energies().total_energies()[ total_score ] << std::endl << std::endl;
	}

	// do a minimization (using the previous code)
	//	do_minimization( pose, allow_moving, repack_task, scorefxn, repack_scorefxn );
	TR << "Score after final minimization: " << pose.energies().total_energies()[ total_score ] << std::endl;

	return;

}


void do_backrub(
										 core::pose::Pose & pose,
										 utility::vector1 <bool> & allow_moving,
										 pack::task::PackerTaskCOP repack_task,
										 scoring::ScoreFunctionCOP scorefxn,
										 scoring::ScoreFunctionCOP repack_scorefxn
) {

	(*scorefxn)(pose);
	//	pose.dump_scored_pdb( "init_pose.pdb", *scorefxn );
	TR << "Initial score: " << pose.energies().total_energies()[ total_score ] << std::endl;

	// set up BackrubMover and read from the database
	protocols::moves::BackrubMover backrubmover;
	backrubmover.branchopt().read_database();

	//clear segments and set the input pose
	backrubmover.clear_segments();
	pose::PoseOP tmp_poseOP ( new pose::Pose() );
	*tmp_poseOP = pose;
	backrubmover.set_input_pose( tmp_poseOP );

	// generate segments from allow_moving
	for ( core::Size j = 2; j <= pose.total_residue(); ++j ) {
		if ( allow_moving.at(j) ) {
			core::Size start_res = j-1;
			core::Size end_res = j;
			// find contiguous segments of allow_moving
			for ( core::Size k = j; k <= pose.total_residue(); ++k ) {
				if ( ! allow_moving.at(k) ) {
					break;
				}
				end_res = k;
			}
			// make all possible segments in this range of residues
			for ( core::Size seg_start = start_res; seg_start < end_res; ++seg_start ) {
				for ( core::Size seg_end = seg_start + 2; seg_end <= end_res; ++seg_end ) {
					id::AtomID start_atom_id = id::AtomID( pose.residue( seg_start ).atom_index("CA") , seg_start );
					id::AtomID end_atom_id = id::AtomID( pose.residue( seg_end ).atom_index("CA"), seg_end );
					// add segment to the mover
					backrubmover.add_segment( start_atom_id, end_atom_id, 0 );
				}
			}
		}
	}

	// optimize branch angles and idealize side chains
	backrubmover.optimize_branch_angles( pose );
	(*scorefxn)(pose);

	Size const num_iter = 5;
	for ( Size iter = 1; iter <= num_iter; ++iter ) {

		// reset backrub input pose
		pose::PoseOP tmp_poseOP ( new pose::Pose() );
		*tmp_poseOP = pose;
		backrubmover.set_input_pose( tmp_poseOP );

		// Setup MC, run backrub
		Real const mc_kT = 20;
		//Real const mc_kT = 0.3 * ( num_iter + 1 - iter );
		protocols::moves::MonteCarlo mc(pose, *scorefxn, mc_kT );
		mc.reset(pose);
		(*scorefxn)(pose);
		// run the backrub
		std::string move_type = backrubmover.type();
		for ( int trial = 1, ntrials = 200; trial <= ntrials; ++trial ) {
			backrubmover.apply(pose);
			mc.boltzmann(pose, move_type);
		}
		pose = mc.lowest_score_pose();

		(*scorefxn)(pose);
		TR << "Current score after backrub: " << pose.energies().total_energies()[ total_score ] << std::endl << std::endl;
		// repack
		pack::pack_rotamers( pose, *repack_scorefxn, repack_task );
		(*scorefxn)(pose);
		TR << "Current score after repacking: " << pose.energies().total_energies()[ total_score ] << std::endl << std::endl;
	}

	// do a minimization (using the previous code)
	do_minimization( pose, allow_moving, repack_task, scorefxn, repack_scorefxn );
	TR << "Score after final minimization: " << pose.energies().total_energies()[ total_score ] << std::endl;

	return;

}


OPT_KEY( Integer, mut_resnum )
OPT_KEY( String, mut_chain )
OPT_KEY( Boolean, mut_to_ala )


/// General testing code
int
main( int argc, char * argv [] )
{

	NEW_OPT( mut_resnum, "which residue to mutate to Gly", 1 );
	NEW_OPT( mut_chain, "which chain the residue is on", "A" );
	NEW_OPT( mut_to_ala, "whether to mutate to Ala instead of Gly", false );

	devel::init(argc, argv);

	TR << "Starting chet_trp_to_gly" << std::endl;

	// JK
	// scoring function
	scoring::ScoreFunctionOP scorefxn( get_score_function() );
	scoring::ScoreFunctionOP repack_scorefxn( get_score_function() );
	//	scoring::ScoreFunctionOP repack_scorefxn( ScoreFunctionFactory::create_score_function(SOFT_REP_WTS) );
	//	scoring::ScoreFunctionOP repack_scorefxn( ScoreFunctionFactory::create_score_function(SOFT_REP_DESIGN_WTS) );

	//	repack_scorefxn->set_weight( core::scoring::fa_dun, 0.1 );

	std::string const tmp_chain = option[ mut_chain ];
	if ( tmp_chain.length() != 1 ) {
		TR << "ERROR!! Chain ID should be one character" << std::endl;
		exit(1);
	}
	char const mut_pdb_chain = tmp_chain[0];
	int const mut_pdb_number = option[ mut_resnum ];

	// create pose for native pose from pdb
	pose::Pose wt_pose_init;
	std::string const input_pdb_name ( basic::options::start_file() );
	core::import_pose::pose_from_pdb( wt_pose_init, input_pdb_name );
	(*scorefxn)( wt_pose_init );

	// set mut_rosetta_resnum to Rosetta internal resid for the residue to be mutated
	core::Size mut_rosetta_resnum = 0;
	for ( core::Size j = 1; j <= wt_pose_init.total_residue(); ++j ) {
		if ( ( wt_pose_init.pdb_info()->chain(j) == mut_pdb_chain ) &&
				 ( wt_pose_init.pdb_info()->number(j) == mut_pdb_number ) ) {
					 mut_rosetta_resnum = j;
		}
	}
	if ( mut_rosetta_resnum == 0 ) {
		TR << "ERROR!! Could not find residue/chain" << std::endl;
		exit(1);
	}

	if ( chemical::oneletter_code_from_aa( wt_pose_init.residue(mut_rosetta_resnum).aa() ) != 'W' ) {
		TR << "ERROR!! Mutation site is not a Trp in the WT structure" << std::endl;
		TR << "This is not a requirement for the code, but too serious for just a warning....." << std::endl;
		exit(1);
	}

	// build the glycine (or alanine) mutant
	pose::Pose mut_pose_init = wt_pose_init;
	// setup a packer task
	pack::task::PackerTaskOP mut_task( pack::task::TaskFactory::create_packer_task( mut_pose_init ));
	mut_task->set_bump_check( false );
	mut_task->initialize_from_command_line();
	mut_task->or_include_current( true );
	// restrict packer task to single sequence position of interest
	int const aa_ala = 1;
	int const aa_gly = 6;
	utility::vector1<bool> allow_redesign( wt_pose_init.total_residue(), false );
	allow_redesign.at(mut_rosetta_resnum) = true;
	utility::vector1< bool > mut_site( core::chemical::num_canonical_aas, false );
	if ( option[ mut_to_ala ] ) {
		mut_site.at( aa_ala ) = true;
	} else {
		mut_site.at( aa_gly ) = true;
	}
	mut_task->restrict_to_residues( allow_redesign );
	mut_task->nonconst_residue_task(mut_rosetta_resnum).restrict_absent_canonical_aas( mut_site );
	pack::pack_rotamers( mut_pose_init, *repack_scorefxn, mut_task );
	mut_pose_init.dump_scored_pdb( "start_mut_"+input_pdb_name, *repack_scorefxn );


	// setup packer task to be used in minimization protocols (repack everything within 10 A of mutation site)
	pack::task::PackerTaskOP task_repack_wt( pack::task::TaskFactory::create_packer_task( wt_pose_init ));
	pack::task::PackerTaskOP task_repack_mut( pack::task::TaskFactory::create_packer_task( mut_pose_init ));
	task_repack_wt->set_bump_check( true );
	task_repack_mut->set_bump_check( true );
	task_repack_wt->initialize_from_command_line();
	task_repack_mut->initialize_from_command_line();
	task_repack_wt->or_include_current( true );
	task_repack_mut->or_include_current( true );
	// restrict repacking to the neighbors of the mutation site
	utility::vector1 <bool>	allow_moving( wt_pose_init.total_residue(), false );
	allow_moving.at( mut_rosetta_resnum ) = true;
	core::scoring::TwelveANeighborGraph const & graph = wt_pose_init.energies().twelveA_neighbor_graph();
	for ( core::graph::Graph::EdgeListConstIter
					iter = graph.get_node( mut_rosetta_resnum )->const_edge_list_begin(),
					iter_end = graph.get_node( mut_rosetta_resnum )->const_edge_list_end();
				iter != iter_end; ++iter ) {
		Size const neighbor_id( (*iter)->get_other_ind( mut_rosetta_resnum ) );
		allow_moving.at(neighbor_id) = true;
	}
	for (core::Size ii = 1; ii < wt_pose_init.total_residue(); ++ii ) {
		task_repack_wt->nonconst_residue_task( ii ).restrict_to_repacking();
		task_repack_mut->nonconst_residue_task( ii ).restrict_to_repacking();
	}
	task_repack_wt->restrict_to_residues( allow_moving );
	task_repack_mut->restrict_to_residues( allow_moving );

	for ( int ii = 1; ii <= option[ OptionKeys::out::nstruct ](); ++ii ) {

		pose::Pose wt_pose = wt_pose_init;
		TR << "Carrying out relax of the WT #" << ii << std::endl;
		//		do_KIC( wt_pose, allow_moving, task_repack_wt, scorefxn, repack_scorefxn );
		std::stringstream wt_fname;
		wt_fname << "min_wt_" << ii << "_" << input_pdb_name;
		wt_pose.dump_scored_pdb( wt_fname.str(), *scorefxn );
		TR << "Final score WT (#" << ii << "): " << wt_pose.energies().total_energies()[ total_score ] << std::endl << std::endl;

		pose::Pose mut_pose = mut_pose_init;
		TR << "Carrying out relax of the MUTANT #" << ii << std::endl;
		do_KIC( mut_pose, allow_moving, task_repack_mut, scorefxn, repack_scorefxn );
		std::stringstream mut_fname;
		mut_fname << "min_mut_" << ii << "_" << input_pdb_name;
		mut_pose.dump_scored_pdb( mut_fname.str(), *scorefxn );
		TR << "Final score Mut (#" << ii << "): " << mut_pose.energies().total_energies()[ total_score ] << std::endl << std::endl;

	}

	TR << "Done!!" << std::endl;

	return 0;

}


