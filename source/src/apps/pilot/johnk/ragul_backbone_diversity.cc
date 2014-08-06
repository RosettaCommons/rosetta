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
#include <core/conformation/Residue.hh>
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
#include <basic/Tracer.hh>
#include <core/scoring/Energies.hh>

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

OPT_KEY( Integer, relax_start_resnum )
OPT_KEY( Integer, relax_final_resnum )

static basic::Tracer TR( "apps.pilot.ragul_backbone_diversity.main" );

/// General testing code
int
main( int argc, char * argv [] )
{
    try {
	NEW_OPT( relax_start_resnum, "first residue allowed to move", 0 );
	NEW_OPT( relax_final_resnum, "last residues allowed to move", 0 );

	devel::init(argc, argv);

	TR << "Starting to relax backbone as requested" << std::endl;

	core::Real const mc_kT = 0.6;
	TR << "Using kT = " << mc_kT << " for backrub MC...." << std::endl;

	std::string const output_tag = option[ OptionKeys::out::output_tag ]();

	// create pose
	pose::Pose input_pose;

	//read in pdb file from command line
	std::string const input_pdb_name ( basic::options::start_file() );
	core::import_pose::pose_from_pdb( input_pose, input_pdb_name );

	// This is the stretch we'll backrub
	core::Size start_relax_res = 0;
	core::Size end_relax_res = 0;
	int const relax_start_pdb_number = option[ relax_start_resnum ];
	int const relax_final_pdb_number = option[ relax_final_resnum ];

	for ( int j = 1, resnum = input_pose.total_residue(); j <= resnum; ++j ) {
		if ( input_pose.pdb_info()->number(j) == relax_start_pdb_number ) {
			start_relax_res = j;
		}
		if ( input_pose.pdb_info()->number(j) == relax_final_pdb_number ) {
			end_relax_res = j;
		}
	}
	if ( start_relax_res == 0 ) {
		std::cerr << "ERROR!! Could not find residue/chain" << std::endl;
		exit(1);
	}
	if ( end_relax_res == 0 ) {
		std::cerr << "ERROR!! Could not find residue/chain" << std::endl;
		exit(1);
	}

	scoring::ScoreFunctionOP scorefxn( get_score_function() );
	(*scorefxn)(input_pose);
	TR << "Starting score is: " << input_pose.energies().total_energies()[ total_score ] << std::endl;

	//  Run this many trajectories
	int const ntraj = option[ OptionKeys::out::nstruct ]();
	for ( int traj = 1; traj <= ntraj; ++traj ) {

		TR << "Starting trajectory: " << traj << " of " << ntraj << std::endl;

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

		for ( core::Size seg_start = start_relax_res; seg_start < end_relax_res; ++seg_start ) {
			for ( core::Size seg_end = seg_start + 2; seg_end <= end_relax_res; ++seg_end ) {
				id::AtomID start_atom_id = id::AtomID( input_pose.residue( seg_start ).atom_index("CA") , seg_start );
				id::AtomID end_atom_id = id::AtomID( input_pose.residue( seg_end ).atom_index("CA"), seg_end );
				// add segment to the mover
				backrubmover.add_segment( start_atom_id, end_atom_id, 0 );
			}
		}

		TR << "Score before backrub: " << relax_poseOP->energies().total_energies()[ total_score ] << std::endl;

		(*scorefxn)(*relax_poseOP);
		// optimize branch angles and idealize side chains
		backrubmover.optimize_branch_angles( *relax_poseOP );
		(*scorefxn)(*relax_poseOP);
		mc.reset(*relax_poseOP);

		TR << "Done optimizing branch angles" << std::endl;

		std::string move_type = backrubmover.type();
		for ( int trial = 1, ntrials = 2000; trial <= ntrials; ++trial ) {
			backrubmover.apply(*relax_poseOP);
			mc.boltzmann(*relax_poseOP, move_type);
		}
		*relax_poseOP = mc.lowest_score_pose();
		(*scorefxn)(*relax_poseOP);

		TR << "Score after backrub: " << relax_poseOP->energies().total_energies()[ total_score ] << std::endl;

		// writes out the current structure
		std::ostringstream outPDB_name;
		outPDB_name << "traj." << output_tag << "_" << traj << ".pdb";
		relax_poseOP->dump_scored_pdb( outPDB_name.str(), *scorefxn );

	}

	TR << "Successfully finished relaxing backbone" << std::endl;

    } catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
    }
    return 0;
}



