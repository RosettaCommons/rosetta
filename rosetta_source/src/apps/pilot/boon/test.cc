// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license.
// (c) The Rosetta software is developed by the contributing members of the
// (c) Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org.
// (c) Questions about this can be addressed to University of Washington UW
// (c) TechTransfer, email: license@u.washington.edu.

/// @file   <filename_of_pilot_app>.cc
/// @brief  This is simply a generic pilot app for testing changes.
/// @author <your_name>

// includes
#include <iostream>
#include <string>

#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <protocols/loops/loop_mover/refine/LoopMover_CCD.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/task/PackerTask.hh>

int main(int argc, char *argv[])
{
	using namespace core;
	using namespace import_pose;
	using namespace pose;

	// initialize core
	devel::init(argc, argv);

	// declare variables
	Pose test_pose;

	// import a test pose
	pose_from_pdb(test_pose, "/home/boon/data/1YY9.pdb");

	std::cout << "Hello, Rosetta World!" << std::endl;
	std::cout << "I just imported my first pose into Rosetta." << std::endl;
	std::cout << "It has " << test_pose.total_residue() << " total residues." << std::endl;

	//std::cout << test_pose.energies() << std::endl;
/*	protocols::simple_moves::SmallMoverOP smallmover = new protocols::simple_moves::SmallMover;
	protocols::simple_moves::ShearMoverOP shearmover = new protocols::simple_moves::ShearMover;
	protocols::moves::SequenceMover mover;
	core::Real num = 1.0;
	mover.add_mover(smallmover, num);
	mover.add_mover(shearmover, num);
	std::cout << mover << std::endl;
	core::Size index = 0;
	std::cout << mover.get_mover(index) << std::endl;
	std::cout << "try" << mover.get_mover(index+1) << std::endl;*/

	/*protocols::simple_moves::SwitchResidueTypeSetMover switchmover;
	std::cout << switchmover << std::endl;
	protocols::simple_moves::SwitchResidueTypeSetMover switch2 = protocols::simple_moves::SwitchResidueTypeSetMover("centroid");
	std::cout << switch2 << std::endl;*/

/*	// create a custom scorefxn
	core::scoring::ScoreFunctionOP scorefxn = new core::scoring::ScoreFunction;
	core::Real num = 1.0;
	scorefxn->set_weight("fa_rep", num);
	scorefxn->set_weight("fa_atr", num);

	// create an empty loops object
	protocols::loops::LoopsOP emptyloops = new protocols::loops::Loops;

	// create and print a loopmover
	protocols::loops::loop_mover::refine::LoopMover_Refine_CCD loopmover = protocols::loops::loop_mover::refine::LoopMover_Refine_CCD(emptyloops);
	std::cout << loopmover << std::endl;

	// now set the scorefxn to our custom scorefxn
	protocols::loops::loop_mover::refine::LoopMover_Refine_CCD loopmover2 = protocols::loops::loop_mover::refine::LoopMover_Refine_CCD(emptyloops);
	loopmover2.set_scorefxn(scorefxn);
	std::cout << loopmover2 << std::endl;*/

/*	For packRotamersMover// print empty packmover
	protocols::simple_moves::PackRotamersMover packmover;
	std::cout << packmover << std::endl;

	// create a custom scorefxn
	core::scoring::ScoreFunctionOP scorefxn = new core::scoring::ScoreFunction;
	scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "score3" );

	// create a standard packer task
	core::pack::task::PackerTaskOP fine_task ( core::pack::task::TaskFactory::create_packer_task( test_pose ));
	fine_task->initialize_from_command_line().restrict_to_repacking().or_include_current( true );

	// define nloop
	core::Size nloop = 1;

	// print packmover
	protocols::simple_moves::PackRotamersMover packmover2 = protocols::simple_moves::PackRotamersMover(scorefxn, fine_task, nloop);
	std::cout << packmover2 << std::endl;*/

	protocols::rigid::RigidBodyRandomizeMover mover;
	std::cout << mover << std::endl;
}
