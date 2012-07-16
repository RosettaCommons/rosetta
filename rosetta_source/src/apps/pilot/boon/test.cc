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

#include <core/fragment/ConstantLengthFragSet.hh>
#include <protocols/loops/loop_mover/perturb/LoopMover_CCD.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/docking/DockMCMProtocol.hh>
#include <protocols/simple_moves/RotamerTrialsMinMover.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/kinematics/MoveMap.hh>

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
	pose_from_pdb(test_pose, "/home/boon/data/test.pdb");

	std::cout << "Hello, Rosetta World!" << std::endl;
	std::cout << "I just imported my first pose into Rosetta." << std::endl;
	std::cout << "It has " << test_pose.total_residue() << " total residues." << std::endl;
	
/*	// create a PyMOL mover
	protocols::moves::PyMolMover pmm;
	std::cout << pmm << std::endl;
	
	// create another PyMOL mover
	protocols::moves::PyMolMover pmm2;
	pmm2.update_energy(true);
	pmm2.keep_history(true);
	std::cout << pmm2 << std::endl;*/

/*	// create a score3 scorefxn
	core::scoring::ScoreFunctionCOP scorefxn = new core::scoring::ScoreFunction;
	scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "standard" );
	// read name
	std::cout << "Name: " << scorefxn->get_name() << std::endl;*/

/*	using namespace fragment;
	ConstantLengthFragSetOP fragset = new ConstantLengthFragSet( 3 );
	fragset->read_fragment_file( "home/boon/data/mfr_aa2GB3_03_05.200_v1_3" );
	protocols::simple_moves::ClassicFragmentMover cfmover = protocols::simple_moves::ClassicFragmentMover(fragset);
	std::cout << cfmover << std::endl;*/

	//std::cout << test_pose.energies() << std::endl;

	protocols::simple_moves::SmallMoverOP smallmover = new protocols::simple_moves::SmallMover;
	protocols::simple_moves::ShearMoverOP shearmover = new protocols::simple_moves::ShearMover;
	protocols::moves::SequenceMover seqmover;
	core::Real num1 = 1.0;
	core::Real num2 = 2.0;
	seqmover.add_mover(smallmover, num1);
	seqmover.add_mover(shearmover, num1);
	std::cout << seqmover << std::endl;

	protocols::moves::RandomMover ranmover;
	ranmover.add_mover(shearmover, num1);
	ranmover.add_mover(smallmover, num1);
	ranmover.add_mover(shearmover, num1);
	std::cout << ranmover << std::endl;

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

/*	// print empty rtmover
	protocols::simple_moves::RotamerTrialsMinMover rtmover;
	std::cout << rtmover << std::endl;

	// create a custom scorefxn
	//core::scoring::ScoreFunctionOP scorefxn = new core::scoring::ScoreFunction;
	//scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "score3" );

	// create a standard packer task
	core::pack::task::PackerTaskOP fine_task ( core::pack::task::TaskFactory::create_packer_task( test_pose ));
	fine_task->initialize_from_command_line().restrict_to_repacking().or_include_current( true );

	// define nloop
	//core::Size nloop = 1;

	// print rtmover2
	protocols::simple_moves::RotamerTrialsMinMover rtmover2 = protocols::simple_moves::RotamerTrialsMinMover(scorefxn, *fine_task);
	std::cout << rtmover2 << std::endl;*/

	//protocols::docking::ConformerSwitchMover mover;
	//std::cout << "HOLA" << std::endl;

	//protocols::docking::ConformerSwitchMover mover2 ( protocols::docking::ConformerSwitchMover("True", 1.0) );
	//std::cout << mover2 << std::endl;

/*	// create a loops object
	core::Size start = 15;
	core::Size stop = 24;
	core::Size cutpoint = 19;
	protocols::loops::Loop loop ( protocols::loops::Loop(start, stop, cutpoint) );

	// create a movemap object
	core::kinematics::MoveMapCOP mm;

	// create a CcdLoopClosureMover object
	protocols::loops::loop_closure::ccd::CcdLoopClosureMover ccdmover ( protocols::loops::loop_closure::ccd::CcdLoopClosureMover(loop, mm) );
	std::cout << ccdmover << std::endl;*/

	// create an empty loops object
	//protocols::loops::LoopsOP emptyloops ( new protocols::loops::Loops );

	// create and print a loopmover
	//protocols::loops::loop_mover::refine::LoopMover_Perturb_CCD loopmover;
	//std::cout << loopmover << std::endl;
}
