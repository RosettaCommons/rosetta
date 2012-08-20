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
#include <core/pose/annotated_sequence.hh>

#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <protocols/loops/loop_mover/perturb/LoopMover_CCD.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_CCD.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loop_closure/ccd/CcdLoopClosureMover.hh>
#include <protocols/docking/ConformerSwitchMover.hh>
#include <protocols/docking/DockFilters.fwd.hh>
#include <protocols/docking/DockFilters.hh>
#include <protocols/docking/DockingEnsemble.hh>
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
	std::string sequence = "ASDFG";
	core::pose::make_pose_from_sequence(test_pose, sequence, "fa_standard", true);
	// import a test pose
	//pose_from_pdb(test_pose, "/home/boon/data/test.pdb");

	//std::cout << "Hello, Rosetta World!" << std::endl;
	//std::cout << "I just imported my first pose into Rosetta." << std::endl;
	std::cout << "It has " << test_pose.total_residue() << " total residues." << std::endl;

	protocols::docking::ConformerSwitchMover mover;
	std::cout << mover << std::endl;
	
/*	// setup a movemap object
	core::kinematics::MoveMapOP mm ( new core::kinematics::MoveMap );
	Size bb_begin = 3, bb_end = 6, chi_begin = 5, chi_end = 7, begin2 = 101;
	mm->set_bb_true_range(bb_begin, bb_end);
	mm->set_chi_true_range(chi_begin, chi_end);
	mm->set_chi(begin2, true);
  mm->set_bb(begin2, false);
	mm->set_jump(2, true);
	mm->set_jump(5, false);
	// create a standard scorefxn 
	core::scoring::ScoreFunctionCOP scorefxn = new core::scoring::ScoreFunction;
	scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "standard" );
	// setup other inputs
	std::string const min_type_in = "linmin";
	Real tolerance_in = 0.01;
	bool use_nb_list_in = true;
	bool deriv_check_in = false;
	bool deriv_check_verbose_in = false;
	// create a MinMover
	protocols::simple_moves::MinMover minmover = protocols::simple_moves::MinMover(mm, scorefxn, min_type_in, tolerance_in, true, false, false) ;
	// print
	std::cout << "print MinMover:" << std::endl << minmover << std::endl;*/

/*// setup a movemap object
	core::kinematics::MoveMapOP mm ( new core::kinematics::MoveMap );
	Size bb_begin = 3, bb_end = 6, chi_begin = 5, chi_end = 7, begin2 = 101;
	mm->set_bb_true_range(bb_begin, bb_end);
	mm->set_chi_true_range(chi_begin, chi_end);
	mm->set_chi(begin2, true);
	mm->set_jump(2, true);


// DockMCMProtocol
	protocols::docking::DockMCMProtocol dmp;
	std::cout << dmp << std::endl;

	protocols::docking::DockMCMProtocol dmp2;
	dmp2.set_move_map(mm);
	std::cout << dmp2 << std::endl;*/

/*	// setup kT and nmoves
	Real kT = 1.0;
	Size nmoves = 5;
	// create a ShearMover
	protocols::simple_moves::ShearMover shearmover ( protocols::simple_moves::ShearMover(mm, kT, nmoves) );
	// print
	std::cout << "print ShearMover:" << std::endl << shearmover << std::endl;

	// create a SmallMover
	protocols::simple_moves::SmallMover smallmover ( protocols::simple_moves::SmallMover(mm, kT, nmoves) );
	// print
	std::cout << "print SmallMover:" << std::endl << smallmover << std::endl;*/


/*	// create a PyMOL mover
	protocols::moves::PyMolMover pmm;
	std::cout << pmm << std::endl;
	
	// create another PyMOL mover
	protocols::moves::PyMolMover pmm2;
	pmm2.update_energy(true);
	pmm2.keep_history(true);
	std::cout << pmm2 << std::endl;*/

/*	// create a score12 scorefxn
	core::scoring::ScoreFunctionOP scorefxn0 = new core::scoring::ScoreFunction;
	scorefxn0 = core::scoring::ScoreFunctionFactory::create_score_function( "score12" );
	std::cout << "Score function 0's name: " << scorefxn0->get_name() << std::endl;

	// create a standard scorefxn and read fa_sol weight
	core::scoring::ScoreFunctionOP scorefxn = new core::scoring::ScoreFunction;
	scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "standard" );
	enum core::scoring::ScoreType fa_sol = core::scoring::score_type_from_name("fa_sol");
	std::cout << "fa_sol standard weight = " << scorefxn->get_weight(fa_sol) << std::endl;
	// read name
	//std::cout << "Score function's name: " << scorefxn->get_name() << std::endl;

	// now change fa_sol weight
	core::Real num = 1.40;
	scorefxn->set_weight(fa_sol, num);
	std::cout << "fa_sol new weight = " << scorefxn->get_weight(fa_sol) << std::endl;
	// read name again
	std::cout << "Score function's new name: " << scorefxn->get_name() << std::endl;*/

/*	using namespace fragment;
	ConstantLengthFragSetOP fragset = new ConstantLengthFragSet( 3 );
	fragset->read_fragment_file( "home/boon/data/mfr_aa2GB3_03_05.200_v1_3" );
	protocols::simple_moves::ClassicFragmentMover cfmover = protocols::simple_moves::ClassicFragmentMover(fragset);
	std::cout << cfmover << std::endl;*/

	//std::cout << test_pose.energies() << std::endl;

/*	protocols::simple_moves::SmallMoverOP smallmover = new protocols::simple_moves::SmallMover;
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
	std::cout << ranmover << std::endl;*/

	/*protocols::simple_moves::SwitchResidueTypeSetMover switchmover;
	std::cout << switchmover << std::endl;
	protocols::simple_moves::SwitchResidueTypeSetMover switch2 = protocols::simple_moves::SwitchResidueTypeSetMover("centroid");
	std::cout << switch2 << std::endl;*/

/*	// create a custom scorefxn
	core::scoring::ScoreFunctionOP scorefxn = new core::scoring::ScoreFunction;
	core::Real num = 1.0;
	scorefxn->set_weight("fa_rep", num);
	scorefxn->set_weight("fa_atr", num);*/

/*	// create an empty loops object
	protocols::loops::LoopsOP emptyloops = new protocols::loops::Loops;

	// create and print a loopmover
	protocols::loops::loop_mover::refine::LoopMover_Refine_CCD loopmover = protocols::loops::loop_mover::refine::LoopMover_Refine_CCD(emptyloops);
	std::cout << loopmover << std::endl;

	// create a loops object
	core::Size start = 15, start2 = 51;
	core::Size stop = 24, stop2 = 60;
	core::Size cutpoint = 19; 
	protocols::loops::Loop loop ( protocols::loops::Loop(start, stop, cutpoint) );
	protocols::loops::Loop loop2 ( protocols::loops::Loop(start2, stop2) );
	protocols::loops::LoopsOP loops = new protocols::loops::Loops;
	loops->add_loop(loop);
	loops->add_loop(loop2);	

	// create an print the new loopmover
	protocols::loops::loop_mover::refine::LoopMover_Refine_CCD loopmover2 = protocols::loops::loop_mover::refine::LoopMover_Refine_CCD(loops);
	loopmover2.move_map(mm);
	std::cout << loopmover2 << std::endl;*/

/*	// print empty rtmover
	protocols::simple_moves::RotamerTrialsMinMover rtmover;
	std::cout << rtmover << std::endl;

	// create a custom scorefxn
	//core::scoring::ScoreFunctionOP scorefxn = new core::scoring::ScoreFunction;
	//scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "score3" );*/

/*	// setup a packer task
	core::pack::task::PackerTaskOP task ( core::pack::task::TaskFactory::create_packer_task( test_pose ));
	task->initialize_from_command_line().or_include_current( true );
	Size resid = 2;
	std::cout << "show residue task for resid 2:" << std::endl;
	task->show_residue_task(resid);
	std::cout << "show all residue tasks" << std::endl;
	task->show_all_residue_tasks();
	std::cout << "print packer task" << std::endl << *task << std::endl;	*/

/*	// define nloop
	//core::Size nloop = 1;

	// print rtmover2
	protocols::simple_moves::RotamerTrialsMinMover rtmover2 = protocols::simple_moves::RotamerTrialsMinMover(scorefxn, *fine_task);
	std::cout << rtmover2 << std::endl;*/

	//protocols::docking::ConformerSwitchMover mover;
	//std::cout << "HOLA" << std::endl;

	//protocols::docking::ConformerSwitchMover mover2 ( protocols::docking::ConformerSwitchMover("True", 1.0) );
	//std::cout << mover2 << std::endl;

  /*// create a loops object
	core::Size start = 15;
	core::Size stop = 24;
	core::Size cutpoint = 19;
	protocols::loops::Loop loop ( protocols::loops::Loop(start, stop, cutpoint) );

  // setup a movemap object
	core::kinematics::MoveMapOP mm (new core::kinematics::MoveMap );
	Size bb_begin = 3, bb_end = 6, chi_begin = 5, chi_end = 7;
	mm->set_bb_true_range(bb_begin, bb_end);
	mm->set_chi_true_range(chi_begin, chi_end);
	mm->set_jump(2, true);

	// create a CcdLoopClosureMover object
	protocols::loops::loop_closure::ccd::CcdLoopClosureMover ccdmover ( protocols::loops::loop_closure::ccd::CcdLoopClosureMover(loop, mm) );
	std::cout << ccdmover << std::endl;*/

	// create an empty loops object
	// protocols::loops::LoopsOP emptyloops ( new protocols::loops::Loops );

	// create and print a loopmover
	//protocols::loops::loop_mover::refine::LoopMover_Perturb_CCD loopmover;
	//std::cout << loopmover << std::endl;
}
