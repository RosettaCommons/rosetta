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
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AA.hh>
#include <core/types.hh>
#include <core/conformation/Residue.hh>

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
#include <core/pack/pack_rotamers.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/docking/DockMCMProtocol.hh>
#include <protocols/farna/RNA_ProtocolUtil.hh>

int main(int argc, char *argv[])
{
    try {

	using namespace core;
	using namespace import_pose;
	using namespace pose;

	// initialize core
	devel::init(argc, argv);

	// declare variables
	Pose test_pose;
	//std::string sequence = "ASDFG";
	//core::pose::make_pose_from_sequence(test_pose, sequence, "fa_standard", true);
	// import a test pose
	pose_from_pdb(test_pose, "/home/boon/data/test.pdb");

	//std::cout << "Hello, Rosetta World!" << std::endl;
	//std::cout << "I just imported my first pose into Rosetta." << std::endl;
	//std::cout << "It has " << test_pose.total_residue() << " total residues." << std::endl;
	/*chemical::ResidueTypeSetCAP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );
	pose_from_pdb(test_pose, *rsd_set, "/home/boon/data/test_rna.pdb");
	protocols::farna::make_phosphate_nomenclature_matches_mini( test_pose );
	// setup a packer task
	pack::task::PackerTaskOP task ( core::pack::task::TaskFactory::create_packer_task( test_pose ) );
	task->initialize_from_command_line().or_include_current( true );
	for (Size ii = 1; ii <= test_pose.total_residue(); ++ii ) {
		task->nonconst_residue_task(ii).allow_aa( na_rad );
		task->nonconst_residue_task(ii).allow_aa( na_ura );
		task->nonconst_residue_task(ii).allow_aa( na_rgu );
		task->nonconst_residue_task(ii).allow_aa( na_rcy );
		assert( task->design_residue(ii) );
	}
  // create a farna/rna_lores scorefxn
	scoring::ScoreFunction scorefxn;
	scorefxn = scoring::ScoreFunctionFactory::create_score_function( "farna/rna_lores" );
	Size const nloop = 3;
	utility::vector1< Real > score_list;
	utility::vector1< std::string > seq_list;
	utility::vector1< pose::PoseOP > pose_list;
	// apply pack rotamers loop
	pack::pack_rotamers_loop( pose, scorefxn, task, nloop, score_list, seq_list, pose_list );

	for (Size n = 1; n <= pose_list.size() ; n++ ){
		std::cout << score_list[n] << " " << seq_list[n] << std::endl;
	}*/

/*	// create a FoldTree
	kinematics::FoldTree ft;
	ft.add_edge(1,13,-1);
	ft.add_edge(13,19,-1);
	ft.add_edge(13,26,1);
	ft.add_edge(26,20,-1);
	ft.add_edge(26,116,-1);
	std::cout << ft << std::endl;
	ft.show();

	// create a FoldTree (new method)
	kinematics::FoldTree ft2;
	ft2.add_edge(1,13);
	ft2.add_edge(13,19);
	ft2.add_jump(13,26,1);
	ft2.add_edge(26,20);
	ft2.add_edge(26,116);
	std::cout << ft2 << std::endl;
	ft2.show();*/

/*	//Reformatting

	// PyMOL mover
	protocols::moves::PyMolMover pmm;
	pmm.update_energy(true);
	pmm.keep_history(true);
	std::cout << pmm << std::endl;

	// ReturnSidechainMover
	protocols::simple_moves::ReturnSidechainMover returnmover = protocols::simple_moves::ReturnSidechainMover(test_pose,10,20) ;
	std::cout << returnmover << std::endl;

	// RigidBodyPerturbMover
	protocols::rigid::RigidBodyPerturbMover pmover;
	std::cout << pmover << std::endl;

	// RigidBodyRandomizeMover
	protocols::rigid::RigidBodyRandomizeMover rmover;
	std::cout << rmover << std::endl;

	// RigidBodySpinMover
	protocols::rigid::RigidBodySpinMover smover (protocols::rigid::RigidBodySpinMover(1));
	std::cout << smover << std::endl;

	// RigidBodyTransMover
	int jump_num = 2;
	protocols::rigid::RigidBodyTransMover mover2 (protocols::rigid::RigidBodyTransMover(test_pose, jump_num));
	std::cout << mover2 << std::endl;

	// TrialMover
	protocols::moves::TrialMover trialmover;
	std::cout << trialmover << std::endl;

	// DockingSlideIntoContact
	protocols::docking::DockingSlideIntoContact dmover (protocols::docking::DockingSlideIntoContact(2) );
	std::cout << dmover << std::endl;

	// FaDockingSlideIntoContact
	protocols::docking::FaDockingSlideIntoContact famover (protocols::docking::FaDockingSlideIntoContact(1) );
	std::cout << famover << std::endl;*/

/*	protocols::docking::ConformerSwitchMover mover;
	std::cout << mover << std::endl;*/

 // setup a movemap object
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
		scorefxn = core::scoring::getScoreFunctionLegacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );
	// setup other inputs
	std::string const min_type_in = "linmin";
	Real tolerance_in = 0.01;
	//bool use_nb_list_in = true;  // unused ~Labonte
	//bool deriv_check_in = false;  // unused ~Labonte
	//bool deriv_check_verbose_in = false;  // unused ~Labonte
	// create a MinMover
	protocols::simple_moves::MinMover minmover = protocols::simple_moves::MinMover(mm, scorefxn, min_type_in, tolerance_in, true, false, false) ;
	// print
	std::cout << "print MinMover:" << std::endl << minmover << std::endl;

/* // setup a movemap object
	core::kinematics::MoveMapOP mm ( new core::kinematics::MoveMap );
	Size bb_begin = 3, bb_end = 6, chi_begin = 5, chi_end = 7, begin2 = 101;
	mm->set_bb_true_range(bb_begin, bb_end);
	mm->set_chi_true_range(chi_begin, chi_end);
	mm->set_chi(begin2, true);
	mm->set_jump(2, true);*/


/*// DockMCMProtocol
	protocols::docking::DockMCMProtocol dmp;
	std::cout << dmp << std::endl;

	protocols::docking::DockMCMProtocol dmp2;
	dmp2.set_move_map(mm);
	std::cout << dmp2 << std::endl;*/

	// setup kT and nmoves
	Real kT = 1.0;
	Size nmoves = 5;
	// create a ShearMover
	protocols::simple_moves::ShearMover shearmover ( protocols::simple_moves::ShearMover(mm, kT, nmoves) );
	// print
	std::cout << "print ShearMover:" << std::endl << shearmover << std::endl;

	// create a SmallMover
	protocols::simple_moves::SmallMover smallmover ( protocols::simple_moves::SmallMover(mm, kT, nmoves) );
	// print
	std::cout << "print SmallMover:" << std::endl << smallmover << std::endl;


/*	// create a PyMOL mover
	protocols::moves::PyMolMover pmm;
	std::cout << pmm << std::endl;

	// create another PyMOL mover
	protocols::moves::PyMolMover pmm2;
	pmm2.update_energy(true);
	pmm2.keep_history(true);
	std::cout << pmm2 << std::endl;*/

/*	// create a scorefxn
	core::scoring::ScoreFunctionOP scorefxn0 = core::scoring::getScoreFunction();
	std::cout << "Score function 0's name: " << scorefxn0->get_name() << std::endl;

	// create a standard scorefxn and read fa_sol weight
	core::scoring::ScoreFunctionOP scorefxn = new core::scoring::ScoreFunction;
	scorefxn = core::scoring::getScoreFunctionLegacy( "pre_talaris_2013_standard" );
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

/*	Size frag_length = 3;
	core::fragment::ConstantLengthFragSetOP fragset = new core::fragment::ConstantLengthFragSet( frag_length, "/home/boon/data/test_in.frag3" );
	protocols::simple_moves::ClassicFragmentMover cfmover = protocols::simple_moves::ClassicFragmentMover(fragset, mm);
	std::cout << cfmover << std::endl;

	protocols::simple_moves::ClassicFragmentMover cfmover2 = protocols::simple_moves::ClassicFragmentMover(fragset);
	std::cout << cfmover2 << std::endl;*/

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

/*  // create an empty loops object
	protocols::loops::LoopsOP emptyloops = new protocols::loops::Loops;

	// create and print a loopmover
	protocols::loops::loop_mover::refine::LoopMover_Refine_CCD loopmover = protocols::loops::loop_mover::refine::LoopMover_Refine_CCD(emptyloops);
	std::cout << loopmover << std::endl;*/

/*  // create a loops object
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
	//scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "score3" );

	// setup a packer task
	core::pack::task::PackerTaskOP task ( core::pack::task::TaskFactory::create_packer_task( test_pose ));
	task->initialize_from_command_line().or_include_current( true );
	Size resid = 2;
	std::cout << "show residue task for resid 2:" << std::endl;
	task->show_residue_task(resid);
	std::cout << "show all residue tasks" << std::endl;
	task->show_all_residue_tasks();
	std::cout << "print packer task" << std::endl << *task << std::endl;*/

/*	// define nloop
	//core::Size nloop = 1;

	// print rtmover2
	protocols::simple_moves::RotamerTrialsMinMover rtmover2 = protocols::simple_moves::RotamerTrialsMinMover(scorefxn, *fine_task);
	std::cout << rtmover2 << std::endl;*/

	//protocols::docking::ConformerSwitchMover mover;
	//std::cout << "HOLA" << std::endl;

	//protocols::docking::ConformerSwitchMover mover2 ( protocols::docking::ConformerSwitchMover("True", 1.0) );
	//std::cout << mover2 << std::endl;

// create a loops object
	//core::Size start = 15;
	//core::Size stop = 24;
	//core::Size cutpoint = 19;
	//protocols::loops::Loop loop ( protocols::loops::Loop(start, stop, cutpoint) );

  // setup a movemap object
	//core::kinematics::MoveMapOP mm (new core::kinematics::MoveMap );
	//Size bb_begin = 3, bb_end = 6, chi_begin = 5, chi_end = 7;
	//mm->set_bb_true_range(bb_begin, bb_end);
	//mm->set_chi_true_range(chi_begin, chi_end);
	//mm->set_jump(2, true);

/*	// create a CcdLoopClosureMover object
	protocols::loops::loop_closure::ccd::CcdLoopClosureMover ccdmover ( protocols::loops::loop_closure::ccd::CcdLoopClosureMover(loop, mm) );
	std::cout << ccdmover << std::endl;*/

	// create an empty loops object
	// protocols::loops::LoopsOP emptyloops ( new protocols::loops::Loops );
	// create a scorefxn
/*	core::scoring::ScoreFunctionOP scorefxn = core::scoring::getScoreFunction();
	// create a loops object
	core::Size start = 15, start2 = 51;
	core::Size stop = 24, stop2 = 60;
	core::Size cutpoint = 19, cutpoint2 = 55;
	protocols::loops::Loop loop ( protocols::loops::Loop(start, stop, cutpoint) );
	protocols::loops::Loop loop2 ( protocols::loops::Loop(start2, stop2, cutpoint2) );
	protocols::loops::LoopsOP loops = new protocols::loops::Loops;
	loops->add_loop(loop);
	loops->add_loop(loop2);
	// create and print a loopmover
	protocols::loops::loop_mover::perturb::LoopMover_Perturb_CCD loopmover;
	std::cout << loopmover << std::endl;
	// create another loopmover (add loops)
	protocols::loops::loop_mover::perturb::LoopMover_Perturb_CCD loopmover2 (protocols::loops::loop_mover::perturb::LoopMover_Perturb_CCD(loops, scorefxn));
	std::cout << loopmover2 << std::endl;*/

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
}
