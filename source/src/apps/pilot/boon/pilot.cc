// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test.cc
/// @brief  This is simply a generic pilot app for testing changes.
/// @author Boon Uranukul

// includes
#include <iostream>
#include <string>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/sequence/Sequence.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AA.hh>
#include <core/types.hh>
#include <core/conformation/Residue.hh>

#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/LK_PolarNonPolarEnergy.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/farna/util.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_KIC.hh>
#include <protocols/loops/loop_mover/perturb/LoopMover_KIC.hh>
#include <protocols/farna/RNA_DeNovoProtocol.hh>
#include <protocols/farna/RNA_Minimizer.hh>

int main(int argc, char *argv[])
{
	try {

	using namespace core;
	using namespace import_pose;
	using namespace pose;

	// initialize core
	devel::init(argc, argv);

	// declare variables
	//Pose test_pose;

	///////////////////////////////Arguments Setup////////////////////////////////////////////////////////
	// create a score12 scorefxn
	scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
	scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "score12" );
	// create a loops object
	Size start = 15, start2 = 51;
	Size stop = 24, stop2 = 60;
	Size cutpoint = 19, cutpoint2 = 55;
	protocols::loops::Loop loop ( protocols::loops::Loop(start, stop, cutpoint) );
	protocols::loops::Loop loop2 ( protocols::loops::Loop(start2, stop2, cutpoint2) );
	protocols::loops::LoopsOP loops( new protocols::loops::Loops );
	loops->add_loop(loop);
	loops->add_loop(loop2);

	/////////////////////////////////////RNA_Minimizer///////////////////////////////////////////////////
	// minimizer setup
	protocols::farna::RNA_Minimizer rna_minimizer;
	rna_minimizer.deriv_check( true );
	rna_minimizer.use_coordinate_constraints( false );
	rna_minimizer.skip_o2prime_trials( true );
	rna_minimizer.vary_bond_geometry( true );
	std::cout << "\nPrint RNA_Minimizer:" << std::endl;
	std::cout << rna_minimizer << std::endl;

	/////////////////////////////////////RNA_DeNovoProtocol//////////////////////////////////////////////
	Size const nstruct = 1;
	std::string const silent_file = "output.txt";
	bool const heat_structure( true );
	bool const minimize_structure = false;
	bool const relax_structure = false;
	bool const is_allow_bulge = false;

  protocols::farna::RNA_DeNovoProtocol rna_de_novo_protocol( nstruct,
                                                             silent_file,
                                                             heat_structure,
                                                             minimize_structure,
                                                             relax_structure,
                                                             is_allow_bulge );

	std::cout << "\nPrint RNA_DeNovoProtocol:" << std::endl;
  std::cout << rna_de_novo_protocol << std::endl;

	//////////////////////////////LoopMover_Perturb_KIC//////////////////////////////////////////////////
	// create and print a KIC perturb loopmover
	protocols::loops::loop_mover::perturb::LoopMover_Perturb_KIC loopmover;
	std::cout << "\nPrint LoopMover Perturb KIC (w/o argument):" << std::endl;
	std::cout << loopmover << std::endl;
	std::cout << "\nPrint LoopMover Perturb KIC (with argument):" << std::endl;
	// create another KIC perturb loopmover (add loops)
	protocols::loops::loop_mover::perturb::LoopMover_Perturb_KIC loopmover2 (protocols::loops::loop_mover::perturb::LoopMover_Perturb_KIC(loops, scorefxn));
	std::cout << loopmover2 << std::endl;

	//////////////////////////////LoopMover_Refine_KIC//////////////////////////////////////////////////
	// create and print a KIC refine loopmover
	protocols::loops::loop_mover::refine::LoopMover_Refine_KIC loopmover3;
	std::cout << "\nPrint LoopMover Refine KIC (w/o argument):" << std::endl;
	std::cout << loopmover3 << std::endl;
	std::cout << "\nPrint LoopMover Refine KIC (with argument):" << std::endl;
	// create another KIC refine loopmover (add loops)
	protocols::loops::loop_mover::refine::LoopMover_Refine_KIC loopmover4 (protocols::loops::loop_mover::refine::LoopMover_Refine_KIC(loops, scorefxn));
	std::cout << loopmover4 << std::endl;

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
