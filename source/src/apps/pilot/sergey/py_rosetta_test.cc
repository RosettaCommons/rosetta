// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
///
/// @brief Various small tests thats mirror PyRosetta unit tests
/// @author Sergey Lyskov


// minmover = rosetta.protocols.simple_moves.MinMover();  minmover.score_function(scorefxn)
// movemap = rosetta.MoveMap()
// movemap.set_bb(True)
// #minmover.movemap(movemap)

// minmover.apply(pose)

#include <basic/Tracer.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/MoveMap.hh>

#include <protocols/simple_moves/MinMover.hh>

#include <devel/init.hh>

// C++ implementation of PyRosetta T400_Refinement test
#include <core/fragment/ConstantLengthFragSet.hh>
#include <protocols/simple_moves/FragmentMover.hh>
void T400_Refinement()
{
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap() );
	movemap->set_bb(true);
	movemap->set_bb(10, false);

	core::pose::PoseOP pose_frag = core::import_pose::pose_from_file("src/python/bindings/test/data/test_fragments.pdb", core::import_pose::PDB_file); // pose_frag = pose_from_file("../test/data/test_fragments.pdb")

	core::fragment::ConstantLengthFragSetOP fragset3mer( new core::fragment::ConstantLengthFragSet(3, "src/python/bindings/test/data/test3_fragments") );  // fragset3mer = ConstantLengthFragSet(3, "../test/data/test3_fragments")# "aatestA03_05.200_v1_3")
	core::fragment::ConstantLengthFragSetOP fragset9mer( new core::fragment::ConstantLengthFragSet(9, "src/python/bindings/test/data/test9_fragments") );  // fragset9mer = ConstantLengthFragSet(9, "../test/data/test9_fragments")# "aatestA09_05.200_v1_3")

	movemap->set_bb(1);
	protocols::simple_moves::ClassicFragmentMoverOP mover_3mer( new protocols::simple_moves::ClassicFragmentMover(fragset3mer, movemap) );
	mover_3mer->apply( *pose_frag.get() );
}


int main( int argc, char * argv [] )
{
	using namespace core;

	devel::init(argc, argv);

	T400_Refinement();

	core::pose::Pose pose;
	core::import_pose::pose_from_file(pose, "src/python/bindings/test/data/test_in.pdb", core::import_pose::PDB_file);

	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function_legacy( scoring::PRE_TALARIS_2013_STANDARD_WTS );

	T("Score:") << scorefxn->score(pose)  << std::endl;

	protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover() );
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap() );
	//min_mover->movemap(movemap);
	min_mover->apply(pose);

	T("Scoring done!") << "---------------------" << std::endl;

	return 0;
}
