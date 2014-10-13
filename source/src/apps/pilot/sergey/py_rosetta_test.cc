// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

int main( int argc, char * argv [] )
{
	using basic::T;
	using namespace core;

	devel::init(argc, argv);

	core::pose::Pose pose;
	core::import_pose::pose_from_pdb(pose, "src/python/bindings/test/data/test_in.pdb");

	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function_legacy( scoring::PRE_TALARIS_2013_STANDARD_WTS );

	T("Score:") << scorefxn->score(pose)  << std::endl;

	protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover() );
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap() );
	//min_mover->movemap(movemap);
	min_mover->apply(pose);

	T("Scoring done!") << "---------------------" << std::endl;

	return 0;
}
