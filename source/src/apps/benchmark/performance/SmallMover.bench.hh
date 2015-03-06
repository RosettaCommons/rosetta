// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   rosetta/benchmark/SmallMover.bench.cc
///
/// @brief  Varios moves benchmark
/// @author Sergey Lyskov


#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <protocols/simple_moves/BackboneMover.hh>

#include <core/kinematics/MoveMap.hh>

#include <apps/benchmark/performance/performance_benchmark.hh>

#include <utility/vector1.hh>


using namespace core;

class SmallMoverBenchmark : public PerformanceBenchmark
{
public:
	pose::PoseOP pose;
	kinematics::MoveMapOP movemap;
	protocols::simple_moves::SmallMover small_mover;

	SmallMoverBenchmark(std::string name) : PerformanceBenchmark(name) {};

	virtual void setUp() {
		pose = pose::PoseOP( new pose::Pose() );
		core::import_pose::pose_from_pdb(*pose, "test_in.pdb");

		movemap = kinematics::MoveMapOP( new kinematics::MoveMap() );
		movemap->set_chi( true );
		movemap->set_bb( true );

		//small_mover.temperature(300.);
		small_mover.nmoves(100);
		small_mover.movemap(movemap);

		small_mover.angle_max( 'H', 2.0 );
		small_mover.angle_max( 'E', 2.0 );
		small_mover.angle_max( 'L', 3.0 );

		// run once to trigger Ramachandran score calculation.
		small_mover.apply(*pose);
	}

	virtual void run(core::Real scaleFactor) {
		core::Size reps( (core::Size)(50*scaleFactor) ); // amw 500 to 50
		if( reps == 0 ) { reps = 1; } // do at least one rep, regardless of scale factor
		for(core::Size i=0; i<reps; i++) {
			small_mover.apply(*pose);
		}
	}

	virtual void tearDown() {
		pose.reset();
		movemap.reset();
	}
};
