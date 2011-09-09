// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
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

#include <protocols/moves/BackboneMover.hh>

#include <core/kinematics/MoveMap.hh>

#include "benchmark.hh"

using namespace core;

class SmallMoverBenchmark : public Benchmark
{
public:
	pose::PoseOP pose;
	kinematics::MoveMapOP movemap;
	protocols::moves::SmallMover small_mover;

	SmallMoverBenchmark(std::string name) : Benchmark(name) {};

	virtual void setUp() {
		pose = new pose::Pose();
		core::import_pose::pose_from_pdb(*pose, "test_in.pdb");

		movemap = new kinematics::MoveMap();
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
	};

	virtual void run(core::Real scaleFactor) {
		for(int i=0; i<2000*scaleFactor; i++) {
			small_mover.apply(*pose);
		}
	};

	virtual void tearDown() {};
};
