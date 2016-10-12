// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   rosetta/benchmark/ShearMover.bench.cc
///
/// @brief  Varios moves benchmark
/// @author Sergey Lyskov


#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <protocols/simple_moves/BackboneMover.hh>

#include <core/kinematics/MoveMap.hh>

#include <apps/benchmark/performance/performance_benchmark.hh>

#include <utility/vector1.hh>


class ShearMoverBenchmark : public PerformanceBenchmark
{
public:
	core::pose::PoseOP pose;
	core::kinematics::MoveMapOP movemap;
	protocols::simple_moves::ShearMover shear_mover;

	//ShearMoverBenchmark(std::string name) : PerformanceBenchmark(name), shear_mover(300., 100) {};
	ShearMoverBenchmark(std::string name) : PerformanceBenchmark(name) {};

	virtual void setUp() {
		using namespace core;
		pose = pose::PoseOP( new pose::Pose() );
		core::import_pose::pose_from_file(*pose, "test_in.pdb", core::import_pose::PDB_file);

		movemap = kinematics::MoveMapOP( new kinematics::MoveMap() );
		movemap->set_chi( true );
		movemap->set_bb( true );

		//shear_mover.temperature(300.);
		shear_mover.nmoves(100);
		shear_mover.movemap(movemap);

		shear_mover.angle_max( 'H', 2.0 );
		shear_mover.angle_max( 'E', 2.0 );
		shear_mover.angle_max( 'L', 3.0 );

		// run once to trigger Ramachandran score calculation.
		shear_mover.apply(*pose);
	};

	virtual void run(core::Real scaleFactor) {
		//protocols::simple_moves::ShearMover shear_mover(movemap, 300., 100);
		core::Size reps( (core::Size)(5*scaleFactor) ); // amw 500 to 50
		if ( reps == 0 ) { reps = 1; } // Do at least one rep, regardless of scaling factor.
		for ( core::Size i=0; i<reps; i++ ) {
			shear_mover.apply(*pose);
		}
	};

	virtual void tearDown() {
		pose.reset();
		movemap.reset();
	}
};
