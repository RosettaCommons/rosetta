// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   rosetta/benchmark/Design.bench.hh
///
/// @brief Perform a complete redesign of 1A30 (HIV protease/inhibitor).
///  (99 res dimer bound to 3 res. peptide= 201 residue)
/// takes about 1 minute on my machine
/// @author Gordon Lemmon

#include <protocols/simple_moves/PackRotamersMover.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <apps/benchmark/performance/performance_benchmark.hh>

//Auto Headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>

class DesignBenchmark : public PerformanceBenchmark
{
public:
	DesignBenchmark(std::string name) : PerformanceBenchmark(name) {};

	protocols::simple_moves::PackRotamersMover pack_mover;
	core::pose::Pose design_pose;

	virtual void setUp() {
		//core_init_with_additional_options( "-ex1" );// I can't get this to work
		core::import_pose::pose_from_pdb(design_pose, "design_in.pdb");
	};

	virtual void run(core::Real scaleFactor) {
		int reps( 1 * scaleFactor );
		if( reps == 0 ) { reps = 1; } // Do at least one repetition, regardless of scaling factor
		for(int i=0; i< reps; i++) {
			pack_mover.apply(design_pose);
		}
		//design_pose.dump_scored_pdb("design_out.pdb", *pack_mover.scorefxn());// write out for debug
	};

	virtual void tearDown() {};
};
