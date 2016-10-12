// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   rosetta/benchmark/score.bench.cc
///
/// @brief  Scoring benchmark
/// @author Sergey Lyskov

#include <apps/benchmark/performance/performance_benchmark.hh>

#include <core/pose/Pose.hh>
//#include <core/import_pose/import_pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>

#include <utility/vector1.hh>

class ScoreBenchmark : public PerformanceBenchmark
{
public:
	core::pose::Pose pose;

	ScoreBenchmark(std::string name) : PerformanceBenchmark(name) {};

	virtual void setUp() {
		core::import_pose::pose_from_file(pose, "test_in.pdb", core::import_pose::PDB_file);
	}

	virtual void run(core::Real scaleFactor) {
		core::scoring::ScoreFunction scorefxn;
		for ( int i=0; i<150*scaleFactor; i++ ) { // /100
			scorefxn(pose);
			pose.energies().clear();
		}
	};

	virtual void tearDown() {};
};
