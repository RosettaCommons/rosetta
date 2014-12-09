// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   rosetta/benchmark/Minimizer.bench.cc
///
/// @brief  Varios moves benchmark
/// @author Sergey Lyskov


#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <apps/benchmark/performance/performance_benchmark.hh>

#include <utility/vector1.hh>

using namespace core;

enum  ScoreFnType {SFT_dfpmin, SFT_dfpmin_armijo, SFT_dfpmin_armijo_nonmonotone};

template  <ScoreFnType sft, int TScale>
class MinimizerBenchmark : public PerformanceBenchmark
{
public:
	pose::PoseOP start_pose;
	kinematics::MoveMap mm;
	core::scoring::ScoreFunctionOP scorefxn;
	core::optimization::AtomTreeMinimizer minimizer;

	MinimizerBenchmark(std::string name) : PerformanceBenchmark(name) {};

	virtual void setUp() {
		start_pose = pose::PoseOP( new pose::Pose() );
		core::import_pose::pose_from_pdb(*start_pose, "test_in.pdb");

		scorefxn = core::scoring::get_score_function();

				//kinematics::MoveMap mm;
				for ( int i=30; i<= 35; ++i ) {
						mm.set_bb ( i, true );
						mm.set_chi( i, true );
				}

		(*scorefxn)( *start_pose ); // to triger dunbrack loading/calcualtion
	};

	virtual void run(core::Real scaleFactor) {
		core::Size reps( (core::Size)(TScale*scaleFactor) );
		if( reps == 0 ) { reps = 1; } // do at least one rep, regardless of scale factor
		for(core::Size i=0; i<reps; i++) {
			std::string stype = "unknow";
			if( sft == SFT_dfpmin ) stype = "dfpmin";
			if( sft == SFT_dfpmin_armijo ) stype = "dfpmin_armijo";
			if( sft == SFT_dfpmin_armijo_nonmonotone ) stype = "dfpmin_armijo_nonmonotone";
			core::optimization::MinimizerOptions options( stype/*"dfpmin"*/, 1e-1, true, true );
			core::pose::Pose pose;
			pose = *start_pose;
			minimizer.run( pose, mm, *scorefxn, options );
		}
	};

	virtual void tearDown() {};
};

typedef MinimizerBenchmark<SFT_dfpmin, 12> MinimizerBenchmark_dfpmin;
typedef MinimizerBenchmark<SFT_dfpmin_armijo, 1> MinimizerBenchmark_dfpmin_armijo;
typedef MinimizerBenchmark<SFT_dfpmin_armijo_nonmonotone, 1> MinimizerBenchmark_dfpmin_armijo_nonmonotone;


//class MinimizerBenchmark_dfpmin : public MinimizerBenchmark

/*
				{ // armijo
						MinimizerOptions options( "dfpmin_armijo", 1e-1, true, true );
						Pose pose;
						pose = start_pose;
						minimizer.run( pose, mm, *scorefxn, options );
						core::Real score = (*scorefxn)( pose );
						TR << "dfpmin_armijo/standard: " << score << "\n";
						TS_ASSERT_DELTA(score, 385.767, err_tol);
				}

				{ // non-monotone
						MinimizerOptions options( "dfpmin_armijo_nonmonotone", 1e-1, true, true );
						Pose pose;
						pose = start_pose;
						minimizer.run( pose, mm, *scorefxn, options );
						core::Real score = (*scorefxn)( pose );
						TR << "dfpmin_armijo_nonmonotone/standard: " << score << "\n";
						TS_ASSERT_DELTA(score, 385.945, err_tol);
				} */

