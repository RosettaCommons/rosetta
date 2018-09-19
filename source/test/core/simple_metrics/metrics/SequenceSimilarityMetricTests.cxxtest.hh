// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  test/core/simple_metrics/metrics/SequenceSimilarityMetricTests.cxxtest.hh
/// @brief  Unit tests for the SequenceSimilarityMetric.
/// @author Jack Maguire, jackmaguire1444@gmail.com

// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <core/simple_metrics/metrics/SequenceSimilarityMetric.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

// Protocols Headers

// Utility, etc Headers
#include <basic/Tracer.hh>
#include <utility/pointer/memory.hh>

#include <cmath>

static basic::Tracer TR( "test.core.simple_metrics.metrics.SequenceSimilarityMetricTests.cxxtest" );

class SequenceSimilarityMetricTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init_with_additional_options("-native core/simple_metrics/metrics/testpose.pdb");
	}

	void tearDown(){

	}

	void test_simple_case(){
		TR << "Running test_simple_case()" << std::endl;
		//Make the following mutations:
		//A -> K, score = -1
		//P -> P, score = 7
		//G -> N, score = 0
		//R -> F, score = -2
		//Total score = 3
		//Total normalized score = 0.75

		std::string const seq1 = "APGR";
		std::string const seq2 = "KPNF";

		core::simple_metrics::metrics::SequenceSimilarityMetric metric;
		core::Real const score = metric.score( seq1, seq2, false );
		core::Real const norm_score = metric.score( seq1, seq2, true );

		TR << "score: " << score << std::endl;
		TR << "norm_score: " << norm_score << std::endl;

		TS_ASSERT( std::abs( score - 3.0 ) < 0.001 );
		TS_ASSERT( std::abs( norm_score - 0.75 ) < 0.001 );

		TR.flush();
	}


	void test_selector(){
		TR << "Running test_selector()" << std::endl;

		core::pose::PoseOP pose( new core::pose::Pose );
		core::import_pose::pose_from_file( * pose, "core/simple_metrics/metrics/testpose.pdb", false, core::import_pose::PDB_file );

		core::select::residue_selector::ResidueIndexSelectorOP index_selector( utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >() );
		index_selector->set_index("1-3,5");//first 3 native residues are A,P,G, 5 is T

		//A-A score is 4
		//P-P score is 7
		//G-G score is 6
		//T-T score is 5
		//Total score should be 22

		core::simple_metrics::metrics::SequenceSimilarityMetric metric;
		metric.set_residue_selector( index_selector );
		metric.set_native_pose( pose );

		metric.set_normalize( false );
		core::Real const score = metric.calculate( * pose );
		TR << "score: " << score << std::endl;
		TS_ASSERT( std::abs( score - 22.0 ) < 0.001 );

		metric.set_normalize( true );
		core::Real const norm_score = metric.calculate( * pose );
		TR << "norm_score: " << norm_score << std::endl;
		TS_ASSERT( std::abs( norm_score - 5.5 ) < 0.001 );

		TR.flush();
	}

};
