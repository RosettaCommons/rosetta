// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/simple_filters/TotalSasaFilter.cxxtest.hh
/// @brief
/// @author Rocco Moretti (rmoretti@u.washington.edu)
///     (Based on code by Ron Jacak)

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <protocols/simple_filters/TotalSasaFilter.hh>

#include <core/scoring/ScoreFunction.hh>

#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>

// Package Headers
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

//Auto Headers
#include <core/pose/util.hh>
#include <utility/vector1.hh>

// --------------- Test Class --------------- //

// using declarations
using namespace core;

class TotalSasaFilterTests : public CxxTest::TestSuite {

public:

	Real TOLERATED_ERROR;
	pose::Pose pose;
	scoring::ScoreFunctionOP sf;

	// --------------- Fixtures --------------- //

	// Shared initialization goes here.
	void setUp() {
		core_init();
		TOLERATED_ERROR = 0.1;

		using namespace core::pose::metrics;
		if ( !CalculatorFactory::Instance().check_calculator_exists( "sasa" ) ) {
			PoseMetricCalculatorOP sasa_calculator( new core::pose::metrics::simple_calculators::SasaCalculatorLegacy() );
			CalculatorFactory::Instance().register_calculator( "sasa", sasa_calculator );
		}

	}

	void tearDown() {}

	// --------------- Test Cases --------------- //

	void test_calc_total_sasa_glycine() {

		core::import_pose::pose_from_file( pose, "core/scoring/nonideal_glycine.pdb" , core::import_pose::PDB_file);

		for ( Size ii=1; ii <= pose.size(); ii+=3 ) {
			pose.set_phi( ii, -150.0 );
			pose.set_psi( ii, 150.0 );
			pose.set_omega( ii, 180.0 );
		}

		protocols::simple_filters::TotalSasaFilter test;

		// Note 223.3912 is different from what Ron got in his tests (223.8206) - assuming small differences are due to different settings, etc.
		TS_ASSERT_DELTA( test.report_sm(pose), 221.0551, TOLERATED_ERROR );

		core::pack::task::TaskFactoryOP factory( new core::pack::task::TaskFactory ); //Single residue, will repack just the single residue
		test.task_factory(factory);

		TS_ASSERT_DELTA( test.report_sm(pose), 221.0551, TOLERATED_ERROR );
	}

	void test_calc_total_sasa_smallprotein() {

		pose = create_1ten_pdb_pose();

		protocols::simple_filters::TotalSasaFilter test;

		//TS_ASSERT_DELTA( test.report_sm(pose), 5236.9648, TOLERATED_ERROR );
		//TS_ASSERT_DELTA( test.report_sm(pose), 5135.0014, TOLERATED_ERROR );
		TS_ASSERT_DELTA( test.report_sm(pose), 5112.5681, TOLERATED_ERROR );

		core::pack::task::TaskFactoryOP factory( new core::pack::task::TaskFactory );
		core::pack::task::operation::PreventRepackingOP prt( new core::pack::task::operation::PreventRepacking );
		for ( core::Size ii(1); ii <= pose.size(); ++ii ) {
			if ( ii != 1 && ii != 11 && ii != 33 ) {
				prt->include_residue(ii);
			}
		}
		factory->push_back(prt);
		test.task_factory(factory);

		TS_ASSERT_DELTA( test.report_sm(pose), 166.2781, TOLERATED_ERROR );
		//TS_ASSERT_DELTA( test.report_sm(pose), 180.6644, TOLERATED_ERROR );
		//TS_ASSERT_DELTA( rsd_sasa[1], 63.4461, TOLERATED_ERROR );
		//TS_ASSERT_DELTA( rsd_sasa[11], 117.8897, TOLERATED_ERROR );
		//TS_ASSERT_DELTA( rsd_sasa[33], 0.0, TOLERATED_ERROR );

	}

	void test_calc_per_res_hydrophobic_sasa() {

		pose = create_1ten_pdb_pose();
		pose.dump_pdb( "1ten_from_total_sasa_filter.pdb" );
		//std::cout << "Nresidues: " << pose.size() << std::endl;

		protocols::simple_filters::TotalSasaFilter test(0, /*hydrophobic=*/ true, /*polar=*/ false );

		//Note 3302.5455 is different from what Ron got in his tests (3212.4579) - assuming small differences are due to slighlty different settings, etc.
		//TS_ASSERT_DELTA( test.report_sm(pose), 3311.7659, TOLERATED_ERROR );
		//TS_ASSERT_DELTA( test.report_sm(pose), 2772.0655, TOLERATED_ERROR );
		TS_ASSERT_DELTA( test.report_sm(pose), 2767.7070, TOLERATED_ERROR );

		core::pack::task::TaskFactoryOP factory( new core::pack::task::TaskFactory );
		core::pack::task::operation::PreventRepackingOP prt( new core::pack::task::operation::PreventRepacking );
		for ( core::Size ii(1); ii <= pose.size(); ++ii ) {
			if ( ii != 1 && ii != 11 && ii != 33 ) {
				prt->include_residue(ii);
			}
		}
		factory->push_back(prt);
		test.task_factory(factory);

		// Note: 95.2054 is different from what Ron got in his tests (90.1560) - assuming difference is due to differences in settings.
		//TS_ASSERT_DELTA( test.report_sm(pose), 95.7852, TOLERATED_ERROR );
		TS_ASSERT_DELTA( test.report_sm(pose), 55.0774, TOLERATED_ERROR );
	}

};


