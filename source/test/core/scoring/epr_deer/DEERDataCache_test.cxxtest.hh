// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/DEERDataCache_test.cxxtest.hh
/// @brief  Tests for DEERDataCache and DEERIO class
/// @author Diego del Alamo

#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/epr_deer/DEERDataCache.hh>
#include <core/types.hh>

#include <map>


class DEERDataCacheTest : public CxxTest::TestSuite {
public:

	void setUp() {
		std::string data_file = "-epr_deer:input_files core/scoring/epr_deer/input_file_epr_deer.txt";
		std::string coords_file = "-epr_deer:coords_files core/scoring/epr_deer/coords_file.txt";
		core_init_with_additional_options( data_file + " " + coords_file );
	}

	void tearDown() {
		// NULL
	}

	void test_sorted_data_direct() {
		core::pose::Pose test_pose;
		core::io::pdb::build_pose_from_pdb_as_is( test_pose, "core/scoring/epr_deer/2lzm.pdb" );

		core::scoring::epr_deer::DEERDataCache datacache;

		// TESTING DECAY METHOD

		core::scoring::epr_deer::DEERDataOP base_op_1 = core::scoring::epr_deer::DEERDataOP(
			new core::scoring::epr_deer::DEERDecayData() );
		core::scoring::epr_deer::DEERDecayDataOP derv_op_1 = utility::pointer::dynamic_pointer_cast< core::scoring::epr_deer::DEERDecayData >( base_op_1 );
		derv_op_1->append_trace_data_and_calculate( 0.000, 1.001561, 200 );
		derv_op_1->append_trace_data_and_calculate( 0.008, 0.999551, 200 );
		derv_op_1->append_trace_data_and_calculate( 0.016, 0.995977, 200 );
		derv_op_1->append_trace_data_and_calculate( 0.024, 0.991484, 200 );
		derv_op_1->append_trace_data_and_calculate( 0.032, 0.985923, 200 );
		derv_op_1->append_trace_data_and_calculate( 0.040, 0.977163, 200 );
		derv_op_1->append_trace_data_and_calculate( 0.048, 0.969097, 200 );
		derv_op_1->append_trace_data_and_calculate( 0.056, 0.962278, 200 );
		derv_op_1->append_trace_data_and_calculate( 0.064, 0.954237, 200 );
		derv_op_1->append_trace_data_and_calculate( 0.072, 0.946555, 200 );
		derv_op_1->append_trace_data_and_calculate( 0.080, 0.93929, 200 );
		derv_op_1->append_trace_data_and_calculate( 0.088, 0.933728, 200 );
		derv_op_1->append_trace_data_and_calculate( 0.096, 0.928721, 200 );
		derv_op_1->append_trace_data_and_calculate( 0.104, 0.923318, 200 );
		derv_op_1->append_trace_data_and_calculate( 0.112, 0.920038, 200 );

		datacache.append( derv_op_1 );
		// Spot check to make sure the object is valid
		TS_ASSERT_EQUALS( datacache.size(), 1 );

		core::scoring::epr_deer::DEERDecayDataOP test_1 = utility::pointer::dynamic_pointer_cast< core::scoring::epr_deer::DEERDecayData >( datacache.at( 1 ) );

		TS_ASSERT_EQUALS( test_1->trace().size(), 15 );
		TS_ASSERT_EQUALS( test_1->k_fit(), 0.0 );
		TS_ASSERT_EQUALS( test_1->noise(), 1.0 );

		// TESTING DISTANCE METHOD

		core::scoring::epr_deer::DEERDataOP base_op_2 = core::scoring::epr_deer::DEERDataOP(
			new core::scoring::epr_deer::DEERDistanceDistribution() );
		core::scoring::epr_deer::DEERDistanceDistributionOP derv_op_2 = utility::pointer::dynamic_pointer_cast< core::scoring::epr_deer::DEERDistanceDistribution >( base_op_2 );

		std::map< core::Size, core::Real > lower_test = { {40 , 0.014}, {41 , 0.017}, {42 , 0.019}, {43 , 0.022}, {44 , 0.026}, {45 , 0.032} };
		std::map< core::Size, core::Real > mid_test = { {40 , 0.016}, {41 , 0.020}, {42 , 0.025}, {43 , 0.031}, {44 , 0.040}, {45 , 0.045} };
		std::map< core::Size, core::Real > upper_test = { {40 , 0.031}, {41 , 0.035}, {42 , 0.038}, {43 , 0.041}, {44 , 0.051}, {45 , 0.063} };

		// Initialize fictitious data
		derv_op_2->lower_bound( lower_test );
		derv_op_2->best_fit( mid_test );
		derv_op_2->upper_bound( upper_test );

		datacache.append( derv_op_2 );

		TS_ASSERT_EQUALS( datacache.size(), 2 );

		// Spot check
		core::scoring::epr_deer::DEERDistanceDistributionOP test_2 = utility::pointer::dynamic_pointer_cast< core::scoring::epr_deer::DEERDistanceDistribution >( datacache.at( 2 ) );
		TS_ASSERT_EQUALS( test_2->best_fit().size(), 6 );
		TS_ASSERT_EQUALS( test_2->lower_bound().size(), 6 );
		TS_ASSERT_EQUALS( test_2->upper_bound().size(), 6 );

		test_2->append( 46, 0.04, 0.07 );
		TS_ASSERT_EQUALS( test_2->lower_bound().size(), 7 );
		TS_ASSERT_EQUALS( test_2->upper_bound().size(), 7 );

		// TESTING SINGLE DISTANCE METHOD

		core::scoring::epr_deer::DEERDataOP base_op_3 = core::scoring::epr_deer::DEERDataOP(
			new core::scoring::epr_deer::DEERDistanceBounds() );
		core::scoring::epr_deer::DEERDistanceBoundsOP derv_op_3 = utility::pointer::dynamic_pointer_cast< core::scoring::epr_deer::DEERDistanceBounds >( base_op_3 );

		derv_op_3->bounds( 50.0, 55.0 );
		derv_op_3->step( 1.5 );

		datacache.append( derv_op_3 );
		TS_ASSERT_EQUALS( datacache.size(), 3 );

		core::scoring::epr_deer::DEERDistanceBoundsOP test_3 = utility::pointer::dynamic_pointer_cast< core::scoring::epr_deer::DEERDistanceBounds >( datacache.at( 3 ) );
		TS_ASSERT_DELTA( test_3->bounds().first,  50.0, 0.00001 );
		TS_ASSERT_DELTA( test_3->bounds().second, 55.0, 0.00001 );
		TS_ASSERT_DELTA( test_3->step(), 1.5, 0.00001 );
	}

	void test_sorted_data_from_options() {
		core::pose::Pose test_pose;
		core::io::pdb::build_pose_from_pdb_as_is( test_pose, "core/scoring/epr_deer/2lzm.pdb" );

		core::scoring::epr_deer::DEERDataCache datacache;
		datacache.fetch_and_organize_data( test_pose );

		TS_ASSERT_EQUALS( datacache.size(), 6 );
		TS_ASSERT_EQUALS( datacache.edge( 60, 94 ).size(), 5 );
		TS_ASSERT_EQUALS( datacache.edge( 94, 123 ).size(), 1 );

		core::scoring::epr_deer::DEERDataOP & dataset_1 = datacache[ 1 ];
		core::scoring::epr_deer::DEERDataOP & dataset_2 = datacache[ 2 ];
		core::scoring::epr_deer::DEERDataOP & dataset_3 = datacache[ 3 ];
		core::scoring::epr_deer::DEERDataOP & dataset_4 = datacache[ 4 ];
		core::scoring::epr_deer::DEERDataOP & dataset_5 = datacache[ 5 ];
		core::scoring::epr_deer::DEERDataOP & dataset_6 = datacache[ 6 ];

		TS_ASSERT_EQUALS( dataset_1->bins_per_angstrom(), 2 );
		TS_ASSERT_DELTA( dataset_1->relative_weight(), 1.0, 0.00001 );

		core::scoring::epr_deer::DEERDecayDataOP data_1 = utility::pointer::dynamic_pointer_cast< core::scoring::epr_deer::DEERDecayData >( dataset_1 );
		TS_ASSERT_EQUALS( data_1->trace().size(), 195 );
		TS_ASSERT_EQUALS( data_1->fit_stdev(), true );
		TS_ASSERT_DELTA( data_1->trace().at( -0.112 ), 0.920038, 0.00001 );
		TS_ASSERT_DELTA( data_1->mod_depth_bounds().first,  0.1, 0.00001 );
		TS_ASSERT_DELTA( data_1->mod_depth_bounds().second, 0.5, 0.00001 );
		TS_ASSERT_DELTA( data_1->noise(), 0.00057976, 0.00001 );

		core::scoring::epr_deer::DEERDecayDataOP data_2 = utility::pointer::dynamic_pointer_cast< core::scoring::epr_deer::DEERDecayData >( dataset_2 );
		TS_ASSERT_EQUALS( data_2->trace().size(), 128 );
		TS_ASSERT_DELTA( data_2->trace().at( 0.00 ), 1.000000, 0.00001 );
		TS_ASSERT_DELTA( data_2->trace().at( 1.936 ), 0.610182, 0.00001 );
		TS_ASSERT_EQUALS( data_2->fit_stdev(), false );
		TS_ASSERT_DELTA( data_2->mod_depth_bounds().first,  0.02, 0.00001 );
		TS_ASSERT_DELTA( data_2->mod_depth_bounds().second, 0.75, 0.00001 );
		TS_ASSERT_DELTA( data_2->noise(), 0.00425869, 0.00001 );

		core::scoring::epr_deer::DEERDistanceDistributionOP data_3 = utility::pointer::dynamic_pointer_cast< core::scoring::epr_deer::DEERDistanceDistribution >( dataset_3 );
		TS_ASSERT_EQUALS( data_3->best_fit().size(), 45 );
		TS_ASSERT_DELTA( data_3->best_fit().at( 33 ), 0.001, 0.00001 ); // two bins per angstrom

		core::scoring::epr_deer::DEERDistanceDistributionOP data_4 = utility::pointer::dynamic_pointer_cast< core::scoring::epr_deer::DEERDistanceDistribution >( dataset_4 );
		TS_ASSERT_EQUALS( data_4->bins_per_angstrom(), 5 );
		// do a spot check on a gaussian distribution

		core::scoring::epr_deer::DEERDistanceBoundsOP data_5 = utility::pointer::dynamic_pointer_cast< core::scoring::epr_deer::DEERDistanceBounds >( dataset_5 );
		TS_ASSERT_DELTA( data_5->bounds().first, 24.0, 0.00001 );
		TS_ASSERT_DELTA( data_5->bounds().second, 26.0, 0.00001 );
		TS_ASSERT_DELTA( data_5->step(), 2.5, 0.00001 );

		core::scoring::epr_deer::DEERDistanceBoundsOP data_6 = utility::pointer::dynamic_pointer_cast< core::scoring::epr_deer::DEERDistanceBounds >( dataset_6 );
		TS_ASSERT_DELTA( data_6->bounds().first, 25.0, 0.00001 );
		TS_ASSERT_DELTA( data_6->bounds().second, 25.0, 0.00001 );
		TS_ASSERT_DELTA( data_6->step(), 1.0, 0.00001 );

		TS_ASSERT_EQUALS( datacache.sl_weights().size(), 2 );
		TS_ASSERT_EQUALS( datacache.labels().size(), 2 );
		TS_ASSERT_DELTA( datacache.sl_weights().at( 1 ), 0.5231429, 0.00001 );
	}

};
