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
#include <core/scoring/epr_deer/metrics/DEERData.hh>
#include <core/scoring/epr_deer/metrics/DEERData.fwd.hh>
#include <core/scoring/epr_deer/metrics/DEERDecayData.hh>
#include <core/scoring/epr_deer/metrics/DEERDecayData.fwd.hh>
#include <core/scoring/epr_deer/metrics/DEERDistanceDistribution.hh>
#include <core/scoring/epr_deer/metrics/DEERDistanceDistribution.fwd.hh>
#include <core/scoring/epr_deer/metrics/DEERDistanceBounds.hh>
#include <core/scoring/epr_deer/metrics/DEERDistanceBounds.fwd.hh>
#include <core/types.hh>

#include <map>

#include <iostream> // AUTO IWYU For operator<<, basic_ostream, basic_os...


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

		core::scoring::epr_deer::metrics::DEERDataOP base_op_1(
			new core::scoring::epr_deer::metrics::DEERDecayData() );
		core::scoring::epr_deer::metrics::DEERDecayDataOP derv_op_1
			= utility::pointer::dynamic_pointer_cast<
			core::scoring::epr_deer::metrics::DEERDecayData >( base_op_1 );
		derv_op_1->fit_stdev( true );
		derv_op_1->bckg( "3D" );
		derv_op_1->max_dist( 100 );

		utility::vector1< core::Real > exp_data{ 1.001561, 0.000551, 0.995977,
			0.991484, -.985823, 0.977163, 0.969097, 0.962278, 0.953237, 0.946555,
			0.939290, 0.933728, 0.928721, 0.923318, 0.920038 };
		utility::vector1< core::Real > time_pts;
		for ( core::Size i = 0; i < exp_data.size(); ++i ) {
			time_pts.push_back( 0.008 * i );
		}

		derv_op_1->init_factory( exp_data, time_pts );
		datacache.append( derv_op_1 );
		// Spot check to make sure the object is valid
		TS_ASSERT_EQUALS( datacache.size(), 1 );

		core::scoring::epr_deer::metrics::DEERDecayDataOP test_1
			= utility::pointer::dynamic_pointer_cast<
			core::scoring::epr_deer::metrics::DEERDecayData >( datacache.at( 1 ) );

		TS_ASSERT_EQUALS( test_1->factory().trace().size(), 15 );
		TS_ASSERT_EQUALS( test_1->noise(), 1.0 );

		// TESTING DISTANCE METHOD

		core::scoring::epr_deer::metrics::DEERDataOP base_op_2(
			new core::scoring::epr_deer::metrics::DEERDistanceDistribution() );
		core::scoring::epr_deer::metrics::DEERDistanceDistributionOP derv_op_2
			= utility::pointer::dynamic_pointer_cast<
			core::scoring::epr_deer::metrics::DEERDistanceDistribution >( base_op_2 );

		std::map< core::Size, core::Real > lower_test = { {40 , 0.014},
			{41 , 0.017}, {42 , 0.019}, {43 , 0.022}, {44 , 0.026}, {45 , 0.032} };
		std::map< core::Size, core::Real > mid_test = { {40 , 0.016},
			{41 , 0.020}, {42 , 0.025}, {43 , 0.031}, {44 , 0.040}, {45 , 0.045} };
		std::map< core::Size, core::Real > upper_test = { {40 , 0.031},
			{41 , 0.035}, {42 , 0.038}, {43 , 0.041}, {44 , 0.051}, {45 , 0.063} };

		// Initialize fictitious data
		derv_op_2->lower_bound( lower_test );
		derv_op_2->best_fit( mid_test );
		derv_op_2->upper_bound( upper_test );

		datacache.append( derv_op_2 );

		TS_ASSERT_EQUALS( datacache.size(), 2 );

		// Spot check
		core::scoring::epr_deer::metrics::DEERDistanceDistributionOP test_2
			= utility::pointer::dynamic_pointer_cast<
			core::scoring::epr_deer::metrics::DEERDistanceDistribution >( datacache.at( 2 ) );
		TS_ASSERT_EQUALS( test_2->best_fit().size(), 6 );
		TS_ASSERT_EQUALS( test_2->lower_bound().size(), 6 );
		TS_ASSERT_EQUALS( test_2->upper_bound().size(), 6 );

		// TESTING SINGLE DISTANCE METHOD

		core::scoring::epr_deer::metrics::DEERDataOP base_op_3
			= core::scoring::epr_deer::metrics::DEERDataOP(
			new core::scoring::epr_deer::metrics::DEERDistanceBounds() );
		core::scoring::epr_deer::metrics::DEERDistanceBoundsOP derv_op_3
			= utility::pointer::dynamic_pointer_cast<
			core::scoring::epr_deer::metrics::DEERDistanceBounds >( base_op_3 );

		derv_op_3->bounds( 50.0, 55.0 );
		derv_op_3->step( 1.5 );

		datacache.append( derv_op_3 );
		TS_ASSERT_EQUALS( datacache.size(), 3 );

		core::scoring::epr_deer::metrics::DEERDistanceBoundsOP test_3
			= utility::pointer::dynamic_pointer_cast<
			core::scoring::epr_deer::metrics::DEERDistanceBounds >( datacache.at( 3 ) );
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
		TS_ASSERT_EQUALS( datacache.edge( 60, 94 ).size(), 4 );
		TS_ASSERT_EQUALS( datacache.edge( 94, 123 ).size(), 1 );

		core::scoring::epr_deer::metrics::DEERDataOP & dataset_1 = datacache[ 1 ];
		core::scoring::epr_deer::metrics::DEERDataOP & dataset_2 = datacache[ 2 ];
		core::scoring::epr_deer::metrics::DEERDataOP & dataset_3 = datacache[ 3 ];
		core::scoring::epr_deer::metrics::DEERDataOP & dataset_4 = datacache[ 4 ];
		core::scoring::epr_deer::metrics::DEERDataOP & dataset_5 = datacache[ 5 ];
		core::scoring::epr_deer::metrics::DEERDataOP & dataset_6 = datacache[ 6 ];

		TS_ASSERT_EQUALS( dataset_1->bins_per_a(), 2 );

		core::scoring::epr_deer::metrics::DEERDecayDataOP data_1
			= utility::pointer::dynamic_pointer_cast<
			core::scoring::epr_deer::metrics::DEERDecayData >( dataset_1 );

		TS_ASSERT_EQUALS( data_1->factory().trace().size(), 196 );
		TS_ASSERT_EQUALS( data_1->fit_stdev(), true );
		TS_ASSERT_DELTA( data_1->factory().trace().at( 2 ), 0.920038, 1e-6 );
		TS_ASSERT_DELTA( data_1->noise(), 0.00057976, 1e-6 );

		core::scoring::epr_deer::metrics::DEERDecayDataOP data_2
			= utility::pointer::dynamic_pointer_cast<
			core::scoring::epr_deer::metrics::DEERDecayData >( dataset_2 );

		TS_ASSERT_EQUALS( data_2->factory().trace().size(), 128 );
		TS_ASSERT_DELTA( data_2->factory().trace().at( 1 ), 1.000000, 0.00001 );
		TS_ASSERT_DELTA( data_2->factory().trace().at( 128 ), 0.610182, 0.00001 );
		TS_ASSERT_EQUALS( data_2->fit_stdev(), false );
		TS_ASSERT_DELTA( data_2->noise(), 0.00425869, 0.00001 );

		core::scoring::epr_deer::metrics::DEERDistanceDistributionOP data_3
			= utility::pointer::dynamic_pointer_cast<
			core::scoring::epr_deer::metrics::DEERDistanceDistribution >( dataset_3 );

		for ( auto const & dist_amp : data_3->best_fit() ) {
			std::cout << "\t" << dist_amp.first << "\t" << dist_amp.second << std::endl;
		}

		TS_ASSERT_EQUALS( data_3->best_fit().size(), 45 );
		TS_ASSERT_DELTA( data_3->best_fit().at( 33 ), 0.001, 1e-6 ); // two bins per angstrom

		core::scoring::epr_deer::metrics::DEERDistanceDistributionOP data_4
			= utility::pointer::dynamic_pointer_cast<
			core::scoring::epr_deer::metrics::DEERDistanceDistribution >( dataset_4 );

		dataset_4->bins_per_a( 5 );
		TS_ASSERT_EQUALS( dataset_4->bins_per_a(), 5 );
		// do a spot check on a gaussian distribution

		core::scoring::epr_deer::metrics::DEERDistanceDistributionOP data_5
			= utility::pointer::dynamic_pointer_cast<
			core::scoring::epr_deer::metrics::DEERDistanceDistribution >( dataset_5 );

		TS_ASSERT_LESS_THAN( 0.0, data_5->best_fit().at( 80 ) );

		core::scoring::epr_deer::metrics::DEERDistanceBoundsOP data_6
			= utility::pointer::dynamic_pointer_cast<
			core::scoring::epr_deer::metrics::DEERDistanceBounds >( dataset_6 );

		TS_ASSERT_DELTA( data_6->bounds().first, 25.0, 1e-6 );
		TS_ASSERT_DELTA( data_6->bounds().second, 25.0, 1e-6 );
		TS_ASSERT_DELTA( data_6->step(), 1.0, 1e-6 );

		TS_ASSERT_EQUALS( datacache.sl_weights().size(), 2 );
		TS_ASSERT_EQUALS( datacache.labels().size(), 2 );
		TS_ASSERT_DELTA( datacache.sl_weights().at( 1 ), 0.5231429, 1e-6 );
	}

};
