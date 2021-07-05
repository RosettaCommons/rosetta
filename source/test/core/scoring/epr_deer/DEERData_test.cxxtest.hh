// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/DEERData_test.cxxtest.hh
/// @brief  Tests for DEERData base and derived classes
/// @author Diego del Alamo

#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <core/scoring/epr_deer/Simulated4PDEERTraceFactory.hh>
#include <core/scoring/epr_deer/metrics/DEERData.hh>
#include <core/scoring/epr_deer/metrics/DEERData.fwd.hh>
#include <core/scoring/epr_deer/metrics/DEERDecayData.hh>
#include <core/scoring/epr_deer/metrics/DEERDecayData.fwd.hh>
#include <core/scoring/epr_deer/metrics/DEERDistanceDistribution.hh>
#include <core/scoring/epr_deer/metrics/DEEROverlapMethod.hh>
#include <core/scoring/epr_deer/metrics/DEEROverlapMethod.fwd.hh>
#include <core/scoring/epr_deer/metrics/DEERDistanceDistribution.fwd.hh>
#include <core/scoring/epr_deer/metrics/DEERDistanceBounds.hh>
#include <core/scoring/epr_deer/metrics/DEERDistanceBounds.fwd.hh>
#include <core/scoring/epr_deer/metrics/DEERChiSqMethod.fwd.hh>
#include <core/scoring/epr_deer/metrics/DEERChiSqMethod.hh>
#include <core/scoring/epr_deer/metrics/DEERJaccardMethod.fwd.hh>
#include <core/scoring/epr_deer/metrics/DEERJaccardMethod.hh>
#include <core/scoring/epr_deer/metrics/DEERJSMethod.fwd.hh>
#include <core/scoring/epr_deer/metrics/DEERJSMethod.hh>
#include <core/scoring/epr_deer/metrics/DEERMiscMethod.fwd.hh>
#include <core/scoring/epr_deer/metrics/DEERMiscMethod.hh>
#include <core/types.hh>

#include <iostream> // AUTO IWYU For endl, basic_ostream<>::__ostream_type

class DEERDataTest : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init_with_additional_options( "-epr_deer:print_data" );
	}

	void tearDown() {
		// NULL
	}

	void test_base_fxns() {
		utility::vector1< std::pair< core::Size, std::string > > residues_in;
		residues_in.push_back( std::make_pair( 20, "DEFAULT" ) );
		residues_in.push_back( std::make_pair( 25, "DEFAULT" ) );
		residues_in.push_back( std::make_pair( 30, "DEFAULT" ) );
		core::Size bins_per_angstrom_in = 2;
		core::Real score_in = 999.99;

		core::scoring::epr_deer::metrics::DEERDataOP test_item_1(
			new core::scoring::epr_deer::metrics::DEERDistanceDistribution() );
		core::scoring::epr_deer::metrics::DEERDistanceDistributionOP derv_op_1
			= utility::pointer::dynamic_pointer_cast<
			core::scoring::epr_deer::metrics::DEERDistanceDistribution >( test_item_1 );

		test_item_1->residues( residues_in );
		test_item_1->bins_per_a( bins_per_angstrom_in );
		test_item_1->score( score_in );

		auto residues_out = test_item_1->residues();
		auto bins_per_angstrom_out = test_item_1->bins_per_a();
		auto score_out = test_item_1->score();

		TS_ASSERT_EQUALS( residues_out, residues_in );
		TS_ASSERT_EQUALS( bins_per_angstrom_out, bins_per_angstrom_in );
		TS_ASSERT_EQUALS( score_out, score_in );

		std::map< core::Size, core::Real > test_histogram;
		test_histogram[ 41 ] = 0.15;
		test_histogram[ 42 ] = 0.20;
		test_histogram[ 43 ] = 0.25;
		test_histogram[ 44 ] = 0.20;
		test_histogram[ 45 ] = 0.15;

		std::map< core::Size, core::Real > lower_test = { {40 , 0.014},
			{41 , 0.017}, {42 , 0.019}, {43 , 0.022}, {44 , 0.026}, {45 , 0.032} };
		std::map< core::Size, core::Real > mid_test = { {40 , 0.016},
			{41 , 0.020}, {42 , 0.025}, {43 , 0.031}, {44 , 0.040}, {45 , 0.045} };
		std::map< core::Size, core::Real > upper_test = { {40 , 0.031},
			{41 , 0.035}, {42 , 0.038}, {43 , 0.041}, {44 , 0.051}, {45 , 0.063} };

		derv_op_1->lower_bound( lower_test );
		derv_op_1->best_fit( mid_test );
		derv_op_1->upper_bound( upper_test );

		core::Real method_score = test_item_1->get_score( test_histogram );
		TS_ASSERT_LESS_THAN( 0.0, method_score );
		TS_ASSERT_LESS_THAN( 0.0, test_item_1->score() );
	}

	void test_cross_entropy() {
		core::scoring::epr_deer::metrics::DEERDataOP test_item_1(
			new core::scoring::epr_deer::metrics::DEERDistanceDistribution() );
		core::scoring::epr_deer::metrics::DEERDistanceDistributionOP derv_op_1
			= utility::pointer::dynamic_pointer_cast<
			core::scoring::epr_deer::metrics::DEERDistanceDistribution >( test_item_1 );
		utility::vector1< std::pair< core::Size, std::string > > residues_in;

		residues_in.push_back( std::make_pair( 20, "DEFAULT" ) );
		residues_in.push_back( std::make_pair( 25, "DEFAULT" ) );
		residues_in.push_back( std::make_pair( 30, "DEFAULT" ) );

		std::map< core::Size, core::Real > test_histogram;
		test_histogram[ 41 ] = 0.15;
		test_histogram[ 42 ] = 0.20;
		test_histogram[ 43 ] = 0.25;
		test_histogram[ 44 ] = 0.20;
		test_histogram[ 45 ] = 0.15;

		std::map< core::Size, core::Real > lower_test = { {40 , 0.014},
			{41 , 0.017}, {42 , 0.019}, {43 , 0.022}, {44 , 0.026}, {45 , 0.032} };
		std::map< core::Size, core::Real > mid_test = { {40 , 0.016},
			{41 , 0.020}, {42 , 0.025}, {43 , 0.031}, {44 , 0.040}, {45 , 0.045} };
		std::map< core::Size, core::Real > upper_test = { {40 , 0.031},
			{41 , 0.035}, {42 , 0.038}, {43 , 0.041}, {44 , 0.051}, {45 , 0.063} };

		derv_op_1->lower_bound( lower_test );
		derv_op_1->best_fit( mid_test );
		derv_op_1->upper_bound( upper_test );

		core::Real method_score = test_item_1->get_score( test_histogram, true );
		TS_ASSERT_LESS_THAN( 0.0, method_score );
		TS_ASSERT_LESS_THAN( 0.0, test_item_1->score() );

		TS_TRACE( "Turning on bounds!" );

		derv_op_1->bounds( true );

		method_score = test_item_1->get_score( test_histogram, true );
		TS_ASSERT_LESS_THAN( 0.0, method_score );
		TS_ASSERT_LESS_THAN( 0.0, test_item_1->score() );

		TS_TRACE( "Turning on reverse calculation!" );

		derv_op_1->reverse( true );

		method_score = test_item_1->get_score( test_histogram, true );
		TS_ASSERT_LESS_THAN( 0.0, method_score );
		TS_ASSERT_LESS_THAN( 0.0, test_item_1->score() );

		TS_TRACE( "Turning on backbone calculation!" );

		derv_op_1->bb( true );
		method_score = test_item_1->get_score( test_histogram, true );
		TS_ASSERT_LESS_THAN( 0.0, method_score );
		TS_ASSERT_LESS_THAN( 0.0, test_item_1->score() );

		derv_op_1->reverse( false );
		method_score = test_item_1->get_score( test_histogram, true );
		TS_ASSERT_LESS_THAN( 0.0, method_score );
		TS_ASSERT_LESS_THAN( 0.0, test_item_1->score() );
	}

	void test_methods() {
		utility::vector1< std::pair< core::Size, std::string > > residues_in;

		residues_in.push_back( std::make_pair( 20, "DEFAULT" ) );
		residues_in.push_back( std::make_pair( 25, "DEFAULT" ) );
		residues_in.push_back( std::make_pair( 30, "DEFAULT" ) );

		std::map< core::Size, core::Real > test_histogram;
		test_histogram[ 41 ] = 0.15;
		test_histogram[ 42 ] = 0.20;
		test_histogram[ 43 ] = 0.25;
		test_histogram[ 44 ] = 0.20;
		test_histogram[ 45 ] = 0.15;

		std::map< core::Size, core::Real > lower_test = { {40 , 0.14},
			{41 , 0.17}, {42 , 0.19}, {43 , 0.22}, {44 , 0.26}, {45 , 0.32} };
		std::map< core::Size, core::Real > mid_test = { {40 , 0.016},
			{41 , 0.20}, {42 , 0.25}, {43 , 0.31}, {44 , 0.40}, {45 , 0.45} };
		std::map< core::Size, core::Real > upper_test = { {40 , 0.031},
			{41 , 0.35}, {42 , 0.38}, {43 , 0.41}, {44 , 0.51}, {45 , 0.63} };

		core::Real total = 0.0;
		for ( auto const & bin_amp : mid_test ) {
			total += bin_amp.second;
		}
		for ( auto & bin_amp : lower_test ) {
			bin_amp.second /= total;
		}

		for ( auto & bin_amp : mid_test ) {
			bin_amp.second /= total;
		}

		for ( auto & bin_amp : upper_test ) {
			bin_amp.second /= total;
		}

		/// Overlap method
		core::scoring::epr_deer::metrics::DEERDataOP test_item_1(
			new core::scoring::epr_deer::metrics::DEEROverlapMethod() );
		core::scoring::epr_deer::metrics::DEEROverlapMethodOP derv_op_1
			= utility::pointer::dynamic_pointer_cast<
			core::scoring::epr_deer::metrics::DEEROverlapMethod >( test_item_1 );

		derv_op_1->lower_bound( lower_test );
		derv_op_1->best_fit( mid_test );
		derv_op_1->upper_bound( upper_test );

		core::Real method_score = test_item_1->get_score( test_histogram, true );
		TS_ASSERT_LESS_THAN( method_score, 0.0 );
		TS_ASSERT_LESS_THAN( test_item_1->score(), 0.0 );
		TS_ASSERT_LESS_THAN( derv_op_1->score(), 0.0 );

		derv_op_1->singleval( true );

		derv_op_1->integral( true );

		method_score = test_item_1->get_score( test_histogram, true );
		TS_ASSERT_LESS_THAN( 0.0, method_score );
		TS_ASSERT_LESS_THAN( 0.0, test_item_1->score() );

		derv_op_1->bounds( true );

		method_score = test_item_1->get_score( test_histogram, true );
		TS_ASSERT_LESS_THAN( 0.0, method_score );
		TS_ASSERT_LESS_THAN( 0.0, test_item_1->score() );

		// Misc method
		core::scoring::epr_deer::metrics::DEERDataOP test_item_2(
			new core::scoring::epr_deer::metrics::DEERMiscMethod() );
		core::scoring::epr_deer::metrics::DEERMiscMethodOP derv_op_2
			= utility::pointer::dynamic_pointer_cast<
			core::scoring::epr_deer::metrics::DEERMiscMethod >( test_item_2 );

		derv_op_2->lower_bound( lower_test );
		derv_op_2->best_fit( mid_test );
		derv_op_2->upper_bound( upper_test );

		method_score = test_item_2->get_score( test_histogram, true );
		TS_ASSERT_LESS_THAN( 0.0, method_score );
		TS_ASSERT_LESS_THAN( 0.0, test_item_2->score() );
		TS_ASSERT_LESS_THAN( 0.0, derv_op_2->score() );

		derv_op_2->mode( "HELLINGER" );
		method_score = test_item_2->get_score( test_histogram, true );
		TS_ASSERT_LESS_THAN( 0.0, method_score );
		TS_ASSERT_LESS_THAN( 0.0, test_item_2->score() );
		TS_ASSERT_LESS_THAN( 0.0, derv_op_2->score() );

		derv_op_2->mode( "BHATTACHARYYA" );
		method_score = test_item_2->get_score( test_histogram, true );
		TS_ASSERT_LESS_THAN( 0.0, method_score );
		TS_ASSERT_LESS_THAN( 0.0, test_item_2->score() );
		TS_ASSERT_LESS_THAN( 0.0, derv_op_2->score() );

		// JS Method
		core::scoring::epr_deer::metrics::DEERDataOP test_item_3(
			new core::scoring::epr_deer::metrics::DEERJSMethod() );
		core::scoring::epr_deer::metrics::DEERJSMethodOP derv_op_3
			= utility::pointer::dynamic_pointer_cast<
			core::scoring::epr_deer::metrics::DEERJSMethod >( test_item_3 );

		derv_op_3->lower_bound( lower_test );
		derv_op_3->best_fit( mid_test );
		derv_op_3->upper_bound( upper_test );

		method_score = test_item_3->get_score( test_histogram, true );
		TS_ASSERT_LESS_THAN( 0.0, method_score );
		TS_ASSERT_LESS_THAN( 0.0, test_item_3->score() );
		TS_ASSERT_LESS_THAN( 0.0, derv_op_3->score() );

		// Jaccard Method
		core::scoring::epr_deer::metrics::DEERDataOP test_item_4(
			new core::scoring::epr_deer::metrics::DEERJaccardMethod() );
		core::scoring::epr_deer::metrics::DEERJaccardMethodOP derv_op_4
			= utility::pointer::dynamic_pointer_cast<
			core::scoring::epr_deer::metrics::DEERJaccardMethod >( test_item_4 );

		derv_op_4->lower_bound( lower_test );
		derv_op_4->best_fit( mid_test );
		derv_op_4->upper_bound( upper_test );

		method_score = test_item_4->get_score( test_histogram, true );
		TS_ASSERT_LESS_THAN( method_score, 0.0 );
		TS_ASSERT_LESS_THAN( test_item_4->score(), 0.0 );
		TS_ASSERT_LESS_THAN( derv_op_4->score(), 0.0 );

		// Chi squared Method
		core::scoring::epr_deer::metrics::DEERDataOP test_item_5(
			new core::scoring::epr_deer::metrics::DEERChiSqMethod() );
		core::scoring::epr_deer::metrics::DEERChiSqMethodOP derv_op_5
			= utility::pointer::dynamic_pointer_cast<
			core::scoring::epr_deer::metrics::DEERChiSqMethod >( test_item_5 );

		derv_op_5->lower_bound( lower_test );
		derv_op_5->best_fit( mid_test );
		derv_op_5->upper_bound( upper_test );

		method_score = test_item_5->get_score( test_histogram, true );
		TS_ASSERT_LESS_THAN( 0.0, method_score );
		TS_ASSERT_LESS_THAN( 0.0, test_item_5->score() );
		TS_ASSERT_LESS_THAN( 0.0, derv_op_5->score() );

		derv_op_5->reverse( true );

		method_score = test_item_5->get_score( test_histogram, true );
		TS_ASSERT_LESS_THAN( 0.0, method_score );
		TS_ASSERT_LESS_THAN( 0.0, test_item_5->score() );
		TS_ASSERT_LESS_THAN( 0.0, derv_op_5->score() );

	}

	void test_derived_bounds_fxns() {
		core::Size bins_per_angstrom_in = 1;
		core::Real lower_bound_in = 19.0;
		core::Real upper_bound_in = 21.0;
		core::Real step_in = 1.0;

		core::scoring::epr_deer::metrics::DEERDataOP data_obj(
			core::scoring::epr_deer::metrics::DEERDataOP(
			new core::scoring::epr_deer::metrics::DEERDistanceBounds() ) );
		core::scoring::epr_deer::metrics::DEERDistanceBoundsOP test_item_1
			= utility::pointer::dynamic_pointer_cast<
			core::scoring::epr_deer::metrics::DEERDistanceBounds >( data_obj );

		test_item_1->bins_per_a( bins_per_angstrom_in );
		test_item_1->score( 999.99 );
		test_item_1->bounds( lower_bound_in, upper_bound_in );
		test_item_1->step( step_in );

		std::pair< core::Real, core::Real > bounds_out = test_item_1->bounds();
		core::Real step_out = test_item_1->step();

		TS_ASSERT_EQUALS( bounds_out.first, lower_bound_in );
		TS_ASSERT_EQUALS( bounds_out.second, upper_bound_in );
		TS_ASSERT_EQUALS( step_out, step_in );

		std::map< core::Size, core::Real > test_histogram_1;
		test_histogram_1[ 20 ] = 1.0;
		data_obj->get_score( test_histogram_1, true );

		TS_ASSERT( data_obj->score() == 0.0 );

		std::map< core::Size, core::Real > test_histogram_2;
		test_histogram_2[ 22 ] = 1.0;

		data_obj->get_score( test_histogram_2, true );

		TS_ASSERT( data_obj->score() > 0.0 );

		std::cout << data_obj->score() << std::endl;

	}

	void test_derived_distance_fxns() {
		core::Size bins_per_angstrom_in = 1;
		std::map< core::Size, core::Real > exp_map_lower;
		std::map< core::Size, core::Real > exp_map_best;
		std::map< core::Size, core::Real > exp_map_upper;

		exp_map_lower[ 20 ] = 0.125;
		exp_map_lower[ 21 ] = 0.25;
		exp_map_lower[ 22 ] = 0.125;

		exp_map_best[ 20 ] = 0.25;
		exp_map_best[ 21 ] = 0.5;
		exp_map_best[ 22 ] = 0.25;

		exp_map_upper[ 20 ] = 0.5;
		exp_map_upper[ 21 ] = 1.0;
		exp_map_upper[ 22 ] = 0.5;

		core::scoring::epr_deer::metrics::DEERDataOP data_obj(
			new core::scoring::epr_deer::metrics::DEERDistanceDistribution() );
		core::scoring::epr_deer::metrics::DEERDistanceDistributionOP test_item_1
			= utility::pointer::dynamic_pointer_cast<
			core::scoring::epr_deer::metrics::DEERDistanceDistribution >( data_obj );

		test_item_1->bins_per_a( bins_per_angstrom_in );
		test_item_1->best_fit( exp_map_best );
		test_item_1->lower_bound( exp_map_lower );
		test_item_1->upper_bound( exp_map_upper );

		TS_ASSERT_EQUALS( exp_map_best, test_item_1->best_fit() );
		TS_ASSERT_EQUALS( exp_map_lower, test_item_1->lower_bound() );
		TS_ASSERT_EQUALS( exp_map_upper, test_item_1->upper_bound() );

		std::map< core::Size, core::Real > test_distribution = exp_map_best;

		core::Real expected_score = 0.0;
		for ( auto const & item : test_distribution ) {
			expected_score += item.second * log( 1.0 / item.second );
		}
		TS_ASSERT_DELTA( test_item_1->get_score( test_distribution ),
			expected_score, 1e-6 );
	}

	void test_derived_decay_fxns() {
		// First check all parameters can be accessed accurately
		core::scoring::epr_deer::metrics::DEERDataOP data_obj(
			new core::scoring::epr_deer::metrics::DEERDecayData() );
		core::scoring::epr_deer::metrics::DEERDecayDataOP test_item_1
			= utility::pointer::dynamic_pointer_cast<
			core::scoring::epr_deer::metrics::DEERDecayData >( data_obj );
		data_obj->bins_per_a( 1.0 );
		utility::vector1< core::Real > test_trace, test_time_pts;
		for ( core::Size i = 0; i <= 250; ++i ) {
			core::Real time_pt = 0.0 + ( 0.008 * i );
			core::Real signal = 1.0 - ( 0.001 * i );
			test_trace.push_back( signal );
			test_time_pts.push_back( time_pt );
		}
		for ( core::Size i = 251; i <= 500; ++i ) {
			core::Real time_pt = 0.0 + ( 0.008 * i );
			core::Real signal = 0.75;
			test_trace.push_back( signal );
			test_time_pts.push_back( time_pt );
		}

		test_item_1->init_factory( test_trace, test_time_pts );

		TS_ASSERT_EQUALS( test_trace.size(), test_item_1->factory().trace().size() );
		TS_ASSERT_EQUALS( test_trace, test_item_1->factory().trace() );
		TS_ASSERT_EQUALS( test_time_pts.size(), test_item_1->factory().time_pts().size() );
		TS_ASSERT_EQUALS( test_time_pts, test_item_1->factory().time_pts() );
		TS_ASSERT_EQUALS( "3D", test_item_1->factory().bckg_type() );

		test_item_1->factory().bckg_type( "NON-3D" );

		core::Real test_noise = 0.02;

		test_item_1->noise( test_noise );

		TS_ASSERT_EQUALS( test_noise, test_item_1->noise() );

		std::map< core::Size, core::Real > test_histogram;
		for ( core::Size dist = 30; dist < 35; ++dist ) {
			test_histogram[ dist ] = 0.2;
		}
		auto sim_trace = test_item_1->factory().trace_from_distr(
			test_histogram );

		TS_ASSERT( data_obj->get_score( test_histogram, true ) > 0.0 );
	}
};
