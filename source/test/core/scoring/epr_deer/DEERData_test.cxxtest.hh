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
#include <core/scoring/epr_deer/DEERData.hh>
#include <core/types.hh>

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
		core::Size bins_per_angstrom_in = 5;
		core::Real score_in = 999.99;

		core::scoring::epr_deer::DEERDataOP test_item_1( core::scoring::epr_deer::DEERDataOP( new core::scoring::epr_deer::DEERData() ) );

		test_item_1->residues( residues_in );
		test_item_1->bins_per_angstrom( bins_per_angstrom_in );
		test_item_1->score( score_in );

		auto residues_out = test_item_1->residues();
		auto bins_per_angstrom_out = test_item_1->bins_per_angstrom();
		auto score_out = test_item_1->score();

		TS_ASSERT_EQUALS( residues_out, residues_in );
		TS_ASSERT_EQUALS( bins_per_angstrom_out, bins_per_angstrom_in );
		TS_ASSERT_EQUALS( score_out, score_in );

		std::map< core::Size, core::Real > test_histogram;
		test_histogram[ 20 ] = 0.5;
		test_histogram[ 21 ] = 0.5;

		core::Real method_score = test_item_1->get_score( test_histogram, true );
		TS_ASSERT_EQUALS( method_score, 0.0 );
		TS_ASSERT_EQUALS( test_item_1->score(), 0.0 );
	}

	void test_derived_bounds_fxns() {
		core::Size bins_per_angstrom_in = 1;
		core::Real lower_bound_in = 19.0;
		core::Real upper_bound_in = 21.0;
		core::Real step_in = 1.0;

		core::scoring::epr_deer::DEERDataOP data_obj( core::scoring::epr_deer::DEERDataOP( new core::scoring::epr_deer::DEERDistanceBounds() ) );
		core::scoring::epr_deer::DEERDistanceBoundsOP test_item_1 = utility::pointer::dynamic_pointer_cast< core::scoring::epr_deer::DEERDistanceBounds >( data_obj );

		test_item_1->bins_per_angstrom( bins_per_angstrom_in );
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
		test_item_1->get_score( test_histogram_1, true );

		TS_ASSERT_EQUALS( test_item_1->score(), 0.0 );

		std::map< core::Size, core::Real > test_histogram_2;
		test_histogram_2[ 22 ] = 1.0;

		test_item_1->get_score( test_histogram_2, true );

		TS_ASSERT_EQUALS( test_item_1->score(), 1.0 );

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

		core::scoring::epr_deer::DEERDataOP data_obj( core::scoring::epr_deer::DEERDataOP( new core::scoring::epr_deer::DEERDistanceDistribution() ) );
		core::scoring::epr_deer::DEERDistanceDistributionOP test_item_1 = utility::pointer::dynamic_pointer_cast< core::scoring::epr_deer::DEERDistanceDistribution >( data_obj );

		test_item_1->bins_per_angstrom( bins_per_angstrom_in );
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
		TS_ASSERT_EQUALS( test_item_1->get_score( test_distribution ), expected_score );
	}

	void test_derived_decay_fxns() {
		// First check all parameters can be accessed accurately
		core::scoring::epr_deer::DEERDataOP data_obj( core::scoring::epr_deer::DEERDataOP( new core::scoring::epr_deer::DEERDecayData() ) );
		core::scoring::epr_deer::DEERDecayDataOP test_item_1 = utility::pointer::dynamic_pointer_cast< core::scoring::epr_deer::DEERDecayData >( data_obj );
		data_obj->bins_per_angstrom( 1.0 );
		std::map< core::Real, core::Real > test_trace;
		for ( core::Size i = 0; i <= 250; ++i ) {
			core::Real time_pt = 0.0 + ( 0.002 * i );
			core::Real signal = 1.0 - ( 0.001 * i );
			test_item_1->append_trace_data_and_calculate( time_pt, signal, 100 );
			test_trace[ time_pt ] = signal;
		}
		for ( core::Size i = 1; i <= 750; ++i ) {
			core::Real time_pt = 0.5 + ( 0.002 * i );
			core::Real signal = 0.75;
			test_item_1->append_trace_data_and_calculate( time_pt, signal, 100 );
			test_trace[ time_pt ] = signal;
		}
		TS_ASSERT( test_item_1->spin_val( 50, 0.0 ) > 0.0 );
		TS_ASSERT_EQUALS( test_trace, test_item_1->trace() );

		core::Real test_k = -0.2;
		core::Real test_mod = 0.01;
		core::Real test_noise = 0.02;
		std::pair< core::Real, core::Real > test_mod_depth_bounds = std::make_pair( 0.025, 0.5 );

		test_item_1->k_fit( test_k );
		test_item_1->modulation_depth_fit( test_mod );
		test_item_1->noise( test_noise );
		test_item_1->mod_depth_bounds( test_mod_depth_bounds.first, test_mod_depth_bounds.second );

		TS_ASSERT_EQUALS( test_k, test_item_1->k_fit() );
		TS_ASSERT_EQUALS( test_mod, test_item_1->modulation_depth_fit() );
		TS_ASSERT_EQUALS( test_noise, test_item_1->noise() );
		TS_ASSERT_EQUALS( test_mod_depth_bounds, test_item_1->mod_depth_bounds() );

		std::map< core::Size, core::Real > test_histogram;
		for ( core::Size dist = 20; dist < 50; ++dist ) {
			test_histogram[ dist ] = 0.2;
		}
		TS_ASSERT( data_obj->get_score( test_histogram ) > 0.0 );
		TS_ASSERT_DELTA( test_item_1->k_fit(), -0.032, 0.01 );
		TS_ASSERT_DELTA( test_item_1->modulation_depth_fit(), 0.205, 0.01 );
	}
};
