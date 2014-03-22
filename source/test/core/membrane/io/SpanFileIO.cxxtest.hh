// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/io/SpanFileIO.cxxtest.hh
///
/// @brief 		Test Suite for izstream reader class that reads in OCTOPUS file data
/// @details
///
/// @author 	Rebecca Alford (rfalford12@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Tested Classes
#include <core/membrane/io/SpanFileIO.hh>

#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Exceptions.hh>

// Package Headers
#include <core/types.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <algorithm>

using namespace core::conformation::membrane;

/**
 * The format of these tests is to initialize a data structure from
 * the appropriate file and then compare it to its expected state with some comparator function that I will write
 * Any tests you add should follow thsi convention
 */

/// @brief Test Class: Span File loader
class SpanFileIOTests : public CxxTest::TestSuite {

public: // test methods

	/// @brief SetUp - Runs before each test
	void setUp()
	{

		core_init();

		// Create an Instance of the spanfile class
		sfio_ = new core::membrane::io::SpanFileIO();

		// Initialize Comparison Cases
		initialize_expected();

	}

	/// @brief tearDon - runs after each test
	void tearDown()
	{}

	/// Tests for IO Class ////

	/// @brief Test cases for valid inputs
	void test_emptyInput()
	{

		TS_TRACE("Testing response to empty or null input...");

		// Checking return null
        TS_ASSERT_THROWS_ANYTHING( sfio_->get_topology_from_spanfile("") );

	}

	void test_singleHelix()
	{

		TS_TRACE("Testing single tmh helix case...");

		// Checking single tmh helix case
		SpanningTopologyOP topology = sfio_->get_topology_from_spanfile("core/membrane/io/single_tmh.span");
		TS_ASSERT( compare_topology( topology, sp_single_tmh_ ) );

	}

	void test_multipleHelix()
	{

		TS_TRACE("Testing multiple helix case...");

		// Checking multiple helix case
		SpanningTopologyOP topology = sfio_->get_topology_from_spanfile("core/membrane/io/multiple_tmh.span");
		TS_ASSERT( compare_topology( topology, sp_multiple_tmh_ ) );

	}

	void test_edgeHelix()
	{

		TS_TRACE("Testing edge helix case...");

		// Checking edge helix case
		SpanningTopologyOP topology = sfio_->get_topology_from_spanfile("core/membrane/io/edge_tmh.span");
		TS_ASSERT( compare_topology( topology, sp_edge_tmh_ ) );

	}

	void test_fullHelix()
	{

		TS_TRACE("Testing full spanning helix case...");
		SpanningTopologyOP topology = sfio_->get_topology_from_spanfile("core/membrane/io/full_tmh.span");
		TS_ASSERT( compare_topology( topology, sp_full_tmh_ ) );

	}

	void test_negativeTotalRes()
	{

		TS_TRACE("Testing negative total res object...");
		TS_ASSERT_THROWS_ANYTHING( SpanningTopologyOP topology = sfio_->get_topology_from_spanfile("core/membrane/io/negative_total.span") );
	}

	void test_helixOutofBounds()
	{

		TS_TRACE("Testing helix out of bounds case...");
		TS_ASSERT_THROWS_ANYTHING( SpanningTopologyOP topology = sfio_->get_topology_from_spanfile("core/membrane/io/out_of_bounds.span"); );

	}

private: // helper testing methods

	/// @brief 	Compare Spanning Topology Data
	bool compare_topology( SpanningTopologyOP sp1, SpanningTopologyOP sp2 ) {

		// false until proven true**
		// by convention, test case should always be sp1 becasue sp2 in the above cases
		// are not fully constructed objects (ok by completeness)

		// Compare Total Values
		TS_ASSERT_EQUALS( sp1->total_residue_in_span_file(), sp2->total_residue_in_span_file() );
		TS_ASSERT_EQUALS( sp1->total_tmhelix(), sp2->total_tmhelix() );

		// Quick loop bounds checking
		TS_ASSERT( sp1->total_tmhelix() != 0 );

		// Check equal values in span and full span arrays
		for ( core::Size i = 1; i <= sp1->total_tmhelix(); ++ i ) {

			// Check spans are constructed properly
			TS_ASSERT_EQUALS( sp1->span()(i, 1), sp2->span()(i, 1) );
			TS_ASSERT_EQUALS( sp1->span()(i, 2), sp2->span()(i, 2) );

			// not checking full span because it has no use right now - legacy/dead code
		}

		// Check scoring has been initialized properly
		for ( core::Size i = 1; i <= sp1->total_residue_in_span_file(); ++i ) {
			TS_ASSERT( sp1->allow_scoring()[i]);
		}

		// Check that tmh scoring has been initialized properly
		for ( core::Size i = 1; i <= sp1->total_tmhelix(); ++i ) {
			TS_ASSERT( sp1->allow_tmh_scoring()[i]);
		}

		// Check that tm region array is filled correctly
       for( core::Size i = 1; i <= sp1->total_residue_in_span_file(); ++i ) {
            for ( core::Size reg1 = 1; reg1 <= sp1->total_tmhelix(); ++reg1 ) {

                if( i >= sp1->span()(reg1, 1) && i <= sp1->span()(reg1, 2) ) {
                    TS_ASSERT( sp1->tmregion()[i] );
                    continue;
                }

            }
        }

		return true;
	}

	/// @brief Initialize some structures with expected inputs to compare to
	void initialize_expected() {

		// Single Spanning TMH
		sp_single_tmh_ = new SpanningTopology();
		sp_single_tmh_->set_total_residue_in_spanfile( 185 );
		sp_single_tmh_->set_total_tmhelix( 1 );
        
        ObjexxFCL::FArray2D< Size > span1;
		span1.dimension( sp_single_tmh_->total_tmhelix(), 2 );
		span1(1, 1) = 4;
		span1(1, 2) = 24;
        sp_single_tmh_->set_span( span1 );
        
		// Multi-Spanning TMH (Rhodopsin)
		sp_multiple_tmh_ = new SpanningTopology();
		sp_multiple_tmh_->set_total_residue_in_spanfile( 354 );
		sp_multiple_tmh_->set_total_tmhelix( 7 );
        ObjexxFCL::FArray2D< Size > span2;
		span2.dimension( sp_multiple_tmh_->total_tmhelix(), 2 );
		span2(1, 1) = 38;
		span2(1, 2) = 58;
		span2(2, 1) = 75;
		span2(2, 2) = 95;
		span2(3, 1) = 113;
		span2(3, 2) = 133;
		span2(4, 1) = 152;
		span2(4, 2) = 172;
		span2(5, 1) = 203;
		span2(5, 2) = 223;
		span2(6, 1) = 253;
		span2(6, 2) = 273;
		span2(7, 1) = 286;
		span2(7, 2) = 306;
        sp_multiple_tmh_->set_span( span2 );

		// Helix Edge (Engineered Glycophorin A)
		sp_edge_tmh_ = new SpanningTopology();
		sp_edge_tmh_->set_total_residue_in_spanfile( 30 );
		sp_edge_tmh_->set_total_tmhelix( 1 );
        ObjexxFCL::FArray2D< Size > span3;
		span3.dimension( sp_edge_tmh_->total_tmhelix(), 2 );
		span3(1, 1) = 1;
		span3(1, 2) = 21;
        sp_edge_tmh_->set_span( span3 );
        
		// Full Helix Span (Engineered Glycophorin A)
		sp_full_tmh_ = new SpanningTopology();
		sp_full_tmh_->set_total_residue_in_spanfile( 21 );
		sp_full_tmh_->set_total_tmhelix( 1 );
        ObjexxFCL::FArray2D< Size > span4;
		span4.dimension( sp_full_tmh_->total_tmhelix(), 2 );
		span4(1, 1) = 1;
		span4(1, 2) = 21;
        sp_full_tmh_->set_span( span4 );
	}

private: // edge cases

	// Class to Test
	core::membrane::io::SpanFileIOOP sfio_;

	// Expected cases
	SpanningTopologyOP sp_single_tmh_;
	SpanningTopologyOP sp_multiple_tmh_;
	SpanningTopologyOP sp_edge_tmh_;
	SpanningTopologyOP sp_full_tmh_;

};
