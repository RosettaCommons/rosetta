// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/select/residue_selector/StoredResidueSubsetSelector.cxxtest.hh
/// @brief test suite for core::select::residue_selector::StoredResidueSubsetSelector
/// @author Tom Linsky (tlinsky@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/select/residue_selector/DummySelectors.hh>

// Package headers
#include <protocols/residue_selectors/StoredResidueSubsetSelector.hh>
#include <protocols/residue_selectors/StoreResidueSubsetMover.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/select/residue_selector/ResidueNameSelector.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// C++ headers
#include <boost/assign.hpp>
#include <string>

using namespace core::select::residue_selector;

static THREAD_LOCAL basic::Tracer TR( "test.core.select.StoredResidueSubsetSelectorTests" );

class StoredResidueSubsetSelectorTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void test_StoredResidueSubsetSelector_parse_my_tag() {
		using namespace protocols::residue_selectors;
		std::stringstream ss;
		ss << "<StoredResidueSubset name=\"stored_subset\" />";
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		set_throw_on_next_assertion_failure();
		StoredResidueSubsetSelectorOP rs( new StoredResidueSubsetSelector );
		TS_ASSERT_THROWS_ANYTHING( rs->parse_my_tag( tag, dm ) );

		std::stringstream ssgood;
		ssgood << "<StoredResidueSubset name=\"good_ss\" subset_name=\"stored_subset\" />";
		tag->read( ssgood );
		TS_ASSERT_THROWS_NOTHING( rs->parse_my_tag( tag, dm ) );
	}

	void test_StoreResidueSubsetMover_parse_my_tag() {
		using namespace protocols::residue_selectors;
		std::stringstream ss;
		ss << "<StoreResidueSubset name=\"store_subset\" />";
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;
		protocols::moves::Movers_map movers;
		protocols::filters::Filters_map filters;

		StoreResidueSubsetMoverOP rs( new StoreResidueSubsetMover );
		set_throw_on_next_assertion_failure();
		TS_ASSERT_THROWS_ANYTHING( rs->parse_my_tag( tag, dm, filters, movers, core::pose::Pose() ) );

		std::stringstream ssgood;
		ssgood << "<StoreResidueSubset name=\"good_ss\" subset_name=\"stored_subset\" />";
		tag->read( ssgood );
		TS_ASSERT_THROWS_NOTHING( rs->parse_my_tag( tag, dm, filters, movers, core::pose::Pose() ) );
	}

	void test_store_and_retrieve_subset() {
		using namespace core::select::residue_selector;
		using namespace protocols::residue_selectors;

		// create comparision subset
		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ResidueSelectorCOP sel_trp( new ResidueNameSelector( "TRP" ) );
		ResidueSubset const subset = sel_trp->apply( trpcage );

		std::string const stored_subset_name = "stored_trp_selection";

		// before storing anything, the selector should fail
		StoredResidueSubsetSelector retrieve_stored( stored_subset_name );
		set_throw_on_next_assertion_failure();
		TS_ASSERT_THROWS_ANYTHING( retrieve_stored.apply( trpcage ) );

		// store subset into pose
		StoreResidueSubsetMover store_subset( sel_trp, stored_subset_name, false );
		store_subset.apply( trpcage );

		// seting twice should fail because overwrite=false
		set_throw_on_next_assertion_failure();
		TS_ASSERT_THROWS_ANYTHING( store_subset.apply( trpcage ) );

		// now that something is stored, this should work
		ResidueSubset const stored_subset = retrieve_stored.apply( trpcage );

		// should be fine to run this over and over
		TS_ASSERT_THROWS_NOTHING( retrieve_stored.apply( trpcage ) );

		TS_ASSERT_EQUALS( stored_subset.size(), subset.size() );
		core::Size idx = 1;
		for ( ResidueSubset::const_iterator a=subset.begin(), b=stored_subset.begin(); ( (a!=subset.end()) && (b!=stored_subset.end()) ); ++a, ++b, ++idx ) {
			TR << idx << ": " << *a << " " << *b << std::endl;
			TS_ASSERT( *a == *b );
		}

		// invalid subset name and this should fail
		StoredResidueSubsetSelector retrieve_stored_badname( "bad_name" );
		set_throw_on_next_assertion_failure();
		TS_ASSERT_THROWS_ANYTHING( retrieve_stored_badname.apply( trpcage ) );

		// copy the pose and this should still work great
		core::pose::Pose posecopy = trpcage;
		ResidueSubset const stored_subset2 = retrieve_stored.apply( posecopy );
		TS_ASSERT_EQUALS( stored_subset2.size(), subset.size() );
		idx = 1;
		for ( ResidueSubset::const_iterator a=subset.begin(), b=stored_subset.begin(); ( (a!=subset.end()) && (b!=stored_subset.end()) ); ++a, ++b, ++idx ) {
			TR << idx << ": " << *a << " " << *b << std::endl;
			TS_ASSERT( *a == *b );
		}

		// add residues to the pose and this should fail
		trpcage.append_residue_by_jump( trpcage.residue(2), 1 );
		set_throw_on_next_assertion_failure();
		TS_ASSERT_THROWS_ANYTHING( retrieve_stored.apply( trpcage ) );

	}

};
