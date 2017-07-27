// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/residue_selector/NeighborhoodResidueSelector.cxxtest.hh
/// @brief test suite for core::select::residue_selector::NeighborhoodResidueSelector
/// @author Robert Lindner (rlindner@mpimf-heidelberg.mpg.de)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/select/residue_selector/DummySelectors.hh>
#include <test/core/select/residue_selector/utilities_for_testing.hh>

// Package headers
#include <core/select/residue_selector/NeighborhoodResidueSelector.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <string>

using namespace core::select::residue_selector;

static THREAD_LOCAL basic::Tracer TR("core.select.residue_selector.NeighborhoodResidueSelectorTests");


class NeighborhoodResidueSelectorTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
		trpcage = create_trpcage_ideal_pose();
		core::scoring::ScoreFunctionOP score = core::scoring::get_score_function();
		score->score(trpcage);
	}

	/// @brief Test NotResidueSelector::parse_my_tag
	void test_NeighborhoodResidueSelector_parse_my_tag_selector() {
		std::string tag_string = "<Neighborhood name=neighbor_rs selector=odd distance=5.2/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;
		ResidueSelectorOP odd_rs( new OddResidueSelector );
		dm.add( "ResidueSelector", "odd", odd_rs );

		ResidueSelectorOP neighbor_rs( new NeighborhoodResidueSelector );

		neighbor_rs->parse_my_tag( tag, dm );
		ResidueSubset subset = neighbor_rs->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.size() );

		// check the result
		// 1. generate fake focus
		utility::vector1< core::Size > testFocus(trpcage.size(), false);
		for ( core::Size ii = 1; ii <= trpcage.size(); ii += 2 ) {
			testFocus[ ii ] = true;
		}
		// test
		TS_ASSERT( check_calculation( trpcage, subset, testFocus, 5.2 ) );
	}


	void test_NeighborhoodResidueSelector_parse_my_tag_str() {
		std::string tag_string = "<Neighborhood name=neighbor_rs resnums=2,3,5 distance=5.2/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP neighbor_rs( new NeighborhoodResidueSelector );
		neighbor_rs->parse_my_tag( tag, dm );

		ResidueSubset subset = neighbor_rs->apply( trpcage );

		utility::vector1< core::Size > testFocus(trpcage.size(), false);
		testFocus[2] = true;
		testFocus[3] = true;
		testFocus[5] = true;

		TS_ASSERT( check_calculation( trpcage, subset, testFocus, 5.2 ) );

	}

	// Different sets of code run when we include focus vs not
	void test_NeighborhoodResidueSelector_dont_include_focus_close() {
		std::string tag_string =
			"<Neighborhood name=neighbor_rs resnums=2,3,5 distance=5.2 include_focus_in_subset=false/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP neighbor_rs( new NeighborhoodResidueSelector );
		neighbor_rs->parse_my_tag( tag, dm );

		ResidueSubset subset = neighbor_rs->apply( trpcage );

		utility::vector1< core::Size > testFocus(trpcage.size(), false);
		testFocus[2] = true;
		testFocus[3] = true;
		testFocus[5] = true;

		TS_ASSERT( check_calculation( trpcage, subset, testFocus, 5.2, false ) );

	}

	// Additionally, different code runs when distance > 10
	void test_NeighborhoodResidueSelector_dont_include_focus_far() {
		std::string tag_string =
			"<Neighborhood name=neighbor_rs resnums=2,3,5 distance=10.1 include_focus_in_subset=false/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP neighbor_rs( new NeighborhoodResidueSelector );
		neighbor_rs->parse_my_tag( tag, dm );

		ResidueSubset subset = neighbor_rs->apply( trpcage );

		utility::vector1< core::Size > testFocus(trpcage.size(), false);
		testFocus[2] = true;
		testFocus[3] = true;
		testFocus[5] = true;

		TS_ASSERT( check_calculation( trpcage, subset, testFocus, 10.1, false ) );

	}

	// make sure we fail if neither selector nor focus string are provided
	void test_NeighbohoodResidueSelector_fail_no_focus() {
		std::string tag_string = "<Neighborhood name=neighbor_rs distance=5.2/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP neighbor_rs( new NeighborhoodResidueSelector );
		try {
			neighbor_rs->parse_my_tag( tag, dm );
			TS_ASSERT( false ); //parsing should fail!
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			TS_ASSERT(true == true); //We should always get here.
		}
	}


	// desired behavior is that the most recent call to set_focus or set_focus_selector
	// determines which source of focus residues is used
	void test_NeighborhoodResidueSelector_use_last_provided_source_of_focus() {

		utility::vector1< core::Size > focus_set(trpcage.size(), false);
		focus_set[2] = true;
		focus_set[3] = true;
		NeighborhoodResidueSelectorOP neighbor_rs( new NeighborhoodResidueSelector(focus_set, 5.0) );
		ResidueSelectorOP odd_rs( new OddResidueSelector );

		ResidueSubset subset( trpcage.size(), false );
		TS_ASSERT_EQUALS( subset.size(), trpcage.size() );

		utility::vector1< core::Size > testFocus_odd(trpcage.size(), false);
		for ( core::Size ii = 1; ii <= trpcage.size(); ii += 2 ) {
			testFocus_odd[ ii ] = true;
		}

		try {
			subset = neighbor_rs->apply( trpcage );
			TS_ASSERT( check_calculation( trpcage, subset, focus_set, 5.0 ) );

			neighbor_rs->set_focus_selector( odd_rs );
			subset = neighbor_rs->apply( trpcage );
			TS_ASSERT( check_calculation( trpcage, subset, testFocus_odd, 5.0 ) );

			neighbor_rs->set_focus( focus_set );
			subset = neighbor_rs->apply( trpcage );
			TS_ASSERT( check_calculation( trpcage, subset, focus_set, 5.0 ) );

		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Exception! " << e.msg();
			TS_ASSERT( false );
		}
	}

	void test_NeighborhoodResidueSelector_code_level_acces() {


		NeighborhoodResidueSelectorOP neighbor_rs( new NeighborhoodResidueSelector );
		utility::vector1< core::Size > testFocus(trpcage.size(), false);
		testFocus[2] = true;
		testFocus[3] = true;
		testFocus[5] = true;

		neighbor_rs->set_focus(testFocus);


		utility::vector1< core::Real > distances;
		distances.push_back(2.0);
		distances.push_back(4.0);
		distances.push_back(6.0);
		distances.push_back(8.0);
		distances.push_back(10.0);
		distances.push_back(12.0);


		neighbor_rs->set_include_focus_in_subset(false);
		for ( core::Real d : distances ) {
			//TR << "Testing distance " << d << " no focus in subset " << std::endl;
			neighbor_rs->set_distance(d);
			ResidueSubset subset = neighbor_rs->apply( trpcage );
			TS_ASSERT( check_calculation( trpcage, subset, testFocus, d, false ) );
		}

		neighbor_rs->set_include_focus_in_subset(true);
		for ( core::Real d : distances ) {
			//TR << "Testing distance "  << d << " focus in subset " << std::endl;
			neighbor_rs->set_distance(d);
			ResidueSubset subset = neighbor_rs->apply( trpcage );
			TS_ASSERT( check_calculation( trpcage, subset, testFocus, d, true ) );
		}


	}

	bool
	check_calculation( core::pose::Pose const & pose,
		ResidueSubset const & subset,
		ResidueSubset const & focus,
		core::Real distance,
		bool include_focus_in_subset = true )
	{

		if ( focus.size() != subset.size() ) {
			return false;
		}

		//JAB - rewrite to simplify logic and change to ResidueSubset as focus.
		/// We measure neighbors to the focus and this is the ctrl_subset.
		///  We then compare this subset to our actual subset.

		ResidueSubset ctrl_subset(subset.size(), false);
		core::Real const dst_squared = distance * distance;

		for ( core::Size i = 1; i <= focus.size(); ++i ) {
			//If we have a focus residue, we will calculate its neighbors.
			if ( !focus[i] ) continue;

			ctrl_subset[i] = true;
			core::conformation::Residue const & r1( pose.residue( i ) );

			//Go through all Subset residues
			for ( core::Size x = 1; x <= subset.size(); ++x ) {

				//If this is already true, we don't need to recalculate.
				if ( ctrl_subset[ x ] ) continue;

				//Measure the distance, set the ctrl_subset.
				core::conformation::Residue const & r2( pose.residue( x ) );
				core::Real const d_sq( r1.xyz( r1.nbr_atom() ).distance_squared( r2.xyz( r2.nbr_atom() ) ) );
				if ( d_sq <= dst_squared ) {
					ctrl_subset[ x ] = true;
				}

			}


		}

		// bcov - Putting this after rather than inside the above loop
		//  to be absolutely sure we get it correct

		if ( ! include_focus_in_subset ) {
			for ( core::Size i = 1; i <= ctrl_subset.size(); i++ ) {
				if ( focus[i] ) {
					ctrl_subset[i] = false;
				}
			}
		}


		TR<< "focus " << utility::to_string(focus) << std::endl;
		TR<< "subset" << utility::to_string(subset) << std::endl;
		TR<< "ctrl  " << utility::to_string(ctrl_subset) << std::endl;

		//Compare subset to control subset. Return False if they do not match.
		for ( core::Size i = 1; i <= pose.size(); ++i ) {
			if ( ctrl_subset[ i ]  != subset[ i ] ) {
				TR << "Resnum "<< i <<" "<< ctrl_subset[ i ] << "!=" << subset[ i ] << std::endl;
				return false;
			}
		}

		// no mismatches found
		return true;
	}

private:
	core::pose::Pose trpcage;

};
