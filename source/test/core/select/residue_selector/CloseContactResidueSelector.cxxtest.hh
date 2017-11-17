// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/residue_selector/CloseContactResidueSelector.cxxtest.hh
/// @brief test suite for core::select::residue_selector::CloseContactResidueSelector
/// @author Robert Lindner (rlindner@mpimf-heidelberg.mpg.de)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/select/residue_selector/utilities_for_testing.hh>

// Package headers
#include <core/select/residue_selector/CloseContactResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
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

static basic::Tracer TR("core.select.residue_selector.CloseContactResidueSelectorTests");


class CloseContactResidueSelectorTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
		trpcage = create_trpcage_ideal_poseop();
		//core::scoring::ScoreFunctionOP score = core::scoring::get_score_function();
		//score->score(trpcage);
	}

	void tearDown() {
		trpcage = 0;
	}

	/// @brief Test NotResidueSelector::parse_my_tag
	void test_CloseContactResidueSelector_parse_my_tag_selector() {
		std::string tag_string = "<CloseContact name=\"close_contact_w_trp6\" residue_selector=\"trp6\" contact_threshold=\"2.5\"/>";
		utility::tag::TagOP tag = utility::tag::Tag::create( tag_string );
		basic::datacache::DataMap dm;
		ResidueIndexSelectorOP  trp6_rs( new ResidueIndexSelector );
		trp6_rs->set_index( "6" );
		dm.add( "ResidueSelector", "trp6", trp6_rs );

		ResidueSelectorOP cc_rs( new CloseContactResidueSelector );

		cc_rs->parse_my_tag( tag, dm );
		ResidueSubset subset = cc_rs->apply( *trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage->size() );

		// check the result
		// 1. generate fake focus
		utility::vector1< core::Size > test_central_residues(trpcage->size(), false );
		test_central_residues[ 6 ] = true;
		// false
		TS_ASSERT( check_calculation( *trpcage, subset, test_central_residues, 2.5 ) );
	}


	//void dont_test_CloseContactResidueSelector_parse_my_tag_str() {
	// std::string tag_string = "<CloseContact name=neighbor_rs resnums=2,3,5 distance=5.2/>";
	// std::stringstream ss( tag_string );
	// utility::tag::TagOP tag( new utility::tag::Tag() );
	// tag->read( ss );
	// basic::datacache::DataMap dm;
	//
	// ResidueSelectorOP neighbor_rs( new CloseContactResidueSelector );
	// neighbor_rs->parse_my_tag( tag, dm );
	//
	// ResidueSubset subset = neighbor_rs->apply( trpcage );
	//
	// utility::vector1< core::Size > testFocus(trpcage.size(), false);
	// testFocus[2] = true;
	// testFocus[3] = true;
	// testFocus[5] = true;
	//
	// TS_ASSERT( check_calculation( trpcage, subset, testFocus, 5.2 ) );
	//
	//}
	//
	//// make sure we fail if neither selector nor focus string are provided
	//void dont_test_NeighbohoodResidueSelector_fail_no_focus() {
	// std::string tag_string = "<CloseContact name=neighbor_rs distance=5.2/>";
	// std::stringstream ss( tag_string );
	// utility::tag::TagOP tag( new utility::tag::Tag() );
	// tag->read( ss );
	// basic::datacache::DataMap dm;
	//
	// ResidueSelectorOP neighbor_rs( new CloseContactResidueSelector );
	// try {
	//  neighbor_rs->parse_my_tag( tag, dm );
	//  TS_ASSERT( false ); //parsing should fail!
	// } catch (utility::excn::Exception e ) {
	//  TS_ASSERT(true == true); //We should always get here.
	// }
	//}
	//
	//
	//// desired behavior is that the most recent call to set_focus or set_focus_selector
	//// determines which source of focus residues is used
	//void dont_test_CloseContactResidueSelector_use_last_provided_source_of_focus() {
	//
	// utility::vector1< core::Size > focus_set(trpcage.size(), false);
	// focus_set[2] = true;
	// focus_set[3] = true;
	// CloseContactResidueSelectorOP neighbor_rs( new CloseContactResidueSelector(focus_set, 5.0) );
	// ResidueSelectorOP odd_rs( new OddResidueSelector );
	//
	// ResidueSubset subset( trpcage.size(), false );
	// TS_ASSERT_EQUALS( subset.size(), trpcage.size() );
	//
	// utility::vector1< core::Size > testFocus_odd(trpcage.size(), false);
	// for ( core::Size ii = 1; ii <= trpcage.size(); ii += 2 ) {
	//  testFocus_odd[ ii ] = true;
	// }
	//
	// try {
	//  subset = neighbor_rs->apply( trpcage );
	//  TS_ASSERT( check_calculation( trpcage, subset, focus_set, 5.0 ) );
	//
	//  neighbor_rs->set_focus_selector( odd_rs );
	//  subset = neighbor_rs->apply( trpcage );
	//  TS_ASSERT( check_calculation( trpcage, subset, testFocus_odd, 5.0 ) );
	//
	//  neighbor_rs->set_focus( focus_set );
	//  subset = neighbor_rs->apply( trpcage );
	//  TS_ASSERT( check_calculation( trpcage, subset, focus_set, 5.0 ) );
	//
	// } catch (utility::excn::Exception e ) {
	//  std::cerr << "Exception! " << e.msg();
	//  TS_ASSERT( false );
	// }
	//}

	bool
	check_calculation(
		core::pose::Pose const & pose,
		ResidueSubset const & subset,
		ResidueSubset const & focus,
		core::Real distance_cutoff
	)
	{

		if ( focus.size() != subset.size() ) {
			return false;
		}

		// JAB - rewrite to simplify logic and change to ResidueSubset as focus.
		// We measure neighbors to the focus and this is the ctrl_subset.
		// We then compare this subset to our actual subset.

		ResidueSubset ctrl_subset(subset.size(), false);
		core::Real const cut2 = distance_cutoff * distance_cutoff;

		for ( core::Size ii = 1; ii <= focus.size(); ++ii ) {
			//If we have a focus residue, we will calculate its neighbors.
			if ( !focus[ii] ) continue;

			ctrl_subset[ii] = true;
			core::conformation::Residue const & ii_res( pose.residue( ii ) );

			//Go through all Subset residues
			for ( core::Size jj = 1; jj <= subset.size(); ++jj ) {

				//If this is already true, we don't need to recalculate.
				if ( ctrl_subset[ jj ] ) continue;

				core::conformation::Residue const & jj_res( pose.residue( jj ) );
				//Measure the distance for all atom pairs, set the ctrl_subset.
				for ( core::Size kk = 1; kk <= ii_res.natoms(); ++kk ) {
					for ( core::Size ll = 1; ll <= jj_res.natoms(); ++ll ) {
						if ( ii_res.xyz( kk ).distance_squared( jj_res.xyz( ll ) ) < cut2 ) {
							ctrl_subset[ jj ] = true;
							break;
						}
					}
					if ( ctrl_subset[ jj ] ) break;
				}
			}
		}

		TR << "focus " << utility::to_string(focus) << std::endl;
		TR << "subset" << utility::to_string(subset) << std::endl;
		TR << "ctrl  " << utility::to_string(ctrl_subset) << std::endl;

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
	core::pose::PoseOP trpcage;

};
