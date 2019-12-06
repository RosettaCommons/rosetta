// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/residue_selector/SliceResidueSelector.cxxtest.hh
/// @brief test suite for core::select::residue_selector::SliceResidueSelector
/// @author Brian Coventry


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/select/residue_selector/DummySelectors.hh>
#include <test/core/select/residue_selector/utilities_for_testing.hh>

// Package headers
#include <core/select/residue_selector/SliceResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/conformation/Residue.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
#include <utility/pointer/memory.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <string>

using namespace core::select::residue_selector;

static basic::Tracer TR("core.select.residue_selector.SliceResidueSelectorTests");


class SliceResidueSelectorTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
		core::pose::make_pose_from_sequence( pose, "AAAAA", "fa_standard", true );
	}

	void select_13_test( SliceResidueSelector const & sel, core::Size testno ) {
		ResidueSubset result = sel.apply( pose );

		TR << "Test " << testno << std::endl;
		TS_ASSERT( result[1] );
		TS_ASSERT( !result[2] );
		TS_ASSERT( result[3] );
		TS_ASSERT( !result[4] );
		TS_ASSERT( !result[5] );


	}

	void select_123_test( SliceResidueSelector const & sel, core::Size testno ) {
		ResidueSubset result = sel.apply( pose );

		TR << "Test " << testno << std::endl;
		TS_ASSERT( result[1] );
		TS_ASSERT( result[2] );
		TS_ASSERT( result[3] );
		TS_ASSERT( !result[4] );
		TS_ASSERT( !result[5] );

	}

	void select_24_test( SliceResidueSelector const & sel, core::Size testno ) {
		ResidueSubset result = sel.apply( pose );

		TR << "Test " << testno << std::endl;
		TS_ASSERT( !result[1] );
		TS_ASSERT( result[2] );
		TS_ASSERT( !result[3] );
		TS_ASSERT( result[4] );
		TS_ASSERT( !result[5] );

	}

	void select_234_test( SliceResidueSelector const & sel, core::Size testno ) {
		ResidueSubset result = sel.apply( pose );

		TR << "Test " << testno << std::endl;
		TS_ASSERT( !result[1] );
		TS_ASSERT( result[2] );
		TS_ASSERT( result[3] );
		TS_ASSERT( result[4] );
		TS_ASSERT( !result[5] );

	}

	void except_test( SliceResidueSelector const & sel, core::Size testno ) {

		TR << "Test " << testno << std::endl;
		try {
			ResidueSubset result = sel.apply( pose );

			TS_ASSERT( false );
		} catch ( utility::excn::Exception const & e ) {

		}
	}

	/// This function tests all the functionality
	void test_SliceResidueSelector_constructors() {

		ResidueIndexSelectorCOP select_13 = utility::pointer::make_shared<ResidueIndexSelector>( "1,3" );

		SliceResidueSelector sel;
		sel = SliceResidueSelector( select_13, slice_enums::SPARSE, 1, -1 );
		select_13_test( sel, 1 );
		sel = SliceResidueSelector( select_13, slice_enums::SPARSE, 1, 2 );
		select_13_test( sel, 2 );

		sel = SliceResidueSelector( select_13, slice_enums::SPARSE, utility::vector1<int> { 1, -1 } );
		select_13_test( sel, 3 );
		sel = SliceResidueSelector( select_13, slice_enums::SPARSE, utility::vector1<int> { 1, 2 } );
		select_13_test( sel, 4 );

		sel = SliceResidueSelector( select_13, slice_enums::SPARSE, 1, -3 );
		except_test( sel, 5 );
		sel = SliceResidueSelector( select_13, slice_enums::SPARSE, 2, -2 );
		except_test( sel, 6 );
		sel = SliceResidueSelector( select_13, slice_enums::SPARSE, 3, 3 );
		except_test( sel, 7 );

		sel = SliceResidueSelector( select_13, slice_enums::SPARSE, utility::vector1<int> { 1, -3 } );
		except_test( sel, 8 );
		sel = SliceResidueSelector( select_13, slice_enums::SPARSE, utility::vector1<int> { 1, 3 } );
		except_test( sel, 9 );

		ResidueIndexSelectorCOP select_14 = utility::pointer::make_shared<ResidueIndexSelector>( "1,3,4" );
		ResidueIndexSelectorCOP select_134 = utility::pointer::make_shared<ResidueIndexSelector>( "1,3,4" );

		sel = SliceResidueSelector( select_14, slice_enums::CONTIGUOUS, 1, -2 );
		select_123_test( sel, 10 );
		sel = SliceResidueSelector( select_14, slice_enums::CONTIGUOUS, 1, 3 );
		select_123_test( sel, 11 );
		sel = SliceResidueSelector( select_14, slice_enums::CONTIGUOUS, -4, 3 );
		select_123_test( sel, 12 );

		sel = SliceResidueSelector( select_14, slice_enums::CONTIGUOUS, utility::vector1<int> { 1, -2 } );
		select_13_test( sel, 13 );
		sel = SliceResidueSelector( select_14, slice_enums::CONTIGUOUS, utility::vector1<int> { 1, 3 } );
		select_13_test( sel, 14 );

		sel = SliceResidueSelector( select_14, slice_enums::CONTIGUOUS, 5, 5 );
		except_test( sel, 15 );


		// sel = SliceResidueSelector( select_13, slice_enums::CONTIGUOUS_AND, 1, -1 );
		// select_13_test( sel, 16 );
		// sel = SliceResidueSelector( select_134, slice_enums::CONTIGUOUS_AND, 1, -2 );
		// select_13_test( sel, 17 );


		ResidueIndexSelectorCOP select_24 = utility::pointer::make_shared<ResidueIndexSelector>( "2,4" );
		ResidueIndexSelectorCOP select_25 = utility::pointer::make_shared<ResidueIndexSelector>( "2,4,5" );
		// ResidueIndexSelectorCOP select_245 = utility::pointer::make_shared<ResidueIndexSelector>( "2,4,5" );

		sel = SliceResidueSelector( select_24, slice_enums::SPARSE, 1, -1 );
		select_24_test( sel, 18 );
		sel = SliceResidueSelector( select_24, slice_enums::SPARSE, 1, 2 );
		select_24_test( sel, 19 );

		sel = SliceResidueSelector( select_24, slice_enums::SPARSE, utility::vector1<int> { 1, -1 } );
		select_24_test( sel, 20 );
		sel = SliceResidueSelector( select_24, slice_enums::SPARSE, utility::vector1<int> { 1, 2 } );
		select_24_test( sel, 21 );

		sel = SliceResidueSelector( select_25, slice_enums::CONTIGUOUS, 1, -2 );
		select_234_test( sel, 22 );
		sel = SliceResidueSelector( select_25, slice_enums::CONTIGUOUS, 1, 3 );
		select_234_test( sel, 23 );
		sel = SliceResidueSelector( select_25, slice_enums::CONTIGUOUS, -4, 3 );
		select_234_test( sel, 24 );

		sel = SliceResidueSelector( select_25, slice_enums::CONTIGUOUS, utility::vector1<int> { 1, -2 } );
		select_24_test( sel, 25 );
		sel = SliceResidueSelector( select_25, slice_enums::CONTIGUOUS, utility::vector1<int> { 1, 3 } );
		select_24_test( sel, 26 );

		sel = SliceResidueSelector( select_25, slice_enums::CONTIGUOUS, 5, 5 );
		except_test( sel, 27 );


		// sel = SliceResidueSelector( select_24, slice_enums::CONTIGUOUS_AND, 1, -1 );
		// select_24_test( sel, 28 );
		// sel = SliceResidueSelector( select_245, slice_enums::CONTIGUOUS_AND, 1, -2 );
		// select_24_test( sel, 29 );


	}

private:
	core::pose::Pose pose;

};
