// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragment/BBTorsionSRFD.cxxtest.hh
/// @author Christopher Miles (cmiles@uw.edu)

// Maximum allowable delta between a pair of floating point values
// to consider them equal
#define DELTA 0.0001

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Project headers
#include <core/types.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragData.fwd.hh>
#include <core/fragment/FragID.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FrameList.hh>
#include <core/fragment/SingleResidueFragData.fwd.hh>

// C/C++ headers
#include <iostream>

//Auto Headers
#include <core/fragment/FrameIteratorWorker_.hh>
#include <utility/vector1.hh>


using namespace core::fragment;
using namespace std;

class BBTorsionSRFDTest : public CxxTest::TestSuite {
 public:
	void setUp() {
		// contains cartesian coordinates
		core_init();
		fragments_ = FragmentIO().read_data("core/fragment/aat049603_05.200_v1_3.gz");
		fragments_->region_simple( fragments_->min_pos(), fragments_->max_pos(), frames_ );
	}

	void tearDown() {}

	// unit tests
	void test_cartesian_coordinates();
	void test_cartesian_coordinates_iter();

 private:
	/// @brief Fragments that contain cartesian coordinates
	FragSetOP fragments_;

	/// @brief Frames corresponding to <fragments_>
	FrameList frames_;
};

// Note: the terminology used below is nonstandard, but illustrative.
//
// A fragment file contains 1 Frame per residue position (1..N).
// Each Frame contains multiple Groups. One way to think about Groups
// is that they represent alternative fragments for a given position.
// Each Group contains multiple Entries. A fragment file containing
// 3-mers will have 3 Entries per Group. Alternatively, a fragment
// file containing 15-mers, will have 15 Entries per Group. To take a
// concrete example, a hypothetical fragment file containing 3-mers for
// a protein of length 2 (with 2 Groups per position) has structure:
//
//  Frame 1
//      Group 1
//          Entry 1
//          Entry 2
//          Entry 3
//      Group 2
//          Entry 1
//          Entry 2
//          Entry 3
//  Frame 2
//      Group 1
//          Entry 1
//          Entry 2
//          Entry 3
//      Group 2
//          Entry 1
//          Entry 2
//          Entry 3
void BBTorsionSRFDTest::test_cartesian_coordinates_iter() {
	using core::Real;
	using core::Size;

	for (ConstFrameIterator i = fragments_->begin(); i != fragments_->end(); ++i) {
		// foreach frame (P)
		FrameCOP frame = *i;
		Size fragment_length = frame->length();

		for (Size j = 1; j <= frame->nr_frags(); ++j) {
			// foreach group
			const FragData& group = frame->fragment(j);

			for (Size k = 1; k <= fragment_length; ++k) {
				// foreach entry
				SingleResidueFragDataCOP entry = group.get_residue(k);

				// explicit typecast required-- base class lacks accessors
				// for Cartesian coordinates
				BBTorsionSRFDCOP residue = static_cast< BBTorsionSRFD const * > (entry());
        Real x = residue->x();
        Real y = residue->y();
				Real z = residue->z();

				// we're not testing against known values here so much as illustrating the use
				// of the API and ensuring that the desired semantics have not been affected
				// through some other change to the source code. the assertions below simply
				// check that the floating point values are not NaN. According to the IEEE
				// standard, NaN values have the property that comparisons involving them are
				// always false.
				TS_ASSERT(x == x);
				TS_ASSERT(y == y);
				TS_ASSERT(z == z);
			}
		}
	}
  TS_ASSERT(1);
}

void BBTorsionSRFDTest::test_cartesian_coordinates() {
	FragID id = frames_.fragID(1);

	// Be very careful here. We must use fragment_ptr() rather than
	// fragment() to get the polymorphic behavior desired below.
	FragDataCOP data = id.fragment_ptr();
	SingleResidueFragDataCOP residue = data->get_residue(1);
	BBTorsionSRFDCOP derived_residue = static_cast< core::fragment::BBTorsionSRFD const * > (residue());

	// Check pdb id and chain. Generally speaking, you cannot assume that pdbid() and chain()
	// have been populated with their correct values. This owes to the myriad ways one can read
	// in a fragment file. In this case, the calls to FragmentIO::read_data() in the constructor
	// ensure that this data is present.
	TS_ASSERT_EQUALS("1mwm", data->pdbid());
	TS_ASSERT_EQUALS('A', data->chain());

	// has cartesian coords?
	TS_ASSERT(derived_residue->has_coordinates());

	// correct values
	TS_ASSERT_DELTA(derived_residue->x(), 28.590, DELTA);
  TS_ASSERT_DELTA(derived_residue->y(), -2.030, DELTA);
  TS_ASSERT_DELTA(derived_residue->z(), 115.920, DELTA);
}
