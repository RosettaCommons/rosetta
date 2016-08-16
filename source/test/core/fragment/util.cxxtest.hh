// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragment/util.cxxtest.hh
/// @brief  test suite for core::fragment::util.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Project Headers
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/util.hh>

//Auto Headers
#include <core/fragment/FrameIteratorWorker_.hh>
#include <utility/vector1.hh>


namespace {

using core::Size;
using core::fragment::FragSetOP;
using core::fragment::Frame;
using core::fragment::FrameIterator;
using core::fragment::ConstFrameIterator;

class UtilTest : public CxxTest::TestSuite {
public:

	/// @bried Top k fragments
	Size k;

	/// @brief Set of fragments read in from file
	FragSetOP all_fragments_;

	void setUp() {
		core_init();
		k = 10;
		all_fragments_ = core::fragment::FragmentIO().read_data("core/fragment/aat049603_05.200_v1_3.gz");
	}

	/// @brief Verifies the functionality of core::fragment::retain_top()
	void test_filter() {
		FragSetOP fragments = all_fragments_->clone();
		core::fragment::retain_top(k, fragments);

		// identical metadata...
		TS_ASSERT_EQUALS(all_fragments_->min_pos(), fragments->min_pos());
		TS_ASSERT_EQUALS(all_fragments_->max_pos(), fragments->max_pos());
		TS_ASSERT_EQUALS(all_fragments_->max_frag_length(), fragments->max_frag_length());

		// but reduced number of fragments per frame
		for ( ConstFrameIterator i = fragments->begin(); i != fragments->end(); ++i ) {
			TS_ASSERT_LESS_THAN_EQUALS((*i)->nr_frags(), k);
		}
	}
};
}  // anonymous namespace
