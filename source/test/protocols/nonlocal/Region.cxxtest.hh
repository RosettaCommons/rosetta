// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/nonlocal/Region.cxxtest.hh
/// @brief test suite for protocols/nonlocal/Region
/// @author Christopher Miles (cmiles@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>

// Project headers
#include <protocols/nonlocal/Region.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/types.hh>
#include <protocols/nonlocal/Region.fwd.hh>
#include <utility/down_cast.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <numeric/xyzVector.fwd.hh>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <limits>


namespace {
using core::Size;
using protocols::nonlocal::Region;

class RegionTest : public CxxTest::TestSuite {
public:
	Region* region_;
	Region* region_backward_;

	void setUp() {
		region_ = new Region(13, 17);
		region_backward_ = new Region(17, 13);
	}

	void tearDown() {
		delete region_;
		delete region_backward_;
	}

	void test_getters() {
		TS_ASSERT_EQUALS(region_->start(), 13);
		TS_ASSERT_EQUALS(region_->stop(), 17);

		TS_ASSERT_EQUALS(region_backward_->start(), 17);
		TS_ASSERT_EQUALS(region_backward_->stop(), 13);
	}

	void test_length() {
		TS_ASSERT_EQUALS(region_->length(), 5);
		TS_ASSERT_EQUALS(region_backward_->length(), 5);
	}
};
}  // anonymous namespace
