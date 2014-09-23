// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/Chunk.cxxtest.hh
/// @brief test suite for protocols/nonlocal/Chunk
/// @author Christopher Miles (cmiles@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// External headers
#include <boost/scoped_ptr.hpp>

// Project headers
#include <core/types.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/nonlocal/Chunk.hh>
#include <protocols/nonlocal/Region.hh>

//Auto Headers
#include <utility/vector1.hh>


namespace {
using core::Size;
using core::kinematics::MoveMap;
using core::kinematics::MoveMapOP;
using protocols::nonlocal::Chunk;
using protocols::nonlocal::Region;
using protocols::nonlocal::RegionOP;

class ChunkTest : public CxxTest::TestSuite {
 public:
	MoveMapOP movable_;
	boost::scoped_ptr<Chunk> valid_chunk_;
	boost::scoped_ptr<Chunk> invalid_chunk_;

	void setUp() {
    core_init();

    // define the movable backbone degrees of freedom
		movable_ = new MoveMap();
		movable_->set_bb(true);
		for (Size i = 5; i <= 10; ++i)
			movable_->set_bb(i, false);

    valid_chunk_.reset(new Chunk(new Region(1,5), movable_));
		invalid_chunk_.reset(new Chunk(new Region(5,10), movable_));
	}

	void tearDown() {}

	/// @brief Ensure that a Chunk is able to identify situations where it is
	/// unable to make a move. This situation occurs when there are no movable
	/// backbone degrees of freedom in the contiguous stretch of sequence
	/// defining the Chunk.
	void test_validity() {
		TS_ASSERT(valid_chunk_->valid());
		TS_ASSERT(!invalid_chunk_->valid());
	}

	void test_start_stop() {
		TS_ASSERT_EQUALS(valid_chunk_->start(), 1);
		TS_ASSERT_EQUALS(valid_chunk_->stop(), 5);
		TS_ASSERT_EQUALS(invalid_chunk_->start(), 5);
		TS_ASSERT_EQUALS(invalid_chunk_->stop(), 10);
	}

	/// @brief Ensure that all insertion positions selected in choose() are on the
	/// closed interval [start(), stop()].
	void test_samples_in_range() {
    for (int i = 0; i < 10000; ++i) {
      Size insert_pos = valid_chunk_->choose();
			TS_ASSERT(insert_pos >= valid_chunk_->start());
			TS_ASSERT(insert_pos <= valid_chunk_->stop());
		}
	}
};
}  // anonymous namespace
