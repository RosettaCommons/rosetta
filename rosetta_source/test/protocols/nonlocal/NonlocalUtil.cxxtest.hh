// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/NonlocalUtil.cxxtest.hh
/// @brief test suite for protocols/nonlocal/util
/// @author Christopher Miles (cmiles@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// C/C++ headers
#include <string>

// Project headers
#include <core/types.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/util.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/nonlocal/util.hh>
#include <utility/vector1.hh>

namespace {
using core::Size;
using core::sequence::SequenceAlignment;
using protocols::loops::Loop;
using protocols::loops::Loops;
using std::string;
using utility::vector1;

class NonlocalUtilTest : public CxxTest::TestSuite {
 public:
  vector1<SequenceAlignment> alignments_;
  Loops regions;

  void setUp() {
    protocols_init();

    vector1<string> filenames;
    filenames.push_back("protocols/nonlocal/alignment.filt");

    core::sequence::read_all_alignments("grishin",
                                        filenames,
                                        &alignments_);
  }

  void test_find_regions() {
    const SequenceAlignment& alignment = alignments_[1];
    const Size unaligned_region_min_chunk_sz = 3;

    Loops aligned, unaligned;
    protocols::nonlocal::find_regions(alignment,
                                      unaligned_region_min_chunk_sz,
                                      &aligned,
                                      &unaligned);

    // Verify aligned
    TS_ASSERT_EQUALS(aligned.size(), 15);
    TS_ASSERT(aligned[1] == Loop(6, 30));
    TS_ASSERT(aligned[2] == Loop(34, 64));
    TS_ASSERT(aligned[3] == Loop(68, 82));
    TS_ASSERT(aligned[4] == Loop(86, 122));
    TS_ASSERT(aligned[5] == Loop(126, 151));
    TS_ASSERT(aligned[6] == Loop(155, 171));
    TS_ASSERT(aligned[7] == Loop(175, 206));
    TS_ASSERT(aligned[8] == Loop(210, 237));
    TS_ASSERT(aligned[9] == Loop(241, 250));
    TS_ASSERT(aligned[10] == Loop(254, 282));
    TS_ASSERT(aligned[11] == Loop(286, 291));
    TS_ASSERT(aligned[12] == Loop(295, 307));
    TS_ASSERT(aligned[13] == Loop(313, 334));
    TS_ASSERT(aligned[14] == Loop(338, 352));
    TS_ASSERT(aligned[15] == Loop(356, 370));

    // Verify unaligned
    TS_ASSERT_EQUALS(unaligned.size(), 15);
    TS_ASSERT(unaligned[1] == Loop(1, 5));
    TS_ASSERT(unaligned[2] == Loop(31, 33));
    TS_ASSERT(unaligned[3] == Loop(65, 67));
    TS_ASSERT(unaligned[4] == Loop(83, 85));
    TS_ASSERT(unaligned[5] == Loop(123, 125));
    TS_ASSERT(unaligned[6] == Loop(152, 154));
    TS_ASSERT(unaligned[7] == Loop(172, 174));
    TS_ASSERT(unaligned[8] == Loop(207, 209));
    TS_ASSERT(unaligned[9] == Loop(238, 240));
    TS_ASSERT(unaligned[10] == Loop(251, 253));
    TS_ASSERT(unaligned[11] == Loop(283, 285));
    TS_ASSERT(unaligned[12] == Loop(292, 294));
    TS_ASSERT(unaligned[13] == Loop(308, 312));
    TS_ASSERT(unaligned[14] == Loop(335, 337));
    TS_ASSERT(unaligned[15] == Loop(353, 356));
  }
};
}  // anonymous namespace
