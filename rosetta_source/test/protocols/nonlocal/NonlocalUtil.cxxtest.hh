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
  static const int MIN_CHUNK_SZ = 7;
  static const int MAX_CHUNK_SZ = 28;
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

    Loops aligned, unaligned;
    protocols::nonlocal::find_regions_with_minimum_size
        (alignment, MIN_CHUNK_SZ, &aligned, &unaligned);

    // Verify sizes
    for (Loops::const_iterator i = aligned.begin(); i != aligned.end(); ++i)
      TS_ASSERT(i->length() >= MIN_CHUNK_SZ);

    for (Loops::const_iterator i = unaligned.begin(); i != unaligned.end(); ++i)
      TS_ASSERT(i->length() >= MIN_CHUNK_SZ);

    // Verify aligned
    TS_ASSERT_EQUALS(aligned.size(), 13);
    TS_ASSERT(aligned[1] == Loop(8, 28));
    TS_ASSERT(aligned[2] == Loop(36, 62));
    TS_ASSERT(aligned[3] == Loop(70, 80));
    TS_ASSERT(aligned[4] == Loop(88, 120));
    TS_ASSERT(aligned[5] == Loop(128, 149));
    TS_ASSERT(aligned[6] == Loop(157, 169));
    TS_ASSERT(aligned[7] == Loop(177, 204));
    TS_ASSERT(aligned[8] == Loop(212, 235));
    TS_ASSERT(aligned[9] == Loop(256, 280));
    TS_ASSERT(aligned[10] == Loop(297, 306));
    TS_ASSERT(aligned[11] == Loop(314, 332));
    TS_ASSERT(aligned[12] == Loop(340, 350));
    TS_ASSERT(aligned[13] == Loop(358, 370));

    // Verify unaligned
    TS_ASSERT_EQUALS(unaligned.size(), 13);
    TS_ASSERT(unaligned[1] == Loop(1, 7));
    TS_ASSERT(unaligned[2] == Loop(29, 35));
    TS_ASSERT(unaligned[3] == Loop(63, 69));
    TS_ASSERT(unaligned[4] == Loop(81, 87));
    TS_ASSERT(unaligned[5] == Loop(121, 127));
    TS_ASSERT(unaligned[6] == Loop(150, 156));
    TS_ASSERT(unaligned[7] == Loop(170, 176));
    TS_ASSERT(unaligned[8] == Loop(205, 211));
    TS_ASSERT(unaligned[9] == Loop(236, 255));
    TS_ASSERT(unaligned[10] == Loop(281, 296));
    TS_ASSERT(unaligned[11] == Loop(307, 313));
    TS_ASSERT(unaligned[12] == Loop(333, 339));
    TS_ASSERT(unaligned[13] == Loop(351, 358));
  }

  void test_limit_chunk_size() {
    const SequenceAlignment& alignment = alignments_[1];

    Loops aligned, unaligned;
    protocols::nonlocal::find_regions_with_minimum_size
        (alignment, MIN_CHUNK_SZ, &aligned, &unaligned);

    protocols::nonlocal::limit_chunk_size(MIN_CHUNK_SZ, MAX_CHUNK_SZ, &aligned);
    for (Loops::const_iterator i = aligned.begin(); i != aligned.end(); ++i) {
      TS_ASSERT(i->length() >= MIN_CHUNK_SZ);
      TS_ASSERT(i->length() <= MAX_CHUNK_SZ);
    }

    protocols::nonlocal::limit_chunk_size(MIN_CHUNK_SZ, MAX_CHUNK_SZ, &unaligned);
    for (Loops::const_iterator i = unaligned.begin(); i != unaligned.end(); ++i) {
      const Loop& loop = *i;
      TS_ASSERT(i->length() >= MIN_CHUNK_SZ);
      TS_ASSERT(i->length() <= MAX_CHUNK_SZ);
    }
  }
};
}  // anonymous namespace
