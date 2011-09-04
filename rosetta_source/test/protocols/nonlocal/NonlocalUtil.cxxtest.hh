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
#include <iostream>
#include <string>

// Project headers
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/util.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/nonlocal/util.hh>
#include <utility/vector1.hh>

namespace {
using core::Size;
using core::id::SequenceMapping;
using core::sequence::SequenceAlignment;
using protocols::loops::Loop;
using protocols::loops::Loops;
using std::endl;
using std::string;
using utility::vector1;

static basic::Tracer TR("protocols.nonlocal.NonlocalUtilTest");

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

    // Verify minimum lengths
    for (Loops::const_iterator i = aligned.begin(); i != aligned.end(); ++i)
      TS_ASSERT(i->length() >= MIN_CHUNK_SZ);

    for (Loops::const_iterator i = unaligned.begin(); i != unaligned.end(); ++i)
      TS_ASSERT(i->length() >= MIN_CHUNK_SZ);

    // Verify contents
    TS_ASSERT_EQUALS(aligned.size(), 1);
    TS_ASSERT(aligned[1] == Loop(1, 53));

    TS_ASSERT_EQUALS(unaligned.size(), 1);
    TS_ASSERT(unaligned[1] == Loop(54, 63));

  }

  void test_limit_chunk_size() {
    const SequenceAlignment& alignment = alignments_[1];

    Loops aligned, unaligned;
    protocols::nonlocal::find_regions_with_minimum_size
        (alignment, MIN_CHUNK_SZ, &aligned, &unaligned);

    protocols::nonlocal::limit_chunk_size(MIN_CHUNK_SZ, MAX_CHUNK_SZ, &aligned);
    protocols::nonlocal::limit_chunk_size(MIN_CHUNK_SZ, MAX_CHUNK_SZ, &unaligned);

    TR << "Aligned: " << aligned << endl;
    TR << "Unaligned: " << unaligned << endl;

    for (Loops::const_iterator i = aligned.begin(); i != aligned.end(); ++i) {
      TS_ASSERT(i->length() >= MIN_CHUNK_SZ);
      TS_ASSERT(i->length() <= MAX_CHUNK_SZ);
    }

    for (Loops::const_iterator i = unaligned.begin(); i != unaligned.end(); ++i) {
      TS_ASSERT(i->length() >= MIN_CHUNK_SZ);
      TS_ASSERT(i->length() <= MAX_CHUNK_SZ);
    }
  }
};
}  // anonymous namespace
