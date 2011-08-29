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
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/nonlocal.OptionKeys.gen.hh>
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
  vector1<SequenceAlignment> alignments;
  Loops regions;

  void setUp() {
    using namespace basic::options;
    using namespace basic::options::OptionKeys;

    // Ensure consistent gap extension
    protocols_init();
    option[OptionKeys::nonlocal::gap_sampling_extension](3);

    vector1<string> filenames;
    filenames.push_back("protocols/nonlocal/alignment.filt");

    core::sequence::read_all_alignments("grishin",
                                        filenames,
                                        &alignments);
  }

  void test_find_regions() {
    const SequenceAlignment& alignment = alignments[1];

    Loops aligned, unaligned;
    protocols::nonlocal::find_regions(alignment,
                                      alignment.length(),
                                      &aligned,
                                      &unaligned);

    // Verify aligned
    TS_ASSERT_EQUALS(aligned.size(), 5);
    TS_ASSERT(aligned[1] == Loop(8, 80));
    TS_ASSERT(aligned[2] == Loop(88, 280));
    TS_ASSERT(aligned[3] == Loop(288, 305));
    TS_ASSERT(aligned[4] == Loop(315, 350));
    TS_ASSERT(aligned[5] == Loop(359, 367));

    // Verify unaligned
    TS_ASSERT_EQUALS(unaligned.size(), 6);
    TS_ASSERT(unaligned[1] == Loop(1, 7));
    TS_ASSERT(unaligned[2] == Loop(81, 87));
    TS_ASSERT(unaligned[3] == Loop(281, 287));
    TS_ASSERT(unaligned[4] == Loop(306, 314));
    TS_ASSERT(unaligned[5] == Loop(351, 358));
    TS_ASSERT(unaligned[6] == Loop(368, 390));
  }

  void test_limit_chunk_size() {
    const SequenceAlignment& alignment = alignments[1];

    Loops aligned, unaligned;
    protocols::nonlocal::find_regions(alignment,
                                      alignment.length(),
                                      &aligned,
                                      &unaligned);

    Loops combined;
    for (Loops::const_iterator i = aligned.begin(); i != aligned.end(); ++i)
      combined.add_loop(*i);
    for (Loops::const_iterator i = unaligned.begin(); i != unaligned.end(); ++i)
      combined.add_loop(*i);
    TS_ASSERT_EQUALS(combined.num_loop(), aligned.num_loop() + unaligned.num_loop());

    Size min_chunk_sz = 7;
    Size max_chunk_sz = 21;
    protocols::nonlocal::limit_chunk_size(min_chunk_sz,
                                          max_chunk_sz,
                                          &combined);

    TS_ASSERT(combined.num_loop() >= aligned.num_loop() + unaligned.num_loop());
    for (Loops::const_iterator i = combined.begin(); i != combined.end(); ++i) {
      const Loop& loop = *i;
      TS_ASSERT(loop.length() >= min_chunk_sz);
      TS_ASSERT(loop.length() <= max_chunk_sz);
    }
  }
};
}  // anonymous namespace
