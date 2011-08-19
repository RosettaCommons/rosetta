// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/NonlocalAbinitio.cxxtest.hh
/// @brief test suite for protocols/nonlocal/NonlocalAbinitio
/// @author Christopher Miles (cmiles@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Utility headers
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/nonlocal.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
#include <protocols/nonlocal/NLGrouping.hh>
#include <protocols/nonlocal/NonlocalAbinitio.hh>
#include <protocols/nonlocal/NonlocalAbinitioReader.hh>

namespace {
using protocols::nonlocal::NLGrouping;
using protocols::nonlocal::NonlocalAbinitio;
using protocols::nonlocal::NonlocalAbinitioReader;

class NonlocalAbinitioTest : public CxxTest::TestSuite {
 public:
    void setUp() {
        core_init();
    }
       
	void test_constructors() {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		// Specify options as if done from the command line.
		// The fasta file isn't actually read, so the file
		// does not have to exist.
		option[OptionKeys::nonlocal::moves]("protocols/nonlocal/test.pairings");
		option[OptionKeys::in::file::fasta]("protocols/nonlocal/does_not_exist");
		option[OptionKeys::in::file::frag3]("protocols/nonlocal/mfr_aa2GB3_03_05.200_v1_3");
		option[OptionKeys::in::file::frag9]("protocols/nonlocal/mfr_aa2GB3_09_05.200_v1_3");

		// constructed from options
		NonlocalAbinitio m1;

		// constructed explicitly
		utility::vector1<NLGrouping> groupings;
    NonlocalAbinitioReader::read("protocols/nonlocal/test.pairings", &groupings);
		NonlocalAbinitio m2(groupings);

		TS_ASSERT_EQUALS(m1.groupings(), m2.groupings());
	}
};
}  // anonymous namespace
