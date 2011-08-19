// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/NonlocalAbinitioReader.cxxtest.hh
/// @brief test suite for protocols/nonlocal/NonlocalAbinitioReader
/// @author Christopher Miles (cmiles@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>

// Utility headers
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
#include <protocols/nonlocal/NLFragment.hh>
#include <protocols/nonlocal/NLFragmentGroup.hh>
#include <protocols/nonlocal/NLGrouping.hh>
#include <protocols/nonlocal/NonlocalAbinitioReader.hh>

// C/C++ headers
#include <string>

#define MARGIN 0.0001
#define UNUSED 0

namespace {
using core::Real;
using core::Size;
using protocols::nonlocal::NLFragment;
using protocols::nonlocal::NLFragmentGroup;
using protocols::nonlocal::NonlocalAbinitioReader;
using protocols::nonlocal::NLGrouping;
using std::string;

class NonlocalAbinitioReaderTest : public CxxTest::TestSuite {
 public:
  void test_reader() {
    utility::vector1<NLGrouping> groupings;
    NonlocalAbinitioReader::read("protocols/nonlocal/test.pairings", &groupings);

		// correct number of NonlocalGroupings?
    TS_ASSERT_EQUALS(groupings.size(), 1);
    const NLGrouping& grouping = groupings[1];

		// correct number of FragmentGroups?
		TS_ASSERT_EQUALS(grouping.num_groups(), 4);

    // correct number of FragmentGroupEntries in each FragmentGroup?
    Size entries_per_group[] = { UNUSED, 5, 6, 15, 5 };
		for (Size i = 1; i <= grouping.num_groups(); ++i)
      TS_ASSERT_EQUALS(grouping.groups(i).num_entries(), entries_per_group[i]);

		// check the contents of the second FragmentGroup
		NLFragmentGroup g2 = grouping.groups(2);
		check_entry(g2.entries(1), 49, -143.132, 128.827, 176.942, 94.852, 53.324, -6.631);
		check_entry(g2.entries(2), 50, -118.006, 112.700, -171.815, 93.025, 50.630, -4.667);
		check_entry(g2.entries(3), 51, -105.538, 117.641, 178.295, 92.132, 51.602, -1.106);
		check_entry(g2.entries(4), 52, -102.343, 117.860, 179.175, 91.170, 48.933, 1.425);
		check_entry(g2.entries(5), 53, -107.322, 91.600, 179.354, 89.445, 50.179, 4.575);
		check_entry(g2.entries(6), 54, -71.162, 135.072, 176.466, 89.341, 47.329, 7.090);
  }

	/// @brief Checks the start/stop positions of a FragmentGroup
	void test_start_stop() {
		utility::vector1<NLGrouping> groupings;
		NonlocalAbinitioReader::read("protocols/nonlocal/test.pairings", &groupings);
		TS_ASSERT_EQUALS(groupings.size(), 1);

		const NLGrouping& grouping = groupings[1];
		Size start_positions[] = { UNUSED, 2, 49, 66, 84 };
		Size stop_positions [] = { UNUSED, 6, 54, 80, 88 };

		for (Size i = 1; i <= grouping.num_groups(); ++i) {
			const NLFragmentGroup& group = grouping.groups(i);
			TS_ASSERT_EQUALS(group.start(), start_positions[i]);
			TS_ASSERT_EQUALS(group.stop(), stop_positions[i]);
		}
  }

	// @brief Compare result of pairings files that differ only in use of
	// syntactic sugar.
	void test_check_syntactic_sugar() {
		utility::vector1<NLGrouping> groupings_1, groupings_2;
		NonlocalAbinitioReader::read("protocols/nonlocal/test.pairings", &groupings_1);
		NonlocalAbinitioReader::read("protocols/nonlocal/test2.pairings", &groupings_2);

		TS_ASSERT_EQUALS(groupings_1.size(), groupings_2.size());
		for (Size i = 1; i <= groupings_1.size(); ++i) {
			const NLGrouping& grouping_1 = groupings_1[i];
			const NLGrouping& grouping_2 = groupings_2[i];
			TS_ASSERT_EQUALS(grouping_1.num_groups(), grouping_2.num_groups());

			for (Size j = 1; j <= grouping_1.num_groups(); ++j) {
				const NLFragmentGroup& group_1 = grouping_1.groups(j);
				const NLFragmentGroup& group_2 = grouping_2.groups(j);
				TS_ASSERT_EQUALS(group_1.num_entries(), group_2.num_entries());
				TS_ASSERT_EQUALS(group_1.start(), group_2.start());
				TS_ASSERT_EQUALS(group_1.stop(), group_2.stop());

				for (Size k = 1; k <= group_1.num_entries(); ++k) {
					const NLFragment& entry_1 = group_1.entries(k);
					const NLFragment& entry_2 = group_2.entries(k);
					TS_ASSERT_EQUALS(entry_1.position(), entry_2.position());
					TS_ASSERT_EQUALS(entry_1.phi(), entry_2.phi());
					TS_ASSERT_EQUALS(entry_1.psi(), entry_2.psi());
					TS_ASSERT_EQUALS(entry_1.omega(), entry_2.omega());
					TS_ASSERT_EQUALS(entry_1.x(), entry_2.x());
					TS_ASSERT_EQUALS(entry_1.y(), entry_2.y());
					TS_ASSERT_EQUALS(entry_1.z(), entry_2.z());
				}
			}
		}
	}

	void test_provenance() {
		utility::vector1<NLGrouping> groupings;
		NonlocalAbinitioReader::read("protocols/nonlocal/test.pairings", &groupings);
		TS_ASSERT_EQUALS(groupings.size(), 1);

		const NLGrouping& grouping = groupings[1];
		TS_ASSERT_EQUALS(grouping.provenance(), "{ 2-6 49-54 66-80 84-88 }");
	}

	void test_ordering() {
		utility::vector1<NLGrouping> groupings;
		NonlocalAbinitioReader::read("protocols/nonlocal/unordered.pairings", &groupings);
		TS_ASSERT_EQUALS(groupings.size(), 1);

		const NLGrouping& grouping = groupings[1];
		TS_ASSERT_EQUALS(grouping.provenance(), "{ 2-6 49-54 66-80 84-88 }");
	}

	/// @brief Utility method for checking the contents of a FragmentGroupEntry
	void check_entry(const NLFragment& e, Size pos, Real phi, Real psi, Real omega, Real x, Real y, Real z) {
		TS_ASSERT_EQUALS(e.position(), pos);
		TS_ASSERT_DELTA(e.phi(), phi, MARGIN);
		TS_ASSERT_DELTA(e.psi(), psi, MARGIN);
		TS_ASSERT_DELTA(e.omega(), omega, MARGIN);
		TS_ASSERT_DELTA(e.x(), x, MARGIN);
		TS_ASSERT_DELTA(e.y(), y, MARGIN);
		TS_ASSERT_DELTA(e.z(), z, MARGIN);
	}
};
}  // anonymous namespace
