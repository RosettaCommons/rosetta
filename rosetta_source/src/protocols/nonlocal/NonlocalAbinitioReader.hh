// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/NonlocalAbinitioReader.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef PROTOCOLS_NONLOCAL_NONLOCALABINITIOREADER_HH_
#define PROTOCOLS_NONLOCAL_NONLOCALABINITIOREADER_HH_

// Package headers
#include <protocols/nonlocal/NLFragmentGroup.hh>
#include <protocols/nonlocal/NLGrouping.hh>

// Utility headers
#include <utility/vector1.hh>

// C/C++ headers
#include <string>

namespace protocols {
namespace nonlocal {

class NonlocalAbinitioReader {
  typedef utility::vector1<NLGrouping> NonlocalGroupings;

 public:
  /// @brief Reads the contents of <filename> into <groups>. NonlocalGrouping's
  /// are contained within '{' and '}'. Each NonlocalGrouping contains 0 or more
  /// FragmentGroups, which are delimited by empty lines. FragmentGroups consist
  /// of 0 or more FragmentGroupEntries. The data within a FragmentGroupEntry is
  /// comma-delimited.
  ///
  /// See test/protocols/nonlocal/test.pairings for a concrete example.
  static void read(const std::string& filename, NonlocalGroupings* groups);

 private:
  static void emit_group(NLFragmentGroup* group, NLGrouping* grouping);
};

}  // namespace nonlocal
}  // namespace protocols

#endif  // PROTOCOLS_NONLOCAL_NONLOCALABINITIOREADER_HH_
