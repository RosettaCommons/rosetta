// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/NLFragmentGroup.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef PROTOCOLS_NONLOCAL_NLFRAGMENTGROUP_HH_
#define PROTOCOLS_NONLOCAL_NLFRAGMENTGROUP_HH_

// Package headers
#include <protocols/nonlocal/NLFragment.hh>

// Project headers
#include <core/types.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace nonlocal {

/// @class A data structure containing information about a consecutive
/// sequence of residues, each element of which having type NLFragment.
class NLFragmentGroup {
  typedef core::Size Size;

 public:
  /// @brief Constructs a new, and initially empty, NLFragmentGroup
  NLFragmentGroup();

  /// @brief Returns the number of entries
  Size num_entries() const;

  /// @brief Retrieves the ith entry
  const NLFragment& entries(Size i) const;

  /// @brief Adds the specified entry
  void add_entry(const NLFragment& e);

  /// @brief Removes all entries
  void clear_entries();

  /// @brief Returns the lowest-numbered position in this fragment
  Size start() const;

  /// @brief Returns the highest-numbered position in this fragment
  Size stop() const;

  /// @brief Orders this fragment's entries by their <position> field. Must be
  /// called prior to start() or stop(). This is done automatically in
  /// NonlocalAbinitioReader. The warning mainly applies to those who wish to
  /// construct NLGroupings programmatically.
  void sort();

  // Operators
  bool operator==(const NLFragmentGroup& other) const;
  bool operator!=(const NLFragmentGroup& other) const;
  bool operator< (const NLFragmentGroup& other) const;

 private:
  utility::vector1<NLFragment> entries_;
};

}  // namespace nonlocal
}  // namespace protocols

#endif  // PROTOCOLS_NONLOCAL_NLFRAGMENTGROUP_HH_
