// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/NLGrouping.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef PROTOCOLS_NONLOCAL_NLGROUPING_HH_
#define PROTOCOLS_NONLOCAL_NLGROUPING_HH_

// C/C++ headers
#include <string>

// Utility headers
#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
#include <utility/vector1.hh>

// Package headers
#include <protocols/nonlocal/NLFragmentGroup.hh>
#include <protocols/nonlocal/Region.hh>

namespace protocols {
namespace nonlocal {

/// @class A data structure containing information about a set of known or
/// predicted interactions between groups of residues within a given protein.
/// Each element has type NLFragmentGroup.
class NLGrouping {
  typedef core::Size Size;

 public:
  /// @brief Constructs a new, and initially empty, NLGrouping
  NLGrouping();

  /// @brief Returns the number of groups
  Size num_groups() const;

  /// @brief Retrieves the ith group
  const NLFragmentGroup& groups(Size i) const;

  /// @brief Adds the specified group
  void add_group(const NLFragmentGroup& g);

  /// @brief Removes all groups
  void clear_groups();

  /// @brief Orders the groups by starting position. Must be called prior to
  /// serialize(). This is done automatically in NonlocalAbinitioReader. The
  /// warning mainly applies to those who wish to construct NLGroupings
  /// programmatically.
  void sort();

	/// @brief Populate <regions> with details about the start/stop positions
	/// of the rigid chunks in this NLGrouping.
	void as_regions(utility::vector1<Region>* regions) const;

	/// @brief Computes the center of mass, storing the result in <center>.
	void center_of_mass(numeric::xyzVector<core::Real>* center) const;

  // Serialization
  std::string provenance() const;

  // Operators
  bool operator==(const NLGrouping& other) const;
  bool operator!=(const NLGrouping& other) const;
  bool operator <(const NLGrouping& other) const;

 private:
  utility::vector1<NLFragmentGroup> groups_;
};

}  // namespace nonlocal
}  // namespace protocols

#endif  // PROTOCOLS_NONLOCAL_NLGROUPING_HH_
