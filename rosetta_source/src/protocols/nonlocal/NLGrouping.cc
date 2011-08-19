// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/NLGrouping.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <protocols/nonlocal/NLGrouping.hh>

// C/C++ headers
#include <cassert>
#include <iterator>
#include <sstream>
#include <string>

// External headers
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>

// Utility headers
#include <numeric/xyzVector.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>

// Package headers
#include <protocols/nonlocal/NLFragment.hh>
#include <protocols/nonlocal/NLFragmentGroup.hh>

namespace protocols {
namespace nonlocal {

typedef core::Size Size;

NLGrouping::NLGrouping() {}

Size NLGrouping::num_groups() const {
  return groups_.size();
}

const NLFragmentGroup& NLGrouping::groups(Size i) const {
  return groups_[i];
}

void NLGrouping::add_group(const NLFragmentGroup& g) {
  groups_.push_back(g);
}

void NLGrouping::clear_groups() {
  groups_.clear();
}

bool NLGrouping::operator==(const NLGrouping& other) const {
  if (num_groups() != other.num_groups())
    return false;

  for (core::Size i = 1; i <= num_groups(); ++i) {
    if (groups(i) != other.groups(i))
      return false;
  }

  return true;
}

bool NLGrouping::operator!=(const NLGrouping& other) const {
  return !(*this == other);
}

bool NLGrouping::operator<(const NLGrouping& other) const {
  assert(num_groups() > 0);
  assert(other.num_groups() > 0);
  return groups(1).start() < other.groups(1).start();
}

/// @detail Recursively sorts the contents of <groups_>
void NLGrouping::sort() {
  utility::vector1<NLFragmentGroup>::iterator i;
  for (i = groups_.begin(); i != groups_.end(); ++i)
    i->sort();

  std::sort(groups_.begin(), groups_.end());
}

void NLGrouping::as_regions(utility::vector1<Region>* regions) const {
  for (core::Size i = 1; i <= num_groups(); ++i) {
    const NLFragmentGroup& group = groups(i);
    regions->push_back(Region(group.start(), group.stop()));
  }
}

std::string NLGrouping::provenance() const {
  std::stringstream ss;

  ss << "{ ";
  for (core::Size i = 1; i <= num_groups(); ++i)
    ss << groups(i).start() << "-" << groups(i).stop() << " ";
  ss <<  "}";

  return ss.str();
}

void NLGrouping::center_of_mass(numeric::xyzVector<core::Real>* center) const {
  using core::Real;
  using core::Size;

  Real mass_x = 0, mass_y = 0, mass_z = 0;
  Size num_residues = 0;

  for (Size i = 1; i <= num_groups(); ++i) {
    const NLFragmentGroup& group = groups(i);
    for (Size j = 1; j <= group.num_entries(); ++j) {
      const NLFragment& fragment = group.entries(j);
      mass_x += fragment.x();
      mass_y += fragment.y();
      mass_z += fragment.z();
      ++num_residues;
    }
  }

  center->x(mass_x / num_residues);
  center->y(mass_y / num_residues);
  center->z(mass_z / num_residues);
}

}  // namespace nonlocal
}  // namespace protocols
