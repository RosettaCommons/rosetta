// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/NLFragmentGroup.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <protocols/nonlocal/NLFragmentGroup.hh>

// Project headers
#include <core/types.hh>

// C/C++ headers
#include <algorithm>

namespace protocols {
namespace nonlocal {

typedef core::Size Size;

NLFragmentGroup::NLFragmentGroup() {}

Size NLFragmentGroup::num_entries() const {
  return entries_.size();
}

const NLFragment& NLFragmentGroup::entries(Size i) const {
  return entries_[i];
}

void NLFragmentGroup::add_entry(const NLFragment& e) {
  entries_.push_back(e);
}

void NLFragmentGroup::clear_entries() {
  entries_.clear();
}

Size NLFragmentGroup::start() const {
  return entries(1).position();
}

Size NLFragmentGroup::stop() const {
  return entries(num_entries()).position();
}

void NLFragmentGroup::sort() {
  std::sort(entries_.begin(), entries_.end());
}

bool NLFragmentGroup::operator==(const NLFragmentGroup& other) const {
  if (num_entries() != other.num_entries())
    return false;

  for (Size i = 1; i <= num_entries(); ++i) {
    if (entries(i) != other.entries(i))
      return false;
  }

  return true;
}

bool NLFragmentGroup::operator!=(const NLFragmentGroup& other) const {
  return !(*this == other);
}

bool NLFragmentGroup::operator<(const NLFragmentGroup& other) const {
  return start() < other.start();
}

}  // namespace nonlocal
}  // namespace protocols
