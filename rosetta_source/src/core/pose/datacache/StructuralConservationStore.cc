// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/pose/datacache/StructuralConservationStore.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <core/pose/datacache/StructuralConservationStore.hh>

// Project headers
#include <core/types.hh>

namespace core {
namespace pose {
namespace datacache {

StructuralConservationStore::StructuralConservationStore() {}

basic::datacache::CacheableDataOP StructuralConservationStore::clone() const {
  return new StructuralConservationStore(*this);
}

core::Real StructuralConservationStore::conservation_score(core::Size residue) const {
  return (conservation_.find(residue) != conservation_.end()) ? conservation_[residue] : -1;
}

void StructuralConservationStore::set_conservation_score(core::Size residue, core::Real score) {
  conservation_[residue] = score;
}

}  // namespace datacache
}  // namespace pose
}  // namespace core
