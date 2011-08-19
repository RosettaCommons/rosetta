// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/pose/datacache/StructuralConservationStore.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef CORE_POSE_DATACACHE_STRUCTURALCONSERVATIONSTORE_HH_
#define CORE_POSE_DATACACHE_STRUCTURALCONSERVATIONSTORE_HH_

// Unit header
#include <core/pose/datacache/StructuralConservationStore.fwd.hh>

// External headers
#include <boost/unordered/unordered_map.hpp>

// Utility headers
#include <basic/datacache/CacheableData.hh>

// Project headers
#include <core/types.hh>

namespace core {
namespace pose {
namespace datacache {

/// @class A cacheable container for storing/retrieving per-residue structural
/// conservation data. Once stored in a pose, data can be retrieved in the
/// following manner:
///
/// #include <core/pose/util.hh>
/// core::Real cons_score = core::pose::structural_conservation(pose, residue)
//
class StructuralConservationStore : public basic::datacache::CacheableData {
  typedef boost::unordered_map<core::Size, core::Real> ConservationMap;

 public:
  /// @brief Default constructor
  StructuralConservationStore();

  /// @brief Creates a copy of this instance
  basic::datacache::CacheableDataOP clone() const;

  /// @brief Returns the structural conservation score associated with <residue>
  /// if it exists, -1 otherwise.
  core::Real conservation_score(core::Size residue) const;

  /// @brief Updates the structural conservation score associated with <residue>
  void set_conservation_score(core::Size residue, core::Real score);

 private:
  /// @brief Associates real-valued structural conservation scores with residues
  mutable ConservationMap conservation_;
};

}  // namespace datacache
}  // namespace pose
}  // namespace core

#endif  // CORE_POSE_DATACACHE_STRUCTURALCONSERVATIONSTORE_HH_
