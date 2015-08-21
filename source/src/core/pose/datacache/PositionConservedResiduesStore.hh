// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/pose/datacache/PositionConservedResiduesStore.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_CORE_POSE_DATACACHE_POSITIONCONSERVEDRESIDUESSTORE_HH
#define INCLUDED_CORE_POSE_DATACACHE_POSITIONCONSERVEDRESIDUESSTORE_HH

// Unit header
#include <core/pose/datacache/PositionConservedResiduesStore.fwd.hh>

// External headers
#include <boost/unordered/unordered_map.hpp>

// Utility headers
#include <basic/datacache/CacheableData.hh>

// Project headers
#include <core/types.hh>

namespace core {
namespace pose {
namespace datacache {

/// @class A cacheable container for storing/retrieving positional conservation
/// data about each residue. Once stored in a pose, data can be retrieved in the
/// following manner:
///
/// #include <core/pose/util.hh>
/// bool is_conserved = core::pose::is_position_conserved_residue(pose, residue)
//
class PositionConservedResiduesStore : public basic::datacache::CacheableData {
	typedef boost::unordered_map<core::Size, bool> ConservationMap;

public:
	/// @brief Default constructor
	PositionConservedResiduesStore();

	/// @brief Creates a copy of this instance
	basic::datacache::CacheableDataOP clone() const;

	/// @brief Returns true if <residue> is positionally conserved, false otherwise
	bool is_conserved(core::Size residue) const;

	/// @brief Updates the positional conservation of <residue>
	void set_conserved(core::Size residue, bool conserved);

private:
	/// @brief Associates real-valued structural conservation scores with residues
	mutable ConservationMap conservation_;
};

}  // namespace datacache
}  // namespace pose
}  // namespace core

#endif  // CORE_POSE_DATACACHE_POSITIONCONSERVEDRESIDUESSTORE_HH_
