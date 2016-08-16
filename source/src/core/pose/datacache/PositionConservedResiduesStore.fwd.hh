// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pose/datacache/PositionConservedResiduesStore.fwd.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_CORE_POSE_DATACACHE_POSITIONCONSERVEDRESIDUESSTORE_FWD_HH
#define INCLUDED_CORE_POSE_DATACACHE_POSITIONCONSERVEDRESIDUESSTORE_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pose {
namespace datacache {

class PositionConservedResiduesStore;
typedef utility::pointer::shared_ptr<PositionConservedResiduesStore> PositionConservedResiduesStoreOP;
typedef utility::pointer::shared_ptr<PositionConservedResiduesStore const> PositionConservedResiduesStoreCOP;

}  // namespace datacache
}  // namespace pose
}  // namespace core

#endif  // CORE_POSE_DATACACHE_POSITIONCONSERVEDRESIDUESSTORE_FWD_HH_
