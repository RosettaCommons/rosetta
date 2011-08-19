// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/pose/datacache/StructuralConservationStore.fwd.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef CORE_POSE_DATACACHE_STRUCTURALCONSERVATIONSTORE_FWD_HH_
#define CORE_POSE_DATACACHE_STRUCTURALCONSERVATIONSTORE_FWD_HH_

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pose {
namespace datacache {

class StructuralConservationStore;
typedef utility::pointer::owning_ptr<StructuralConservationStore> StructuralConservationStoreOP;
typedef utility::pointer::owning_ptr<StructuralConservationStore const> StructuralConservationStoreCOP;

}  // namespace datacache
}  // namespace pose
}  // namespace core

#endif  // CORE_POSE_DATACACHE_STRUCTURALCONSERVATIONSTORE_FWD_HH_
