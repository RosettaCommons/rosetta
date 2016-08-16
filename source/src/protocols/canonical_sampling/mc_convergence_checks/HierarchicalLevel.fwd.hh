// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author
/// @author

#ifndef INCLUDED_protocols_canonical_sampling_mc_convergence_checks_HierarchicalLevel_fwd_hh
#define INCLUDED_protocols_canonical_sampling_mc_convergence_checks_HierarchicalLevel_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace canonical_sampling {
namespace mc_convergence_checks {

struct PoolData;
typedef utility::pointer::shared_ptr< PoolData > PoolDataOP;
typedef utility::pointer::shared_ptr< PoolData const > PoolDataCOP;

class HierarchicalLevel;
typedef utility::pointer::shared_ptr< HierarchicalLevel > HierarchicalLevelOP;
typedef utility::pointer::shared_ptr< HierarchicalLevel const > HierarchicalLevelCOP;


}
} // moves
} // protocols

#endif
