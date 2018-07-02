// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file   core/pack/task/residue_selector/UnsatSelector.hh
/// @brief  A ResidueSelector that selects hydrogen bond acceptors or selectors that are not satisfied with h-bond
/// @author Parisa Hosseinzadeh (parisah@uw.edu)


#ifndef INCLUDED_protocols_hbnet_UnsatSelector_fwd_hh
#define INCLUDED_protocols_hbnet_UnsatSelector_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/vector1.fwd.hh>

// Forward
namespace protocols {
namespace hbnet {

class UnsatSelector;

typedef utility::pointer::shared_ptr< UnsatSelector > UnsatSelectorOP;
typedef utility::pointer::shared_ptr< UnsatSelector const > UnsatSelectorCOP;

} //hbnet
} //protocols


#endif //INCLUDED_protocols_hbnet_UnsatSelector_fwd_hh





