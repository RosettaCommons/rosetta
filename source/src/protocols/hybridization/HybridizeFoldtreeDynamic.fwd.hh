// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/nonlocal/StarTreeBuilder.fwd.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_PROTOCOLS_HYBRIDIZATION_HYBRIDIZEFOLDTREEDYNAMIC_FWD_HH
#define INCLUDED_PROTOCOLS_HYBRIDIZATION_HYBRIDIZEFOLDTREEDYNAMIC_FWD_HH

#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace hybridization {

class HybridizeFoldtreeMover;
typedef utility::pointer::shared_ptr<HybridizeFoldtreeMover> HybridizeFoldtreeMoverOP;
typedef utility::pointer::shared_ptr<HybridizeFoldtreeMover const> HybridizeFoldtreeMoverCOP;

}  //  namespace hybridization
}  // namespace protocols

#endif  // PROTOCOLS_HYBRIDIZATION_HYBRIDIZEFOLDTREEMOVER_FWD_HH_
