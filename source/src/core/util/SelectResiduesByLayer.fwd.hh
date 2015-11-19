// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/util/SelectResiduesByLayer.fwd.hh.
/// @brief  Forward declarations for a class permitting selection of residues by layer (burial).
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )
/// @author Gabe Rocklin (sidechain neighbour selection)
/// @author Vikram K. Mulligan (vmullig@uw.edu -- moving this class to core and refactoring for noncanonicals)

#ifndef INCLUDED_core_util_SelectResiduesByLayer_fwd_hh
#define INCLUDED_core_util_SelectResiduesByLayer_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace util {

class SelectResiduesByLayer;
typedef utility::pointer::shared_ptr< SelectResiduesByLayer > SelectResiduesByLayerOP;
typedef utility::pointer::shared_ptr< SelectResiduesByLayer const > SelectResiduesByLayerCOP;

}
}

#endif
