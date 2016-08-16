// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/metal_interface/ZincHeterodimerMover.fwd.hh
/// @brief  ZincHeterodimerMover protocol-mover forward declarations header
/// @author Steven Lewis (smlewi@gmail.com)


#ifndef INCLUDED_protocols_metal_interface_HeterodimerMover_fwd_hh
#define INCLUDED_protocols_metal_interface_HeterodimerMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace metal_interface {

//Forwards and OP typedefs
class ZincHeterodimerMover;
typedef utility::pointer::shared_ptr< ZincHeterodimerMover > ZincHeterodimerMoverOP;
typedef utility::pointer::shared_ptr< ZincHeterodimerMover const > ZincHeterodimerMoverCOP;

}//metal_interface
}//protocols

#endif //INCLUDED_protocols_metal_interface_ZincHeterodimerMover_FWD_HH
