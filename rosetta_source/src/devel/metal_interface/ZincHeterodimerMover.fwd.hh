// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/*/ZincHeterodimerMover.fwd.hh
/// @brief  ZincHeterodimerMover protocol-mover forward declarations header
/// @author Steven Lewis (smlewi@unc.edu)


#ifndef INCLUDED_devel_metal_interface_ZincHeterodimerMover_fwd_hh
#define INCLUDED_devel_metal_interface_ZincHeterodimerMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace devel{
namespace metal_interface{

//Forwards and OP typedefs
class ZincHeterodimerMover;
typedef utility::pointer::owning_ptr< ZincHeterodimerMover > ZincHeterodimerMoverOP;
typedef utility::pointer::owning_ptr< ZincHeterodimerMover const > ZincHeterodimerMoverCOP;

}//MetalInterface
}//devel

#endif //INCLUDED_devel_MetalInterface_ZincHeterodimerMover_FWD_HH
