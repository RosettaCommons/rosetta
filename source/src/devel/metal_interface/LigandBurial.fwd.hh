// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   devel/metal_interface/LigandBurial.fwd.hh
/// @brief  LigandBurial protocol-mover forward declarations header
/// @author Bryan Der


#ifndef INCLUDED_devel_metal_interface_LigandBurial_FWD_HH
#define INCLUDED_devel_metal_interface_LigandBurial_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace devel{
namespace metal_interface{

//Forwards and OP typedefs
class LigandBurial;
typedef utility::pointer::owning_ptr< LigandBurial > LigandBurialOP;
typedef utility::pointer::owning_ptr< LigandBurial const > LigandBurialCOP;

}//metal_interface
}//devel

#endif //INCLUDED_devel_metal_interface_LigandBurial_FWD_HH
