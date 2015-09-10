// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief  fwd headers for ns jumping
/// @author Oliver Lange


#ifndef INCLUDED_protocols_jumping_PairingLibrary_fwd_hh
#define INCLUDED_protocols_jumping_PairingLibrary_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>


namespace protocols {
namespace jumping {

// Forward
class BasePairingLibrary;
class PairingLibrary;
class StandardPairingLibrary;

// Types
typedef  utility::pointer::shared_ptr< BasePairingLibrary >  BasePairingLibraryOP;
typedef  utility::pointer::shared_ptr< BasePairingLibrary const >  BasePairingLibraryCOP;

// Types
typedef  utility::pointer::shared_ptr< PairingLibrary >  PairingLibraryOP;
typedef  utility::pointer::shared_ptr< PairingLibrary const >  PairingLibraryCOP;

// Types
typedef  utility::pointer::shared_ptr< StandardPairingLibrary >  StandardPairingLibraryOP;
typedef  utility::pointer::shared_ptr< StandardPairingLibrary const >  StandardPairingLibraryCOP;


} // namespace jumping
} // namespace protocols

#endif
