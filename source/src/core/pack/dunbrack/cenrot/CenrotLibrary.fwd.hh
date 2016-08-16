// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/dunbrack/cenrot/CenrotLibrary.fwd.hh
/// @brief  CenrotLibrary forward declarations and typedefs
/// @author


#ifndef INCLUDED_core_pack_dunbrack_cenrot_CenrotLibrary_fwd_hh
#define INCLUDED_core_pack_dunbrack_cenrot_CenrotLibrary_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace core {
namespace pack {
namespace dunbrack {
namespace cenrot {

class CenrotLibrary;
typedef utility::pointer::shared_ptr< CenrotLibrary > CenrotLibraryOP;
typedef utility::pointer::shared_ptr< CenrotLibrary const > CenrotLibraryCOP;

} // cenrot
} // dunbrack
} // pack
} // core

#endif
