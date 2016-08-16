// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/rotamers/SingleResiduePeptoidLibrary.fwd.hh
/// @brief SingleResiduePeptoidLibrary forward declarations and typedefs
/// @author

#ifndef INCLUDED_core_pack_rotamers_SingleResiduePeptoidLibrary_fwd_hh
#define INCLUDED_core_pack_rotamers_SingleResiduePeptoidLibrary_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace core {
namespace pack {
namespace rotamers {

class SingleResiduePeptoidLibrary;
typedef utility::pointer::shared_ptr< SingleResiduePeptoidLibrary > SingleResiduePeptoidLibraryOP;
typedef utility::pointer::weak_ptr< SingleResiduePeptoidLibrary > SingleResiduePeptoidLibraryAP;
typedef utility::pointer::shared_ptr< SingleResiduePeptoidLibrary const > SingleResiduePeptoidLibraryCOP;
typedef utility::pointer::weak_ptr< SingleResiduePeptoidLibrary const > SingleResiduePeptoidLibraryCAP;

} // rotamers
} // pack
} // core

#endif // INCLUDED_core_pack_rotamers_SingleResiduePeptoidLibrary_HH

