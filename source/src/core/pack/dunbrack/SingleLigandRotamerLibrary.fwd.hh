// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/dunbrack/SingleLigandRotamerLibrary.fwd.hh
/// @brief  SingleLigandRotamerLibrary forward declarations and typedefs
/// @author Ian W. Davis

#ifndef INCLUDED_core_pack_dunbrack_SingleLigandRotamerLibrary_fwd_hh
#define INCLUDED_core_pack_dunbrack_SingleLigandRotamerLibrary_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace core {
namespace pack {
namespace dunbrack {

class SingleLigandRotamerLibrary;
typedef utility::pointer::shared_ptr< SingleLigandRotamerLibrary > SingleLigandRotamerLibraryOP;
typedef utility::pointer::shared_ptr< SingleLigandRotamerLibrary const > SingleLigandRotamerLibraryCOP;
typedef utility::pointer::weak_ptr< SingleLigandRotamerLibrary > SingleLigandRotamerLibraryAP;
typedef utility::pointer::weak_ptr< SingleLigandRotamerLibrary const > SingleLigandRotamerLibraryCAP;

} // namespace dunbrack
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_pack_dunbrack_SingleLigandRotamerLibrary_fwd_hh
