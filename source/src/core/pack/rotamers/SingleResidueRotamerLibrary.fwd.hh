// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/rotamers/SingleResidueRotamerLibrary.fwd.hh
/// @brief  SingleResidueRotamerLibrary forward declarations and typedefs
/// @author Andrew Leaver-Fay
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

#ifndef INCLUDE_core_pack_rotamers_SingleResidueRotamerLibrary_fwd_hh
#define INCLUDE_core_pack_rotamers_SingleResidueRotamerLibrary_fwd_hh

#include <core/conformation/Residue.fwd.hh>

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace core {
namespace pack {
namespace rotamers {

class SingleResidueRotamerLibrary;

typedef utility::pointer::shared_ptr< SingleResidueRotamerLibrary > SingleResidueRotamerLibraryOP;
typedef utility::pointer::shared_ptr< SingleResidueRotamerLibrary const > SingleResidueRotamerLibraryCOP;
typedef utility::pointer::weak_ptr< SingleResidueRotamerLibrary > SingleResidueRotamerLibraryAP;
typedef utility::pointer::weak_ptr< SingleResidueRotamerLibrary const > SingleResidueRotamerLibraryCAP;

typedef utility::vector1< conformation::ResidueOP > RotamerVector;

} // namespace rotamers
} // namespace pack
} // namespace core

#endif
