// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/palette/PackerPalette.fwd.hh
/// @brief  PackerPalette: a class for storing the set of ResidueTypes
/// and VariantTypes that the packer uses by default, in the absence of any
/// TaskOperations that limit the set actually used.
/// @details The PackerPalette says, "Here are the types that you're
/// allowed to use, and which are on in the absence of TaskOperations."
/// TaskOperations then prune this, turning OFF types that have been
/// enabled.  This allows users to turn on noncanonicals for design, and
/// then use TaskOperations with the same commutativity rules (turning OFF
/// types only) that are used for canonicals, making mixed design with
/// canonicals and noncanonicals much easier.\nThis was implemented as
/// part of the 2016 Chemical XRW (eXtreme Rosetta Workshop).
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu).

#ifndef INCLUDED_core_pack_palette_PackerPalette_fwd_hh
#define INCLUDED_core_pack_palette_PackerPalette_fwd_hh

//Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace palette {

class PackerPalette;

typedef utility::pointer::shared_ptr< PackerPalette > PackerPaletteOP;
typedef utility::pointer::shared_ptr< PackerPalette const > PackerPaletteCOP;


} // namespace palette
} // namespace pack
} // namespace core

#endif //INCLUDED_core_pack_palette_PackerPalette_fwd_hh
