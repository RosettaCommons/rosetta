// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/palette/CustomBaseTypePackerPalette.fwd.hh
/// @brief  CustomBaseTypePackerPalette: a PackerPalette that allows a user to define additional
/// ResidueTypes with which to design (but not additional VariantTypes, at this time).
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu).

#ifndef INCLUDED_core_pack_palette_CustomBaseTypePackerPalette_fwd_hh
#define INCLUDED_core_pack_palette_CustomBaseTypePackerPalette_fwd_hh

//Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace palette {

class CustomBaseTypePackerPalette;

typedef utility::pointer::shared_ptr< CustomBaseTypePackerPalette > CustomBaseTypePackerPaletteOP;
typedef utility::pointer::shared_ptr< CustomBaseTypePackerPalette const > CustomBaseTypePackerPaletteCOP;


} // namespace palette
} // namespace pack
} // namespace core

#endif //INCLUDED_core_pack_palette_PackerPalette_fwd_hh
