// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/palette/NoDesignPackerPalette.fwd.hh
/// @brief  NoDesignPackerPalette: a PackerPalette that sets up absolutely no design residues.
/// @details Not necessary for repacking (the RestrictToRepacking task operation produces the same effect), but handy
/// for efficiency when you know that you're not doing any design.  (There's no point setting up a list of ResidueTypes)
/// only to prune them all away.
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu).

#ifndef INCLUDED_core_pack_palette_NoDesignPackerPalette_fwd_hh
#define INCLUDED_core_pack_palette_NoDesignPackerPalette_fwd_hh

//Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace palette {

class NoDesignPackerPalette;

typedef utility::pointer::shared_ptr< NoDesignPackerPalette > NoDesignPackerPaletteOP;
typedef utility::pointer::shared_ptr< NoDesignPackerPalette const > NoDesignPackerPaletteCOP;


} // namespace palette
} // namespace pack
} // namespace core

#endif //INCLUDED_core_pack_palette_PackerPalette_fwd_hh
