// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/palette/PackerPaletteFactory.fwd.hh
/// @brief  Forward declarations for factory class for creating instances of PackerPalettes (e.g. for RosettaScripts).
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_core_pack_palette_PackerPaletteFactorywd_hh
#define INCLUDED_core_pack_palette_PackerPaletteFactorywd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace palette {

class PackerPaletteFactory;

typedef utility::pointer::shared_ptr< PackerPaletteFactory > PackerPaletteFactoryOP;
typedef utility::pointer::shared_ptr< PackerPaletteFactory const > PackerPaletteFactoryCOP;

} //namespace palette
} //namespace pack
} //namespace core

#endif
