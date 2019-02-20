// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief  Declaration of the base class for PackerPalette factory registration and creation
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_core_pack_palette_PackerPaletteRegistrator_hh
#define INCLUDED_core_pack_palette_PackerPaletteRegistrator_hh

// Package headers
#include <core/pack/palette/PackerPaletteFactory.hh>
#include <utility/factory/WidgetRegistrator.hh>

namespace core {
namespace pack {
namespace palette {

/// @brief This templated class will register an instance of an
/// PackerPalettteCreator (class T) with the PackerPaletteFactory.  It will ensure
/// that no PackerPaletteCreator is registered twice, and, centralizes
/// this registration logic so that thread safety issues can be handled in
/// one place.
template < class T >
class PackerPaletteRegistrator : public utility::factory::WidgetRegistrator< PackerPaletteFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< PackerPaletteFactory, T > parent;
public:
	PackerPaletteRegistrator() : parent() {}
};

}
}
}

#endif
