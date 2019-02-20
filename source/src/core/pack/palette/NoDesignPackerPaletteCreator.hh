// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/palette/NoDesignPackerPaletteCreator.hh
/// @brief  Declaration of the creator for NoDesignPackerPalette.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_core_pack_palette_NoDesignPackerPaletteCreator_hh
#define INCLUDED_core_pack_palette_NoDesignPackerPaletteCreator_hh

// Unit headers
#include <core/pack/palette/PackerPaletteCreator.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <string>


namespace core {
namespace pack {
namespace palette {

/// @brief The PackerPaletteCreator class's responsibilities are to create
/// on demand a new instance of a PackerPalette class.
/// @details The PackerPaletteCreator must register itself with the PackerPaletteFactory
/// at load time (before main() begins) so that the PackerPaletteFactory is ready
/// to start creating PackerPalettes by the time any protocol
/// requests one.
class NoDesignPackerPaletteCreator : public PackerPaletteCreator
{
public:
	virtual ~NoDesignPackerPaletteCreator() {}

	/// @brief Instantiate a new PackerPalette
	PackerPaletteOP create_packer_palette() const override;
	std::string keyname() const override { return "NoDesignPackerPalette"; }

	/// @brief Describe the allowed XML options for a particular PackerPalette subclass.
	/// @details Pure virtual.  Must be implemented by subclasses.
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

}
}
}

#endif
