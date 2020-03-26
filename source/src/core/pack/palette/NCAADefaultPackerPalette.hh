// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/palette/NCAADefaultPackerPalette.hh
/// @brief  NCAADefaultPackerPalette: a PackerPalette that just sets up
/// default Packer behaviour (design with the canonical 20 amino acids and equivalent
/// sets for ANY backbone.)
// @author Andy Watkins (amw579@stanford.edu).

#ifndef INCLUDED_core_pack_palette_NCAADefaultPackerPalette_hh
#define INCLUDED_core_pack_palette_NCAADefaultPackerPalette_hh

// Unit Headers
#include <core/pack/palette/NCAADefaultPackerPalette.fwd.hh>

// Package Headers

// Project Headers
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/pack/palette/PackerPalette.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/VirtualBase.hh>
#include <utility/vector1.fwd.hh>

// STL Headers
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace palette {

/// @brief  The NCAADefaultPackerPalette class gives instructions to the packer about
/// the set of ResidueTypes and VariantTypes to use by default, in the
/// absence of any TaskOperations that prune the list.  Specifically, the NCAADefaultPackerPalette
/// says, "By default, use the twenty canonical amino acids and whatever is present at a given position -- and
/// nothing else."
class NCAADefaultPackerPalette : public PackerPalette
{

	/// @brief The parent class (PackerPalette).
	///
	typedef PackerPalette parent;

public:

	/// @brief NCAADefault constructor.
	///
	NCAADefaultPackerPalette();

	/// @brief Destructor.
	///
	~NCAADefaultPackerPalette() override;

	/// @brief Clone operator: make a copy and return an owning pointer to the copy.
	/// @details Derived classes MUST implement this.
	PackerPaletteOP clone() const override;

	/// @brief Function to parse XML tags, implemented by derived classes.
	/// @brief Failure to implement this results in a warning message, but does not prevent compilation or
	/// program execution.
	void
	parse_my_tag(
		utility::tag::TagCOP const &tag,
		basic::datacache::DataMap const &datamap
	) override;

	/// @brief Provide information about the XML options available for this PackerPalette.
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief Get the name of this object ("NCAADefaultPackerPalette").
	std::string const & name() const override;

public: //Functions for the CitationManager:

	/// @brief Returns true (this packer palette is unpublished).
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
	bool packer_palette_is_unpublished() const override;

	/// @brief Returns Andy Watkins as the author.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
	utility::vector1< basic::citation_manager::UnpublishedModuleInfoCOP > provide_authorship_info_for_unpublished() const override;

protected:

	/// @brief Generate a list of possible base residue types
	/// @param [in] restypeset The ResidueTypeSet to use as a reference for related types.
	/// @return A list of basename:base residue type pairs
	BaseTypeList
	get_base_residue_types( core::chemical::ResidueTypeSetCOP const & restypeset ) const override;

private: //Private setup functions:

	/// @brief Set up the NCAADefaultPackerPalette with the default set of position-specific behaviours.
	///
	void set_up_behaviours();

private: //Private member variables:



};


} //namespace palette
} //namespace pack
} //namespace core

#endif
