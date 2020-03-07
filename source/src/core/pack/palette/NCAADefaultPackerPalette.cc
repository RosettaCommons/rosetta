// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/palette/NCAADefaultPackerPalette.cc
/// @brief  NCAADefaultPackerPalette: a PackerPalette that just sets up
/// default Packer behaviour (design with the canonical 20 amino acids and equivalent
/// sets for ANY backbone.)
/// @author Andy Watkins (amw579@stanford.edu).

// Unit Headers
#include <core/pack/palette/NCAADefaultPackerPalette.hh>
#include <core/pack/palette/NCAADefaultPackerPaletteCreator.hh>

// Project Headers
#include <core/pack/palette/PackerPalette.hh>
#include <core/pack/palette/xsd_util.hh>

// Basic Headers
#include <basic/Tracer.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

namespace core {
namespace pack {
namespace palette {

static basic::Tracer TR( "core.pack.palette.NCAADefaultPackerPalette" );

/// @brief Creator create_packer_palette function implemented.
PackerPaletteOP
NCAADefaultPackerPaletteCreator::create_packer_palette() const
{
	return utility::pointer::make_shared< NCAADefaultPackerPalette >();
}

// @brief Describe the allowed XML options for a particular PackerPalette subclass.
void
NCAADefaultPackerPaletteCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	NCAADefaultPackerPalette::provide_xml_schema(xsd);
}

/// @brief NCAADefault constructor
NCAADefaultPackerPalette::NCAADefaultPackerPalette() :
	PackerPalette()
	//TODO -- initialize all private member vars here.
{
	set_up_behaviours(); //Set up the default behaviours.
}

/// @brief Destructor
NCAADefaultPackerPalette::~NCAADefaultPackerPalette() {}

/// @brief Clone operator.
PackerPaletteOP
NCAADefaultPackerPalette::clone() const
{
	return utility::pointer::make_shared< NCAADefaultPackerPalette >( *this );
}

/// @brief Function to parse XML tags, implemented by derived classes.
/// @brief Failure to implement this results in a warning message, but does not prevent compilation or
/// program execution.
void
NCAADefaultPackerPalette::parse_my_tag(
	utility::tag::TagCOP const &/*tag*/,
	basic::datacache::DataMap const &/*datamap*/
) {
	TR.Debug << "Parsing XML for NCAADefaultPackerPalette.  (No user-configurable options)." << std::endl;
}

/// @brief Provide information about the XML options available for this PackerPalette.
void
NCAADefaultPackerPalette::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) {
	using namespace utility::tag;
	AttributeList attlist;

	xsd_type_definition_w_attributes( xsd, "NCAADefaultPackerPalette", "Sets up a default packer palette, with no user-configurable options.  This permits design with canonical residue types only.  (Note that this is the default behaviour in the absence of a PackerPalette, too.)", attlist );
}

/// @brief Generate a list of possible base residue types
/// @param [in] restypeset The ResidueTypeSet to use as a reference for related types.
/// @return A list of basename:base residue type pairs
BaseTypeList
NCAADefaultPackerPalette::get_base_residue_types( core::chemical::ResidueTypeSetCOP const & restypeset ) const {
	BaseTypeList base_types;

	if ( restypeset ) {
		parent::set_up_expanded_base_types( *restypeset, base_types );
	}

	return base_types;
}

/// @brief Get the name of this object ("NCAADefaultPackerPalette").
std::string const &
NCAADefaultPackerPalette::name() const {
	static const std::string myname( "NCAADefaultPackerPalette" );
	return myname;
}

/// @brief Returns true (this packer palette is unpublished).
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
bool
NCAADefaultPackerPalette::packer_palette_is_unpublished() const {
	return true;
}

/// @brief Returns Andy Watkins as the author.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
utility::vector1< basic::citation_manager::UnpublishedModuleInfoCOP >
NCAADefaultPackerPalette::provide_authorship_info_for_unpublished() const {
	return utility::vector1< basic::citation_manager::UnpublishedModuleInfoCOP > {
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		name(), basic::citation_manager::CitedModuleType::PackerPalette,
		"Andrew Watkins", "Dept. of Biochemistry, Stanford University", "amw579@stanford.edu"
		)
		};
}

/// @brief Set up the NCAADefaultPackerPalette with the default set of position-specific behaviours.
void
NCAADefaultPackerPalette::set_up_behaviours()
{
	if ( TR.Debug.visible() ) TR.Debug << "Setting up default PackerPalette behaviours for NCAADefaultPackerPalette." << std::endl;
	parent::set_up_default_special_behaviours();
}

} // palette
} // pack
} // core
