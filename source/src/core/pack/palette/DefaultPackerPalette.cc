// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/palette/DefaultPackerPalette.cc
/// @brief  DefaultPackerPalette: a PackerPalette that just sets up
/// default Packer behaviour (design with the canonical 20 amino acids and whatever
/// is present at a position in a pose).  This PackerPalette has no user-configurable
/// options.\nThis was implemented as part of the 2016 Chemical XRW (eXtreme Rosetta
/// Workshop).
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu).

// Unit Headers
#include <core/pack/palette/DefaultPackerPalette.hh>
#include <core/pack/palette/DefaultPackerPaletteCreator.hh>

// Project Headers
#include <core/pack/palette/PackerPalette.hh>
#include <core/pack/palette/xsd_util.hh>

// Basic Headers
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

namespace core {
namespace pack {
namespace palette {

static basic::Tracer TR( "core.pack.palette.DefaultPackerPalette" );

/// @brief Creator create_packer_palette function implemented.
PackerPaletteOP
DefaultPackerPaletteCreator::create_packer_palette() const
{
	return utility::pointer::make_shared< DefaultPackerPalette >();
}

// @brief Describe the allowed XML options for a particular PackerPalette subclass.
void
DefaultPackerPaletteCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	DefaultPackerPalette::provide_xml_schema(xsd);
}

/// @brief Constructor with ResidueTypeSet.
DefaultPackerPalette::DefaultPackerPalette( core::chemical::ResidueTypeSetCOP restypeset) :
	PackerPalette(restypeset) //Default constructor -- assumes fa_standard ResidueTypeSet
	//TODO -- initialize all private member vars here.
{
	set_up_base_types(); //Set up the default residue types as the base types.
	set_up_behaviours(); //Set up the default behaviours.
}


/// @brief Default constructor
DefaultPackerPalette::DefaultPackerPalette() :
	PackerPalette() //Default constructor -- assumes fa_standard ResidueTypeSet
	//TODO -- initialize all private member vars here.
{
	set_up_base_types(); //Set up the default residue types as the base types.
	set_up_behaviours(); //Set up the default behaviours.
}

/// @brief Copy constructor
DefaultPackerPalette::DefaultPackerPalette(
	DefaultPackerPalette const &src
) :
	PackerPalette(src)
	//TODO -- copy all private member vars here.
{}

/// @brief Destructor
DefaultPackerPalette::~DefaultPackerPalette() {}

/// @brief Clone operator.
PackerPaletteOP
DefaultPackerPalette::clone() const
{
	return utility::pointer::make_shared< DefaultPackerPalette >( *this );
}

/// @brief Function to parse XML tags, implemented by derived classes.
/// @brief Failure to implement this results in a warning message, but does not prevent compilation or
/// program execution.
void
DefaultPackerPalette::parse_my_tag(
	utility::tag::TagCOP const &/*tag*/,
	basic::datacache::DataMap const &/*datamap*/
) {
	TR.Debug << "Parsing XML for DefaultPackerPalette.  (No user-configurable options)." << std::endl;
}

/// @brief Provide information about the XML options available for this PackerPalette.
void
DefaultPackerPalette::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) {
	using namespace utility::tag;
	AttributeList attlist;

	xsd_type_definition_w_attributes( xsd, "DefaultPackerPalette", "Sets up a default packer palette, with no user-configurable options.  This permits design with canonical residue types only.  (Note that this is the default behaviour in the absence of a PackerPalette, too.)", attlist );
}

/// @brief Function to allow a different ResidueTypeSet to be set.
/// @details Each PackerPalette derived class must implement this.
/// After setting the new ResidueTypeSet, things need to happen.
void
DefaultPackerPalette::set_residue_type_set(
	core::chemical::ResidueTypeSetCOP new_type_set
) {
	parent::set_residue_type_set( new_type_set );
	set_up_base_types();
	set_up_behaviours();
}

/// @brief Get the name of this object ("DefaultPackerPalette").
std::string const &
DefaultPackerPalette::name() const {
	static const std::string myname( "DefaultPackerPalette" );
	return myname;
}

/// @brief Set up the DefaultPackerPalette with the standard residues.
void
DefaultPackerPalette::set_up_base_types()
{
	if ( TR.Debug.visible() ) TR.Debug << "Setting up default base types (ACDEFGHIKLMNPQRSTVWY/acgt) for DefaultPackerPalette." << std::endl;
	parent::set_up_default_base_types();
}

/// @brief Set up the DefaultPackerPalette with the default set of position-specific behaviours.
void
DefaultPackerPalette::set_up_behaviours()
{
	if ( TR.Debug.visible() ) TR.Debug << "Setting up default PackerPalette behaviours for DefaultPackerPalette." << std::endl;
	parent::set_up_default_special_behaviours();
}

} // palette
} // pack
} // core
