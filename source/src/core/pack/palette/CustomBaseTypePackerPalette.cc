// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/palette/CustomBaseTypePackerPalette.cc
/// @brief  CustomBaseTypePackerPalette: a PackerPalette that allows a user to define additional
/// ResidueTypes with which to design (but not additional VariantTypes, at this time).
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu).

// Unit Headers
#include <core/pack/palette/CustomBaseTypePackerPalette.hh>
#include <core/pack/palette/CustomBaseTypePackerPaletteCreator.hh>

// Project Headres
#include <core/pack/palette/PackerPalette.hh>
#include <core/pack/palette/xsd_util.hh>

// Basic Headers
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

namespace core {
namespace pack {
namespace palette {

static basic::Tracer TR( "core.pack.palette.CustomBaseTypePackerPalette" );

/// @brief Creator create_packer_palette function implemented.
PackerPaletteOP
CustomBaseTypePackerPaletteCreator::create_packer_palette() const
{
	return utility::pointer::make_shared< CustomBaseTypePackerPalette >();
}

/// @brief Describe the allowed XML options for a particular PackerPalette subclass.
void
CustomBaseTypePackerPaletteCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	CustomBaseTypePackerPalette::provide_xml_schema(xsd);
}

/// @brief Default constructor
CustomBaseTypePackerPalette::CustomBaseTypePackerPalette() :
	PackerPalette(), //Default constructor -- assumes fa_standard ResidueTypeSet
	additional_residue_types_()
	//TODO -- initialize all private member vars here.
{
	set_up_base_types(); //Set up the default residue types as the base types.
	set_up_behaviours(); //Set up the default PackerPalette behaviours, too.
}

/// @brief Copy constructor
CustomBaseTypePackerPalette::CustomBaseTypePackerPalette(
	CustomBaseTypePackerPalette const &//src
) = default;

/// @brief Destructor
CustomBaseTypePackerPalette::~CustomBaseTypePackerPalette() {}

/// @brief Clone operator.
PackerPaletteOP
CustomBaseTypePackerPalette::clone() const
{
	return utility::pointer::make_shared< CustomBaseTypePackerPalette >( *this );
}

/// @brief Function to parse XML tags, implemented by derived classes.
/// @brief Failure to implement this results in a warning message, but does not prevent compilation or
/// program execution.
void
CustomBaseTypePackerPalette::parse_my_tag(
	utility::tag::TagCOP const &tag,
	basic::datacache::DataMap const &/*datamap*/
) {
	TR << "Parsing XML for " << name() << "." << std::endl;

	std::string const additional_types( tag->getOption<std::string>("additional_residue_types", "") );
	parse_additional_residue_types( additional_types );
}

/// @brief Add the options for this PackerPalette to the AttributeList.
void
CustomBaseTypePackerPalette::add_xsd_options(
	utility::tag::AttributeList & attlist
) {
	using namespace utility::tag;
	attlist + XMLSchemaAttribute::attribute_w_default( "additional_residue_types", xsct_string_cslist,
		"A comma-separated list of additional residue types (by full base name) to add to the PackerPalette.", "" );
}

/// @brief Provide information about the XML options available for this PackerPalette.
void
CustomBaseTypePackerPalette::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) {
	using namespace utility::tag;
	AttributeList attlist;
	add_xsd_options(attlist);
	xsd_type_definition_w_attributes( xsd, "CustomBaseTypePackerPalette", "Sets up a packer palette that expands the default (canonical) residue type set with user-defined base types or types selected by ResidueProperties.", attlist );
}

/// @brief Function to allow a different ResidueTypeSet to be set.
/// @details Each PackerPalette derived class must implement this.  After setting the new ResidueTypeSet, things need to happen.
void
CustomBaseTypePackerPalette::set_residue_type_set(
	core::chemical::ResidueTypeSetCOP new_type_set
) {
	parent::set_residue_type_set( new_type_set );
	set_up_base_types();
	for ( core::Size i=1, imax=additional_residue_types_.size(); i<=imax; ++i ) {
		parent::add_base_residue_type( additional_residue_types_[i] );
	}
	set_up_behaviours();
}

/// @brief Add a ResidueType (by base type full name -- not 3-letter code) to the set of ResidueTypes
/// being used for design.
void
CustomBaseTypePackerPalette::add_type (
	std::string const &type
) {
	additional_residue_types_.push_back( type );
	parent::add_base_residue_type( type );
}

/// @brief Given a comma-separated list of additional residue types, separate it out and add the
/// additonal types to the types used for design.
/// @details calls CustomBaseTypePackerPalette::add_type().
void
CustomBaseTypePackerPalette::parse_additional_residue_types(
	std::string const &typelist
) {
	if ( typelist.empty() ) return; //Do nothing if given an empty string.
	std::string const typelist_copy( utility::strip(typelist, " \t\n" ) ); //Strip whitespace from start and end.
	utility::vector1 < std::string > types_to_add;
	types_to_add = utility::string_split( typelist_copy, ',' );
	for ( core::Size i=1, imax=types_to_add.size(); i<=imax; ++i ) {
		add_type(types_to_add[i]);
	}
}

/// @brief Given a whitespace-separated list of residue base names from a file, parse the file contents
/// and set up this CustomBaseTypePackerPalette.
/// @details Appends to any other base types already set up.
void
CustomBaseTypePackerPalette::initialize_from_file_contents(
	std::string const & file_contents,
	std::string const & filename
) {
	std::stringstream file_contents_ss( file_contents );
	std::string line;
	utility::vector1< std::string > new_basetypes;
	while ( std::getline( file_contents_ss, line ) ) {
		line = line.substr( 0, line.find("#") ); //Strip anything following a comment.
		utility::trim( line, " \t\n"); //Strip whitespace from the line.
		if ( line.empty() ) continue;
		std::stringstream linestream( line );
		while ( !linestream.eof() && !linestream.bad() ) {
			std::string name;
			linestream >> name;
			if ( !new_basetypes.has_value(name) /*Avoiding duplicates*/ ) {
				new_basetypes.push_back(name);
				add_type( name );
			} else {
				TR.Warning << "Warning!  The base type \"" + name + "\" appears more than once in file \"" + filename + "\".  Ignoring instances after the first." << std::endl;
			}
		}
	}
	TR.Warning.flush();
}

/// @brief Get the name of this object ("CustomBaseTypePackerPalette").
std::string const &
CustomBaseTypePackerPalette::name() const {
	static const std::string myname( "CustomBaseTypePackerPalette" );
	return myname;
}

/// @brief Set up the CustomBaseTypePackerPalette with the standard residues.
void
CustomBaseTypePackerPalette::set_up_base_types() {
	if ( TR.Debug.visible() ) TR.Debug << "Setting up default base types (ACDEFGHIKLMNPQRSTVWY/acgt) for " << name() << "." << std::endl;
	parent::set_up_default_base_types();
}

/// @brief Set up the CustomBaseTypePackerPalette with the default set of position-specific behaviours.
void
CustomBaseTypePackerPalette::set_up_behaviours() {
	if ( TR.Debug.visible() ) TR.Debug << "Setting up default PackerPalette behaviours for " << name() << "." << std::endl;
	parent::set_up_default_special_behaviours();
	set_only_design_polymer_residues(false);
	set_only_design_protein_peptoid_dna_saccharide(false);
}

} // palette
} // pack
} // core
