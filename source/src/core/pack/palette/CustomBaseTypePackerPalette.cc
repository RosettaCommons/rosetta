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
#include <basic/citation_manager/CitationManager.hh>
#include <basic/citation_manager/CitationCollection.hh>

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
	PackerPalette(),
	additional_residue_types_()
	//TODO -- initialize all private member vars here.
{
	set_up_behaviours(); //Set up the default PackerPalette behaviours.
}

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

/// @brief Generate a list of possible base residue types
/// @param [in] restypeset The ResidueTypeSet to use as a reference for related types.
/// @return A list of basename:base residue type pairs
BaseTypeList
CustomBaseTypePackerPalette::get_base_residue_types( core::chemical::ResidueTypeSetCOP const & restypeset ) const {
	BaseTypeList base_types;

	if ( restypeset ) {
		parent::set_up_default_base_types( *restypeset, base_types );

		for ( core::Size i=1, imax=additional_residue_types_.size(); i<=imax; ++i ) {
			parent::add_base_residue_type( additional_residue_types_[i], *restypeset, base_types );
		}
	}

	return base_types;
}

/// @brief Test if this CustomBaseTypePackerPalette has the provided type already.
/// Note that this only tests for explicitly added types. It will not test for default types.
bool
CustomBaseTypePackerPalette::has_type (
	std::string const &type
) const {
	return additional_residue_types_.has_value( type );
}

/// @brief Add a ResidueType (by base type full name -- not 3-letter code) to the set of ResidueTypes
/// being used for design.
void
CustomBaseTypePackerPalette::add_type (
	std::string const &type
) {
	additional_residue_types_.push_back( type );
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

/// @brief This packer palette does provide citation info.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
bool
CustomBaseTypePackerPalette::packer_palette_provides_citation_info() const {
	return true;
}

/// @brief Provide the citation (Mulligan et al. 2020).
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
utility::vector1< basic::citation_manager::CitationCollectionCOP >
CustomBaseTypePackerPalette::provide_citation_info() const {
	basic::citation_manager::CitationCollectionOP my_citation(
		utility::pointer::make_shared< basic::citation_manager::CitationCollection > (
		name(), basic::citation_manager::CitedModuleType::PackerPalette
		)
	);
	my_citation->add_citation( basic::citation_manager::CitationManager::get_instance()->get_citation_by_doi( "Mulligan_2020_underreview" ) );
	return utility::vector1< basic::citation_manager::CitationCollectionCOP > { my_citation };
}

/// @brief Set up the CustomBaseTypePackerPalette with the default set of position-specific behaviours.
void
CustomBaseTypePackerPalette::set_up_behaviours() {
	parent::set_up_default_special_behaviours();
	set_only_design_polymer_residues(false);
	set_only_design_protein_peptoid_dna_saccharide(false);
}

} // palette
} // pack
} // core
