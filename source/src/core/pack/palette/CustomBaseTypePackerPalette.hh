// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/palette/CustomBaseTypePackerPalette.hh
/// @brief  CustomBaseTypePackerPalette: a PackerPalette that allows a user to define additional
/// ResidueTypes with which to design (but not additional VariantTypes, at this time).
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu).

#ifndef INCLUDED_core_pack_palette_CustomBaseTypePackerPalette_hh
#define INCLUDED_core_pack_palette_CustomBaseTypePackerPalette_hh

// Unit Headers
#include <core/pack/palette/CustomBaseTypePackerPalette.fwd.hh>

// Package Headers

// Project Headers
#include <core/pack/palette/PackerPalette.hh>
#include <basic/datacache/DataMap.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>

// STL Headers
#include <utility/vector1.hh>

#include <utility/tag/XMLSchemaGeneration.fwd.hh> // AUTO IWYU For AttributeList


namespace core {
namespace pack {
namespace palette {

/// @brief  CustomBaseTypePackerPalette: a PackerPalette that allows a user to define additional
/// ResidueTypes with which to design (but not additional VariantTypes, at this time).
class CustomBaseTypePackerPalette : public PackerPalette
{

	/// @brief The parent class (PackerPalette).
	///
	typedef PackerPalette parent;

public:
	/// @brief Default constructor.
	///
	CustomBaseTypePackerPalette();

	/// @brief Destructor.
	///
	~CustomBaseTypePackerPalette() override;

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

	/// @brief Add the options for this PackerPalette to the AttributeList.
	static void add_xsd_options( utility::tag::AttributeList & attlist );

	/// @brief Provide information about the XML options available for this PackerPalette.
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief Test if this CustomBaseTypePackerPalette has the provided type already.
	/// Note that this only tests for explicitly added types. It will not test for default types.
	bool has_type( std::string const &type ) const;

	/// @brief Add a ResidueType (by base type full name -- not 3-letter code) to the set of ResidueTypes
	/// being used for design.
	void add_type (std::string const &type);

	/// @brief Given a comma-separated list of additional residue types, separate it out and add the
	/// additonal types to the types used for design.
	/// @details calls CustomBaseTypePackerPalette::add_type().
	void parse_additional_residue_types( std::string const &typelist );

	/// @brief Given a whitespace-separated list of residue base names from a file, parse the file contents
	/// and set up this CustomBaseTypePackerPalette.
	/// @details Appends to any other base types already set up.
	void initialize_from_file_contents( std::string const & file_contents, std::string const & filename );

	/// @brief Get the name of this object ("CustomBaseTypePackerPalette").
	std::string const & name() const override;

public: //Functions for the CitationManager:

	/// @brief Provide the citation.
	void provide_citation_info(basic::citation_manager::CitationCollectionList & ) const override;

protected:

	/// @brief Generate a list of possible base residue types
	/// @param [in] restypeset The ResidueTypeSet to use as a reference for related types.
	/// @return A list of basename:base residue type pairs
	BaseTypeList
	get_base_residue_types( core::chemical::ResidueTypeSetCOP const & restypeset ) const override;

private: //Private setup functions:

	/// @brief Set up the CustomBaseTypePackerPalette with the default set of position-specific behaviours.
	///
	void set_up_behaviours();

private: //Private member variables:

	utility::vector1 < std::string > additional_residue_types_;

};


} //namespace palette
} //namespace pack
} //namespace core

#endif
