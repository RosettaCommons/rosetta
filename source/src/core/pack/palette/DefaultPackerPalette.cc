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
#include <basic/citation_manager/CitationManager.hh>
#include <basic/citation_manager/CitationCollection.hh>

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

/// @brief Default constructor
DefaultPackerPalette::DefaultPackerPalette() :
	PackerPalette()
	//TODO -- initialize all private member vars here.
{
	set_up_behaviours(); //Set up the default behaviours.
}

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

BaseTypeList
DefaultPackerPalette::get_base_residue_types( core::chemical::ResidueTypeSetCOP const & restypeset ) const {
	BaseTypeList base_types;

	if ( restypeset ) {
		parent::set_up_default_base_types( *restypeset, base_types );
	}

	return base_types;
}

/// @brief Get the name of this object ("DefaultPackerPalette").
std::string const &
DefaultPackerPalette::name() const {
	static const std::string myname( "DefaultPackerPalette" );
	return myname;
}

/// @brief This packer palette does provide citation info.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
bool
DefaultPackerPalette::packer_palette_provides_citation_info() const {
	return true;
}

/// @brief Provide the citation (Mulligan et al. 2020).
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
utility::vector1< basic::citation_manager::CitationCollectionCOP >
DefaultPackerPalette::provide_citation_info() const {
	basic::citation_manager::CitationCollectionOP my_citation(
		utility::pointer::make_shared< basic::citation_manager::CitationCollection >(
		name(), basic::citation_manager::CitedModuleType::PackerPalette
		)
	);
	my_citation->add_citation( basic::citation_manager::CitationManager::get_instance()->get_citation_by_doi( "Mulligan_2020_underreview" ) );
	return utility::vector1< basic::citation_manager::CitationCollectionCOP > { my_citation };
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
