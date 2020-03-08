// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/palette/NoDesignPackerPalette.cc
/// @brief  NoDesignPackerPalette: a PackerPalette that sets up absolutely no design residues.
/// @details Not necessary for repacking (the RestrictToRepacking task operation produces the same effect), but handy
/// for efficiency when you know that you're not doing any design.  (There's no point setting up a list of ResidueTypes)
/// only to prune them all away.
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu).

// Unit Headers
#include <core/pack/palette/NoDesignPackerPalette.hh>
#include <core/pack/palette/NoDesignPackerPaletteCreator.hh>

// Project Headers
#include <core/pack/palette/PackerPalette.hh>
#include <core/pack/palette/xsd_util.hh>

// Basic Headers
#include <basic/citation_manager/UnpublishedModuleInfo.hh>
#include <basic/citation_manager/CitationManager.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

namespace core {
namespace pack {
namespace palette {

static basic::Tracer TR( "core.pack.palette.NoDesignPackerPalette" );

/// @brief Creator create_packer_palette function implemented.
PackerPaletteOP
NoDesignPackerPaletteCreator::create_packer_palette() const
{
	return utility::pointer::make_shared< NoDesignPackerPalette >();
}

// @brief Describe the allowed XML options for a particular PackerPalette subclass.
void
NoDesignPackerPaletteCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	NoDesignPackerPalette::provide_xml_schema(xsd);
}

/// @brief Default constructor
NoDesignPackerPalette::NoDesignPackerPalette() :
	PackerPalette()
	//TODO -- initialize all private member vars here.
{
	set_up_behaviours(); //Set up the default behaviours.
}

/// @brief Destructor
NoDesignPackerPalette::~NoDesignPackerPalette() {}

/// @brief Clone operator.
PackerPaletteOP
NoDesignPackerPalette::clone() const
{
	return utility::pointer::make_shared< NoDesignPackerPalette >( *this );
}

/// @brief Function to parse XML tags, implemented by derived classes.
/// @details This parse_my_tag is unusual in that it handles the registration with the CitationManager.  The
/// NoDesignPackerPalette probably doesn't need to be cited unless the user has explicitly scripted it in a
/// RosettaScripts script.  Normally, registration would be handled by provide_authorship_info_for_unpublished().
void
NoDesignPackerPalette::parse_my_tag(
	utility::tag::TagCOP const &/*tag*/,
	basic::datacache::DataMap const &/*datamap*/
) {
	TR.Debug << "Parsing XML for NoDesignPackerPalette.  (No user-configurable options)." << std::endl;

	// Only register this PackerPalette with the CitationManager if it is explicitly RosettaScripts-scripted:
	basic::citation_manager::CitationManager::get_instance()->add_unpublished_modules(
		utility::vector1< basic::citation_manager::UnpublishedModuleInfoCOP > {
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		name(), basic::citation_manager::CitedModuleType::PackerPalette,
		"Vikram K. Mulligan", "Systems Biology, Center for Computational Biology, Flatiron Institute",
		"vmulligan@flatironinstitute.org"
		)
		}
	);
}

/// @brief Provide information about the XML options available for this PackerPalette.
void
NoDesignPackerPalette::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) {
	using namespace utility::tag;
	AttributeList attlist;

	xsd_type_definition_w_attributes( xsd, "NoDesignPackerPalette", "Sets up an empty packer palette, specifying no residues, with no user-configurable options.  This restricts the packer to repacking.  Note that this behaviour can be achieved with the RestrictToRepacking task operation, but doing this with the PackerPalette can be preferable since it means that, under the hood, we're not unnecessarily setting up a list of ResidueTypes only to delete them all.", attlist );
}

/// @brief Get the name of this object ("NoDesignPackerPalette").
std::string const &
NoDesignPackerPalette::name() const {
	static const std::string myname( "NoDesignPackerPalette" );
	return myname;
}

BaseTypeList
NoDesignPackerPalette::get_base_residue_types( core::chemical::ResidueTypeSetCOP const & ) const {
	//Absolutely NO base types for NoDesignPackerPalette.
	return {};
}

/// @brief Set up the NoDesignPackerPalette with the default set of position-specific behaviours.
void
NoDesignPackerPalette::set_up_behaviours()
{
	if ( TR.Debug.visible() ) TR.Debug << "Setting up default PackerPalette behaviours for NoDesignPackerPalette." << std::endl;
	parent::set_up_default_special_behaviours();
	parent::set_force_existing_base_type(true); //Force only the existing base type at each position.
}

} // palette
} // pack
} // core
