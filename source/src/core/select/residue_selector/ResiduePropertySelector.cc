// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/residue_selector/ResiduePropertySelector.hh
/// @brief  A residue selector that selects based on set residue properties.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <core/select/residue_selector/ResiduePropertySelector.hh>
#include <core/select/residue_selector/ResiduePropertySelectorCreator.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh> // for xml schema utility functions

// Project headers
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/util.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <utility/assert.hh>
#include <map>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


static basic::Tracer TR( "core.select.residue_selector.ResiduePropertySelector" );

namespace core {
namespace select {
namespace residue_selector {

/// @brief Constructor.
///
ResiduePropertySelector::ResiduePropertySelector():
	core::select::residue_selector::ResidueSelector()
{
}

ResiduePropertySelector::ResiduePropertySelector( chemical::ResidueProperty property):
	core::select::residue_selector::ResidueSelector()
{
	set_property(property);
}

ResiduePropertySelector::ResiduePropertySelector(utility::vector1< chemical::ResidueProperty > properties ):
	core::select::residue_selector::ResidueSelector()
{
	set_properties(properties);
}

/// @brief Destructor.
///
ResiduePropertySelector::~ResiduePropertySelector() {}

/// @brief Copy Constructor.  Usually not necessary unless you need deep copying (e.g. OPs)
//ResiduePropertySelector::ResiduePropertySelector(ResiduePropertySelector const & src):
// core::select::residue_selector::ResidueSelector( src )
//{
//}

/// @brief Clone function.
/// @details Copy this object and return owning pointer to the copy (created on the heap).
core::select::residue_selector::ResidueSelectorOP
ResiduePropertySelector::clone() const {
	return utility::pointer::make_shared< ResiduePropertySelector >(*this);
}

void
ResiduePropertySelector::set_property(core::chemical::ResidueProperty property){
	properties_.clear();
	properties_.push_back( property );
}

void
ResiduePropertySelector::add_property(core::chemical::ResidueProperty property){
	properties_.push_back( property );
}

void
ResiduePropertySelector::set_properties(utility::vector1<core::chemical::ResidueProperty> properties){
	properties_ = properties;
}

void
ResiduePropertySelector::set_selection_logic(core::select::residue_selector::basic_selection_logic logic){
	logic_ = logic;
}

/// @brief "Apply" function.
/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
/// indicating whether each residue is selected ("true") or not ("false").
ResiduePropertySelector::ResidueSubset
ResiduePropertySelector::apply(
	core::pose::Pose const & pose
) const {

	utility::vector1< bool > selection(pose.size(), false);
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		bool all_properties = true;
		for ( auto const & property : properties_ ) {
			if ( pose.residue_type(i).properties().has_property( property) ) {
				if ( logic_ == or_logic ) {
					selection[i] = true;
				}
			} else if ( logic_ == and_logic ) {
				all_properties = false;
				continue;
			}
		}
		if ( logic_ == and_logic && all_properties == true ) {
			selection[i] = true;
		}
	}
	return selection;
}



std::string ResiduePropertySelector::get_name() const
{
	return ResiduePropertySelector::class_name();
}

std::string ResiduePropertySelector::class_name()
{
	return "ResiduePropertySelector";
}

/// @brief XML parse.
/// @details Parse RosettaScripts tags and set up this mover.
void
ResiduePropertySelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & )
{

	properties_.clear();
	utility::vector1< std::string > properties = utility::string_split(tag->getOption< std::string>("properties"), ',');

	//Set and check properties
	for ( std::string const & prop : properties ) {
		core::chemical::ResidueProperty const prop_enum( core::chemical::ResidueProperties::get_property_from_string( utility::upper(prop ) ));
		if ( prop_enum == core::chemical::NO_PROPERTY ) {
			utility_exit_with_message("ResiduePropertySelector: unknown property: \"" + prop + "\".  Please see documentation for list of available residue type properties.");
		} else {
			add_property(prop_enum);
		}
	}

	set_selection_logic(logic_map.at(utility::lower(tag->getOption< std::string >("logic", "and_logic"))));
}

void ResiduePropertySelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	//Syntax Example:
	using namespace utility::tag;
	AttributeList attributes;

	std::string properties_description = "A comma-separated list of properties for selection. Current choices are: \n";
	for ( core::Size i(1); i<=static_cast<core::Size>( chemical::N_PROPERTIES ); ++i ) {
		//Loop through all properties
		// get_string_from_property() accesses the global string/properties map that was generated once:
		properties_description += chemical::ResidueProperties::get_string_from_property( static_cast< chemical::ResidueProperty >(i) );

		if ( i < static_cast<core::Size>( chemical::N_PROPERTIES ) ) properties_description += ", ";
	}

	attributes + XMLSchemaAttribute(  "properties", xsct_string_cslist, properties_description );

	utility::vector1< std::string > logic_options;
	for ( auto const & opt : logic_map ) {
		logic_options.push_back(opt.first);
	}

	std::string logic_description = "The logic to use for mutliple properties.  Default is OR logic. Current choices are: \n " + utility::to_string( logic_options);

	utility::tag::add_schema_restrictions_for_strings( xsd, "logic_types", logic_options);

	attributes + XMLSchemaAttribute::attribute_w_default(  "logic", "logic_types", logic_description, "and_logic" );

	core::select::residue_selector::xsd_type_definition_w_attributes( xsd, class_name(), "Author: Jared Adolf-Bryfogle (jadolfbr@gmail.com)\n"
		"A residue selector that selects based on set residue properties.  Default is to use OR logic for multiple properties.  This can be set via selection_logic option.", attributes );

}

core::select::residue_selector::ResidueSelectorOP
ResiduePropertySelectorCreator::create_residue_selector() const {
	return core::select::residue_selector::ResidueSelectorOP( new ResiduePropertySelector );
}

std::string
ResiduePropertySelectorCreator::keyname() const {
	return ResiduePropertySelector::class_name();
}

/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
void
ResiduePropertySelectorCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	ResiduePropertySelector::provide_xml_schema( xsd );
}

} //core
} //select
} //residue_selector

#ifdef    SERIALIZATION

// See the serialization documentation here.  There is a script you can run.
//  https://wiki.rosettacommons.org/index.php/SerializationFAQ

template< class Archive >
void
core::select::residue_selector::ResiduePropertySelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( properties_ ) );
	arc( CEREAL_NVP( logic_));
}

template< class Archive >
void
core::select::residue_selector::ResiduePropertySelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( properties_ );
	arc( logic_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::ResiduePropertySelector );
CEREAL_REGISTER_TYPE( core::select::residue_selector::ResiduePropertySelector )

CEREAL_REGISTER_DYNAMIC_INIT( core_select_residue_selector_ResiduePropertySelector )
#endif // SERIALIZATION
