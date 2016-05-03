// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/NotResidueSelector.cc
/// @brief  The NotResidueSelector negates the logic of its loaded ResidueSelector
/// @author Robert Lindner (rlindner@mpimf-heidelberg.mpg.de)
/// @see    core::select::residue_selector::AndResidueSelector

// Unit headers
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>

// Package headers
#include <core/select/residue_selector/util.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>

// C++ headers
#include <utility/assert.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {


NotResidueSelector::NotResidueSelector() {}
NotResidueSelector::~NotResidueSelector() {}

/// @brief Copy constructor
///
NotResidueSelector::NotResidueSelector( NotResidueSelector const &src) :
	selector_( src.selector_ )
{}

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
ResidueSelectorOP NotResidueSelector::clone() const { return ResidueSelectorOP( new NotResidueSelector(*this) ); }

NotResidueSelector::NotResidueSelector( ResidueSelectorCOP selector )
{
	set_residue_selector( selector );
}

ResidueSubset
NotResidueSelector::apply( core::pose::Pose const & pose ) const
{
	debug_assert( selector_ );

	ResidueSubset subset = selector_->apply( pose );
	subset.flip();
	return subset;
}

void NotResidueSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
)
{
	if ( tag->hasOption( "selector" ) ) { // fetch selector from datamap

		if ( tag->size() > 1 ) { // has subtags
			throw utility::excn::EXCN_Msg_Exception( "NotResidueSelector can negate ONE ResidueSelector! Either specify 'selector' option or provide subtags but not BOTH\n" );
		}
		// grab the ResidueSelector to be negated from the selector option
		// and then grab each of the indicated residue selectors from the datamap.
		std::string selector_str;
		try {
			selector_str = tag->getOption< std::string >( "selector" );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::stringstream error_msg;
			error_msg << "Failed to access required option 'selector' from NotResidueSelector::parse_my_tag.\n";
			error_msg << e.msg();
			throw utility::excn::EXCN_Msg_Exception( error_msg.str() );
		}

		try {
			ResidueSelectorCOP selector = datamap.get_ptr< ResidueSelector const >( "ResidueSelector", selector_str );
			set_residue_selector(selector);
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::stringstream error_msg;
			error_msg << "Failed to find ResidueSelector named '" << selector_str << "' from the Datamap from NotResidueSelector::parse_my_tag.\n";
			error_msg << e.msg();
			throw utility::excn::EXCN_Msg_Exception( error_msg.str() );
		}
	} else if ( tag->size() > 1 ) { // attempt reading subtag
		utility::vector0< utility::tag::TagCOP > const & tags = tag->getTags();
		if ( tags.size() > 1 ) {
			throw utility::excn::EXCN_Msg_Exception( "NotResidueSelector takes exactly ONE ResidueSelector! Multiple selectors were specified.\n" );
		}
		ResidueSelectorCOP rs = ResidueSelectorFactory::get_instance()->new_residue_selector(
			tags.front()->getName(),
			tags.front(),
			datamap
		);
		set_residue_selector( rs );
	}

	if ( !selector_ ) {
		std::stringstream error_msg;
		error_msg << "No ResidueSelector given to the NotResidueSelector; NotResidueSelector requires a ResidueSelector as input\n";
		throw utility::excn::EXCN_Msg_Exception( error_msg.str() );
	}
}

void NotResidueSelector::set_residue_selector( ResidueSelectorCOP selector )
{
	selector_ = selector;
}

std::string NotResidueSelector::get_name() const {
	return NotResidueSelector::class_name();
}

std::string NotResidueSelector::class_name() {
	return "Not";
}

void
NotResidueSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attributes;
	attributes.push_back( XMLSchemaAttribute( "selector", xs_string ));
	xsd_type_definition_w_attributes_and_optional_subselector( xsd, class_name(), attributes );
}


ResidueSelectorOP
NotResidueSelectorCreator::create_residue_selector() const {
	return ResidueSelectorOP( new NotResidueSelector );
}

std::string
NotResidueSelectorCreator::keyname() const {
	return NotResidueSelector::class_name();
}

void
NotResidueSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	NotResidueSelector::provide_xml_schema( xsd );
}


} //namespace residue_selector
} //namespace select
} //namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::select::residue_selector::NotResidueSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( selector_ ) ); // ResidueSelectorCOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::select::residue_selector::NotResidueSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	std::shared_ptr< core::select::residue_selector::ResidueSelector > local_selector;
	arc( local_selector ); // ResidueSelectorCOP
	selector_ = local_selector; // copy the non-const pointer(s) into the const pointer(s)
}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::NotResidueSelector );
CEREAL_REGISTER_TYPE( core::select::residue_selector::NotResidueSelector )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_task_residue_selector_NotResidueSelector )
#endif // SERIALIZATION
