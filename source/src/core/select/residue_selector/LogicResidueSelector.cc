// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/residue_selector/LogicResidueSelector.hh
/// @brief  A residue selector that takes arbitrarily many selectors and performs and/or/not boolean logic within the "selectors" option.
/// @author Frances Chu (francesc345@gmail.com)

// Unit headers
#include <core/select/residue_selector/LogicResidueSelector.hh>
#include <core/select/residue_selector/LogicResidueSelectorCreator.hh>

// Basic Headers
#include <basic/datacache/DataMap.fwd.hh>

// Package headers
#include <core/select/residue_selector/util.hh> // for xml schema utility functions

// Project headers
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/pointer/memory.hh>
#include <basic/Tracer.hh>

// C++ headers

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


static basic::Tracer TR( "core.select.residue_selector.LogicResidueSelector" );

namespace core {
namespace select {
namespace residue_selector {

/// @brief Constructor.
///
LogicResidueSelector::LogicResidueSelector():
	core::select::residue_selector::ResidueSelector()
{
}

/// @brief Destructor.
///
LogicResidueSelector::~LogicResidueSelector() {}

/// @brief Copy Constructor.  Usually not necessary unless you need deep copying (e.g. OPs)
//LogicResidueSelector::LogicResidueSelector(LogicResidueSelector const & src):
// core::select::residue_selector::ResidueSelector( src )
//{
//}

/// @brief Clone function.
/// @details Copy this object and return owning pointer to the copy (created on the heap).
core::select::residue_selector::ResidueSelectorOP
LogicResidueSelector::clone() const {

	return utility::pointer::make_shared< LogicResidueSelector >( *this );
}

/// @brief "Apply" function.
/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
/// indicating whether each residue is selected ("true") or not ("false").
LogicResidueSelector::ResidueSubset
LogicResidueSelector::apply(
	core::pose::Pose const & pose
) const {
	utility::vector1< bool > selection = selector_->apply(pose);
	return selection;
}

/// @brief XML parse.
/// @details Parse RosettaScripts tags and set up this mover.
void
LogicResidueSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap)
{
	std::string selector_name;
	if ( (tag->getOptions().size()==2 && tag->hasOption("name")) || (tag->getOptions().size()==1 && !(tag->hasOption("name"))) ) {
		if ( tag->hasOption("selector") ) {
			try {
				selector_name = tag->getOption < std::string >( "selector");
			} catch ( utility::excn::Exception & e ) {
				std::stringstream error_msg;
				error_msg << "Failed to access option \"selector\" from LogicResidueSelector::parse_my_tag.\n";
				throw CREATE_EXCEPTION(utility::excn::Exception,  error_msg.str() );
			}
		} else if ( tag->hasOption("selectors") ) {
			try {
				selector_name = tag->getOption < std::string >( "selectors");
			} catch ( utility::excn::Exception & e ) {
				std::stringstream error_msg;
				error_msg << "Failed to access option \"selectors\" from LogicResidueSelector::parse_my_tag.\n";
				throw CREATE_EXCEPTION(utility::excn::Exception,  error_msg.str() );
			}
		} else if ( tag->hasOption("residue_selector") ) {
			try {
				selector_name = tag->getOption < std::string >( "residue_selector");
			} catch ( utility::excn::Exception & e ) {
				std::stringstream error_msg;
				error_msg << "Failed to access option \"residue_selector\" from LogicResidueSelector::parse_my_tag.\n";
				throw CREATE_EXCEPTION(utility::excn::Exception,  error_msg.str() );
			}
		} else if ( tag->hasOption("residue_selectors") ) {
			try {
				selector_name = tag->getOption < std::string >( "residue_selectors");
			} catch ( utility::excn::Exception & e ) {
				std::stringstream error_msg;
				error_msg << "Failed to access option \"residue_selectors\" from LogicResidueSelector::parse_my_tag.\n";
				throw CREATE_EXCEPTION(utility::excn::Exception,  error_msg.str() );
			}
		}
	} else {
		std::stringstream error_msg;
		error_msg << "LogicResidueSelector must be given exactly one of the following residue selector attributes: \"selector\", \"selectors\", \"residue_selector\", or \"residue_selectors\". The LogicResidueSelector performs the same function no matter the attribute used.\n";
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  error_msg.str() );
	}
	if ( selector_name.find(',') != std::string::npos ) {
		std::stringstream ss;
		ss << "LogicResidueSelector does not take a comma separated list of residue selectors. Residue selectors should be separated by boolean(s) \"and/or/not\" with spaces on both sides." << std::endl;
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  ss.str() );
	}
	selector_ = core::select::residue_selector::get_residue_selector(selector_name, datamap);
	debug_assert( selector_ );
	if ( selector_ == nullptr ) {
		std::stringstream ss;
		ss << "No ResidueSelectors given to the LogicResidueSelector; LogicResidueSelector requires at least one selector to be set through one of the following attributes: \"selector\", \"selectors\", \"residue_selector\", or \"residue_selectors\". The LogicResidueSelector performs the same function no matter the attribute used.\n" << std::endl;
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  ss.str() );
	}
	TR << "Using residue selectors " << selector_name << std::endl;
}

void LogicResidueSelector::set_residue_selector( ResidueSelectorCOP selector ) {
	selector_ = selector;
}

std::string LogicResidueSelector::get_name() const
{
	return LogicResidueSelector::class_name();
}

std::string LogicResidueSelector::class_name()
{
	return "Logic";
}

void LogicResidueSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	//Syntax Example:
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute("selector", xs_string, "Residue selectors separated by and/or/not boolean logic.")
		+ XMLSchemaAttribute("selectors", xs_string, "Residue selectors separated by and/or/not boolean logic.")
		+ XMLSchemaAttribute("residue_selector", xs_string, "Residue selectors separated by and/or/not boolean logic.")
		+ XMLSchemaAttribute("residue_selectors", xs_string, "Residue selectors separated by and/or/not boolean logic.");
	xsd_type_definition_w_attributes(xsd, class_name(), "Author: Frances Chu (francesc345@gmail.com).\n"
		"A residue selector that takes one of the following options: \"selector\", \"selectors\", \"residue_selector\", \"residue_selectors\", and performs the same function no matter the option given. Takes arbitrarily many residue selectors within the option and performs and/or/not boolean logic." , attributes );
}

core::select::residue_selector::ResidueSelectorOP
LogicResidueSelectorCreator::create_residue_selector() const {
	return utility::pointer::make_shared< LogicResidueSelector >();
}

std::string
LogicResidueSelectorCreator::keyname() const {
	return LogicResidueSelector::class_name();
}

/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
void
LogicResidueSelectorCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	LogicResidueSelector::provide_xml_schema( xsd );
}

} //residue_selector
} //select
} //core

#ifdef    SERIALIZATION

// See the serialization documentation here.  There is a script you can run.
//  https://wiki.rosettacommons.org/index.php/SerializationFAQ

template< class Archive >
void
core::select::residue_selector::LogicResidueSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP ( selector_ ));
}

template< class Archive >
void
core::select::residue_selector::LogicResidueSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	std::shared_ptr< core::select::residue_selector::ResidueSelector > local_selector;
	arc( local_selector );
	selector_ = local_selector;
}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::LogicResidueSelector );
CEREAL_REGISTER_TYPE( core::select::residue_selector::LogicResidueSelector )

CEREAL_REGISTER_DYNAMIC_INIT( core_select_residue_selector_LogicResidueSelector )
#endif // SERIALIZATION
