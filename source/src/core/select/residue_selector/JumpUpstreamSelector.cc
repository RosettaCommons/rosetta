// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/JumpUpstreamSelector.hh
/// @brief  The JumpUpstreamSelector selects residues downstream of a given jump in a FoldTree
/// @author Robert Lindner (rlindner@mpimf-heidelberg.mpg.de)

// Unit headers
#include <core/select/residue_selector/JumpUpstreamSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/select/residue_selector/util.hh>

// Project headers
#include <core/pose/selection.hh>
#include <core/kinematics/FoldTree.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>

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

JumpUpstreamSelector::JumpUpstreamSelector():
	jump_(0) {}

/// @brief Copy constructor
///
JumpUpstreamSelector::JumpUpstreamSelector( JumpUpstreamSelector const &src) :
	jump_( src.jump_ )
{}

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
ResidueSelectorOP JumpUpstreamSelector::clone() const { return ResidueSelectorOP( new JumpUpstreamSelector(*this) ); }

JumpUpstreamSelector::JumpUpstreamSelector( int jump )
{
	jump_ = jump;
}


JumpUpstreamSelector::~JumpUpstreamSelector() {}

ResidueSubset
JumpUpstreamSelector::apply( core::pose::Pose const & pose ) const
{
	debug_assert( jump_ > 0 );
	ResidueSubset subset( pose.size(), false );

	ObjexxFCL::FArray1D_bool upstream( pose.size() );
	pose.fold_tree().partition_by_jump( jump_, upstream );

	for ( core::Size ii = 1; ii < upstream.size(); ++ii ) {
		subset[ ii ] = upstream( ii );
	}
	return subset;
}

void
JumpUpstreamSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &
)
{
	try {
		set_jump( tag->getOption< int >( "jump" ) );
	} catch ( utility::excn::EXCN_Msg_Exception & e ) {
		std::stringstream err_msg;
		err_msg << "Failed to access required option 'jump' from JumpUpstreamSelector::parse_my_tag.\n";
		err_msg << e.msg();
		throw utility::excn::EXCN_Msg_Exception( err_msg.str() );
	}
}

void
JumpUpstreamSelector::set_jump( int jump )
{
	jump_ = jump;
}

std::string JumpUpstreamSelector::get_name() const {
	return JumpUpstreamSelector::class_name();
}

std::string JumpUpstreamSelector::class_name() {
	return "JumpUpstream";
}

void
JumpUpstreamSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attributes;
	attributes.push_back( XMLSchemaAttribute::required_attribute( "jump", xs_integer ));
	xsd_type_definition_w_attributes( xsd, class_name(), attributes );
}


ResidueSelectorOP
JumpUpstreamSelectorCreator::create_residue_selector() const {
	return ResidueSelectorOP( new JumpUpstreamSelector );
}

std::string
JumpUpstreamSelectorCreator::keyname() const {
	return JumpUpstreamSelector::class_name();
}

void
JumpUpstreamSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	JumpUpstreamSelector::provide_xml_schema( xsd );
}

} //namespace residue_selector
} //namespace select
} //namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::select::residue_selector::JumpUpstreamSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( jump_ ) ); // int
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::select::residue_selector::JumpUpstreamSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( jump_ ); // int
}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::JumpUpstreamSelector );
CEREAL_REGISTER_TYPE( core::select::residue_selector::JumpUpstreamSelector )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_task_residue_selector_JumpUpstreamSelector )
#endif // SERIALIZATION
