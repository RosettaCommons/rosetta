// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/JumpDownstreamSelector.hh
/// @brief  The JumpDownstreamSelector selects residues downstream of a given jump in a FoldTree
/// @author Robert Lindner (rlindner@mpimf-heidelberg.mpg.de)

// Unit headers
#include <core/select/residue_selector/JumpDownstreamSelector.hh>
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

JumpDownstreamSelector::JumpDownstreamSelector():
	jump_(0) {}

/// @brief Copy constructor
///
JumpDownstreamSelector::JumpDownstreamSelector( JumpDownstreamSelector const &src) :
	jump_( src.jump_ )
{}

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
ResidueSelectorOP JumpDownstreamSelector::clone() const { return ResidueSelectorOP( new JumpDownstreamSelector(*this) ); }

JumpDownstreamSelector::JumpDownstreamSelector( int jump )
{
	jump_ = jump;
}


JumpDownstreamSelector::~JumpDownstreamSelector() {}

ResidueSubset
JumpDownstreamSelector::apply( core::pose::Pose const & pose ) const
{
	debug_assert( jump_ > 0 );

	ResidueSubset subset( pose.total_residue(), false );

	ObjexxFCL::FArray1D_bool upstream( pose.total_residue() );
	pose.fold_tree().partition_by_jump( jump_, upstream );

	for ( core::Size ii = 1; ii < upstream.size(); ++ii ) {
		subset[ ii ] = !upstream( ii );
	}
	return subset;
}

void
JumpDownstreamSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &)
{
	try {
		set_jump( tag->getOption< int >( "jump" ) );
	} catch ( utility::excn::EXCN_Msg_Exception & e ) {
		std::stringstream err_msg;
		err_msg << "Failed to access required option 'jump' from JumpDownstreamSelector::parse_my_tag.\n";
		err_msg << e.msg();
		throw utility::excn::EXCN_Msg_Exception( err_msg.str() );
	}
}

void
JumpDownstreamSelector::set_jump( int jump )
{
	jump_ = jump;
}

std::string JumpDownstreamSelector::get_name() const {
	return JumpDownstreamSelector::class_name();
}

std::string JumpDownstreamSelector::class_name() {
	return "JumpDownstream";
}

void
JumpDownstreamSelector::provide_selector_xsd( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attributes;
	attributes.push_back( XMLSchemaAttribute( "jump", xs_integer, true ));
	xsd_type_definition_w_attributes( xsd, class_name(), attributes );
}

ResidueSelectorOP
JumpDownstreamSelectorCreator::create_residue_selector() const {
	return ResidueSelectorOP( new JumpDownstreamSelector );
}

std::string
JumpDownstreamSelectorCreator::keyname() const {
	return JumpDownstreamSelector::class_name();
}

void
JumpDownstreamSelectorCreator::provide_selector_xsd( utility::tag::XMLSchemaDefinition & xsd ) const {
	return JumpDownstreamSelector::provide_selector_xsd( xsd );
}

} //namespace residue_selector
} //namespace select
} //namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::select::residue_selector::JumpDownstreamSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( jump_ ) ); // int
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::select::residue_selector::JumpDownstreamSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( jump_ ); // int
}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::JumpDownstreamSelector );
CEREAL_REGISTER_TYPE( core::select::residue_selector::JumpDownstreamSelector )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_task_residue_selector_JumpDownstreamSelector )
#endif // SERIALIZATION
