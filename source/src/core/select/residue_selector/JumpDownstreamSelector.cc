// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/JumpDownstreamSelector.hh
/// @brief  The JumpDownstreamSelector selects residues downstream of a given jump in a FoldTree
/// @author Robert Lindner (rlindner@mpimf-heidelberg.mpg.de)

// Unit headers
#include <core/select/residue_selector/JumpDownstreamSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>

// Basic Headers
#include <basic/datacache/DataMap.fwd.hh>

// Package headers
#include <core/select/residue_selector/util.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>
#include <core/select/jump_selector/JumpSelector.hh>
#include <core/select/jump_selector/util.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>

// C++ headers


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {

JumpDownstreamSelector::JumpDownstreamSelector() = default;

/// @brief Copy constructor
///
JumpDownstreamSelector::JumpDownstreamSelector( JumpDownstreamSelector const &) = default;

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
ResidueSelectorOP JumpDownstreamSelector::clone() const { return utility::pointer::make_shared< JumpDownstreamSelector >(*this); }

JumpDownstreamSelector::JumpDownstreamSelector( int jump )
{
	jump_ = jump;
}

JumpDownstreamSelector::JumpDownstreamSelector( core::select::jump_selector::JumpSelectorCOP sele ) :
	jump_selector_( sele )
{}

JumpDownstreamSelector::~JumpDownstreamSelector() = default;

ResidueSubset
JumpDownstreamSelector::apply( core::pose::Pose const & pose ) const
{
	core::Size jump = jump_;
	if ( jump_selector_ != nullptr ) {
		utility::vector1< core::Size > const jumps = jump_selector_->selection_jumps( pose );
		runtime_assert_msg( jumps.size() == 1, "Please make sure the selector passed to JumpDownstreamSelector selects exactly one jump. Currently selects " + std::to_string( jumps.size() ) + " jumps!" );
		jump = jumps.front();
	}

	runtime_assert( jump > 0 );

	ResidueSubset subset( pose.size(), false );

	ObjexxFCL::FArray1D_bool upstream( pose.size() );
	pose.fold_tree().partition_by_jump( jump, upstream );

	for ( core::Size ii = 1; ii <= upstream.size(); ++ii ) {
		subset[ ii ] = !upstream( ii );
	}

	core::Size const downstream_resid_for_jump = pose.fold_tree().jump_edge(jump).stop();
	if ( ! subset[downstream_resid_for_jump] ) {
		//We need to flip the bools
		//This sometimes happens with non-standard foldtrees
		for ( core::Size ii = 1; ii <= upstream.size(); ++ii ) {
			subset[ ii ] = !subset[ ii ];
		}
	}
	runtime_assert( subset[downstream_resid_for_jump] );

	return subset;
}

void
JumpDownstreamSelector::parse_my_tag (
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
) {
	bool jump_or_jump_selector_was_provided = false;

	if ( tag->hasOption( "jump" ) ) {
		set_jump( tag->getOption< int >( "jump" ) );
		jump_or_jump_selector_was_provided = true;
	}

	if ( tag->hasOption( "jump_selector" ) ) {
		std::string const jump_selector_name =
			tag->getOption< std::string >( "jump_selector", "" );
		if ( ! jump_selector_name.empty() ) {
			set_jump_selector( core::select::jump_selector::get_jump_selector( jump_selector_name, data ) );
			jump_or_jump_selector_was_provided = true;
		} else {
			set_jump_selector( nullptr );
		}
	}

	runtime_assert_string_msg( jump_or_jump_selector_was_provided, "Please provide either a jump or a jump_selector to the JumpDownstreamSelector" );
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
JumpDownstreamSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute( "jump", xs_integer, "The integer given for the \"jump\" argument should refer to a Jump that is present in the Pose." )
		+ XMLSchemaAttribute::attribute_w_default( "jump_selector", xs_string, "Jump selector to be used as an alternative to the 'jump' option. This selector should only select one jump." , "" );

	xsd_type_definition_w_attributes(
		xsd, class_name(),
		"The JumpDownstreamSelector sets the positions corresponding to "
		"all of the residues that are downstream of the indicated jump to "
		"true, and all the other positions to false. This selector is "
		"logically equivalent to a NotSelector applied to the "
		"JumpUpstreamSelector for the same jump.",
		attributes );
}

ResidueSelectorOP
JumpDownstreamSelectorCreator::create_residue_selector() const {
	return utility::pointer::make_shared< JumpDownstreamSelector >();
}

std::string
JumpDownstreamSelectorCreator::keyname() const {
	return JumpDownstreamSelector::class_name();
}

void
JumpDownstreamSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	JumpDownstreamSelector::provide_xml_schema( xsd );
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
	arc( CEREAL_NVP( jump_selector_ ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::select::residue_selector::JumpDownstreamSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( jump_ ); // int
	arc( jump_selector_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::JumpDownstreamSelector );
CEREAL_REGISTER_TYPE( core::select::residue_selector::JumpDownstreamSelector )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_task_residue_selector_JumpDownstreamSelector )
#endif // SERIALIZATION
