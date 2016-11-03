// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/environment/claims/EnvLabelSelector.cc
/// @brief  The EnvLabelSelector holds a selection that another object can set inside of it.
/// @author Justin R. Porter

// Unit headers
#include <protocols/environment/claims/EnvLabelSelector.hh>

// Package headers
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>

#include <core/environment/LocalPosition.hh>
#include <protocols/environment/ProtectedConformation.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>

// C++ headers
#include <utility/assert.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace environment {
namespace claims {

EnvLabelSelector::EnvLabelSelector() :
	positions_()
{}

/// @brief Copy constructor
///
EnvLabelSelector::EnvLabelSelector( EnvLabelSelector const &src) :
	positions_( src.positions_ )
{}

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
core::select::residue_selector::ResidueSelectorOP
EnvLabelSelector::clone() const {
	return core::select::residue_selector::ResidueSelectorOP( new EnvLabelSelector(*this) );
}

EnvLabelSelector::EnvLabelSelector( LocalPositions const& positions_in ) {
	this->set_local_positions( positions_in );
}

EnvLabelSelector::EnvLabelSelector( LocalPosition const& local_pos ) {
	LocalPositions local_positions = LocalPositions();
	local_positions.push_back( core::environment::LocalPositionOP( new LocalPosition( local_pos ) ) );

	this->set_local_positions( local_positions );
}

EnvLabelSelector::EnvLabelSelector( std::string const& label,
	std::pair< core::Size, core::Size > const& range ) {
	LocalPositions local_positions = LocalPositions();

	for ( Size i = range.first; i <= range.second; ++i ) {
		local_positions.push_back( core::environment::LocalPositionOP( new LocalPosition( label, i ) ) );
	}

	this->set_local_positions( local_positions );
}


EnvLabelSelector::~EnvLabelSelector() {}


core::select::residue_selector::ResidueSubset
EnvLabelSelector::apply(
	core::pose::Pose const & pose
) const
{
	using core::environment::LocalPositionOP;

	ResidueSubset subset( pose.size(), false );

	ProtectedConformationCOP conf = utility::pointer::dynamic_pointer_cast< protocols::environment::ProtectedConformation const > ( pose.conformation_ptr() );
	core::environment::SequenceAnnotationCOP ann = conf->annotations();

	for ( LocalPositions::const_iterator pos_it = positions_.begin();
			pos_it != positions_.end(); ++pos_it ) {
		core::Size const seqpos = ann->resolve_seq( **pos_it );
		subset[ seqpos ] = true;
	}
	return subset;
}

void EnvLabelSelector::parse_my_tag(
	utility::tag::TagCOP,
	basic::datacache::DataMap & )
{
	throw utility::excn::EXCN_RosettaScriptsOption( "Not to be used in RosettaScripts. For legacy compatibility only." );
}

void EnvLabelSelector::set_local_positions( LocalPositions const& positions_in ){
	using namespace core::environment;

	for ( LocalPositionOP pos : positions_in ) {
		positions_.push_back( core::environment::LocalPositionOP( new LocalPosition( *pos ) ) );
	}
}

std::string EnvLabelSelector::get_name() const {
	return EnvLabelSelector::class_name();
}

std::string EnvLabelSelector::class_name() {
	return "EnvLabel";
}

} //namespace claims
} //namespace environment
} //namespace protocols


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::environment::claims::EnvLabelSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( positions_ ) ); // LocalPositions
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::environment::claims::EnvLabelSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( positions_ ); // LocalPositions
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::environment::claims::EnvLabelSelector );
CEREAL_REGISTER_TYPE( protocols::environment::claims::EnvLabelSelector )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_environment_claims_EnvLabelSelector )
#endif // SERIALIZATION
