// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/ResidueNameSelector.hh
/// @brief  The ResidueNameSelector selects residues using a string containing residue names
/// @author Tom Linsky (tlinsky@uw.edu))

// Unit headers
#include <core/select/residue_selector/ResidueNameSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/pose/selection.hh>
#include <core/conformation/Residue.hh>

// Utility Headers
#include <utility/tag/Tag.hh>

// C++ headers
#include <utility/assert.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {

ResidueNameSelector::ResidueNameSelector():
	res_name_str_(""),
	res_name3_str_("")
{
}

/// @brief Copy constructor
///
ResidueNameSelector::ResidueNameSelector( ResidueNameSelector const &src) :
	res_name_str_( src.res_name_str_ ),
	res_name3_str_( src.res_name3_str_ )
{}

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
ResidueSelectorOP ResidueNameSelector::clone() const { return ResidueSelectorOP( new ResidueNameSelector(*this) ); }

ResidueNameSelector::ResidueNameSelector( std::string const & res_name_str ):
	res_name3_str_()
{
	res_name_str_ = res_name_str;
}


ResidueNameSelector::~ResidueNameSelector() {}

ResidueSubset
ResidueNameSelector::apply( core::pose::Pose const & pose ) const
{
	debug_assert( !res_name3_str_.empty() || !res_name_str_.empty() );

	ResidueSubset subset( pose.total_residue(), false );
	// quit here if there are no residues in the pose
	if ( pose.total_residue() == 0 ) {
		return subset;
	}

	utility::vector1< std::string > const res_name_vec = utility::string_split( res_name_str_, ',' );
	std::set< std::string > res_set;
	for ( utility::vector1< std::string >::const_iterator n=res_name_vec.begin(), endn=res_name_vec.end(); n!=endn; ++n ) {
		if ( n->empty() ) {
			continue;
		}
		// check if the given name is valid
		if ( !pose.residue(1).residue_type_set()->has_name( *n ) ) {
			std::stringstream err;
			err << "ResidueNameSelector: " << *n << " is not a valid residue type name.";
			throw utility::excn::EXCN_BadInput( err.str() );
		}
		res_set.insert( *n );
	}

	utility::vector1< std::string > const res_name3_vec = utility::string_split( res_name3_str_, ',' );
	std::set< std::string > res_name3_set;
	for ( utility::vector1< std::string >::const_iterator n=res_name3_vec.begin(), endn=res_name3_vec.end(); n!=endn; ++n ) {
		if ( n->empty() ) {
			continue;
		}
		// check if the given name is valid
		if ( !pose.residue(1).residue_type_set()->has_name3( *n ) ) {
			std::stringstream err;
			err << "ResidueNameSelector: " << *n << " is not a valid residue type name.";
			throw utility::excn::EXCN_BadInput( err.str() );
		}
		res_name3_set.insert( *n );
	}

	for ( core::Size i=1, endi=pose.total_residue(); i<=endi; ++i ) {
		if ( res_set.find( pose.residue(i).name() ) != res_set.end() ) {
			subset[i] = true;
			continue;
		}
		if ( res_name3_set.find( pose.residue(i).name3() ) != res_name3_set.end() ) {
			subset[i] = true;
		}
	}

	return subset;
}

void
ResidueNameSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &)
{
	if ( tag->hasOption( "residue_name3" ) ) {
		set_residue_name3( tag->getOption< std::string >( "residue_name3" ) );
	}
	if ( tag->hasOption( "residue_names" ) ) {
		set_residue_names( tag->getOption< std::string >( "residue_names" ) );
	}
	if ( res_name_str_.empty() && res_name3_str_.empty() ) {
		std::stringstream err_msg;
		err_msg << "ResidueName selector requires either 'residue_names' or 'residue_name3' to be specified. From ResidueNameSelector::parse_my_tag." << std::endl;
		throw utility::excn::EXCN_BadInput( err_msg.str() );
	}
}

void
ResidueNameSelector::set_residue_names( std::string const &res_name_str )
{
	res_name_str_ = res_name_str;
}

void
ResidueNameSelector::set_residue_name3( std::string const &res_name3_str )
{
	res_name3_str_ = res_name3_str;
}

std::string ResidueNameSelector::get_name() const {
	return ResidueNameSelector::class_name();
}

std::string ResidueNameSelector::class_name() {
	return "ResidueName";
}

ResidueSelectorOP
ResidueNameSelectorCreator::create_residue_selector() const {
	return ResidueSelectorOP( new ResidueNameSelector );
}

std::string
ResidueNameSelectorCreator::keyname() const {
	return ResidueNameSelector::class_name();
}

void ResidueNameSelectorCreator::provide_selector_xsd( utility::tag::XMLSchemaDefinition & ) const
{
	// TO DO!
}

} //namespace residue_selector
} //namespace select
} //namespace core


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::select::residue_selector::ResidueNameSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( res_name_str_ ) ); // std::string
	arc( CEREAL_NVP( res_name3_str_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::select::residue_selector::ResidueNameSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( res_name_str_ ); // std::string
	arc( res_name3_str_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::ResidueNameSelector );
CEREAL_REGISTER_TYPE( core::select::residue_selector::ResidueNameSelector )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_task_residue_selector_ResidueNameSelector )
#endif // SERIALIZATION
