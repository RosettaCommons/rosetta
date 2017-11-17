// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/ResidueNameSelector.hh
/// @brief  The ResidueNameSelector selects residues using a string containing residue names
/// @author Tom Linsky (tlinsky@uw.edu))

// Unit headers
#include <core/select/residue_selector/ResidueNameSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>

// package headers
#include <core/select/residue_selector/util.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/pose/selection.hh>
#include <core/pose/util.hh>
#include <core/conformation/Residue.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

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

	ResidueSubset subset( pose.size(), false );
	// quit here if there are no residues in the pose
	if ( pose.size() == 0 ) {
		return subset;
	}

	utility::vector1< std::string > const res_name_vec = utility::string_split( res_name_str_, ',' );
	std::set< std::string > res_set;
	for ( std::string const & n : res_name_vec ) {
		if ( n.empty() ) {
			continue;
		}
		// check if the given name is valid
		if ( ! pose.residue_type_set_for_pose()->has_name( n ) ) {
			std::stringstream err;
			err << "ResidueNameSelector: " << n << " is not a valid residue type name.";
			throw CREATE_EXCEPTION(utility::excn::BadInput,  err.str() );
		}
		res_set.insert( n );
	}

	utility::vector1< std::string > const res_name3_vec = utility::string_split( res_name3_str_, ',' );
	std::set< std::string > res_name3_set;
	for ( std::string const & n : res_name3_vec ) {
		if ( n.empty() ) {
			continue;
		}
		// check if the given name is valid
		if ( ! pose.residue_type_set_for_pose()->has_name3( n ) ) {
			std::stringstream err;
			err << "ResidueNameSelector: " << n << " is not a valid residue type name.";
			throw CREATE_EXCEPTION(utility::excn::BadInput,  err.str() );
		}
		res_name3_set.insert( n );
	}

	for ( core::Size i=1, endi=pose.size(); i<=endi; ++i ) {
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
		throw CREATE_EXCEPTION(utility::excn::BadInput,  err_msg.str() );
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

void
ResidueNameSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	std::string const req_warning(" Note that one of residue_name3 and residue_names are REQUIRED.");

	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute(
		"residue_name3", xs_string,
		"A comma-separated list of 3-letter Rosetta residue names. "
		"These will be selected regardless of variant type. "
		"For example, 'SER' will select residues named 'SER', "
		"'SER:NtermProteinFull', and 'SER:Phosphorylated'." + req_warning )
		+ XMLSchemaAttribute(
		"residue_names", xs_string,
		"A comma-separated list of Rosetta residue names (including patches). "
		"For example, 'CYD' will select all disulfides, and 'CYD,SER:NTermProteinFull,ALA' "
		"will select all disulfides, alanines, and N-terminal serines -- "
		"all other residues will not be selected (i.e. be false in the ResidueSubset object)." + req_warning );

	xsd_type_definition_w_attributes(
		xsd, class_name(),
		"The ResidueNameSelector selects residues using a string containing residue names" + req_warning,
		attributes );
}

ResidueSelectorOP
ResidueNameSelectorCreator::create_residue_selector() const {
	return ResidueSelectorOP( new ResidueNameSelector );
}

std::string
ResidueNameSelectorCreator::keyname() const {
	return ResidueNameSelector::class_name();
}

void ResidueNameSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ResidueNameSelector::provide_xml_schema( xsd );
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
