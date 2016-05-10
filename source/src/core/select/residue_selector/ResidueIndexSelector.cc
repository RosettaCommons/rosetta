// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/ResidueIndexSelector.hh
/// @brief  The ResidueIndexSelector selects residues using a string containing pose indices
/// @author Robert Lindner (rlindner@mpimf-heidelberg.mpg.de)

// Unit headers
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>

// Package headers
#include <core/select/residue_selector/util.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Project headers
#include <core/pose/selection.hh>
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

ResidueIndexSelector::ResidueIndexSelector():
	index_str_() {}

/// @brief Copy constructor
///
ResidueIndexSelector::ResidueIndexSelector( ResidueIndexSelector const &src) :
	index_str_( src.index_str_ )
{}

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
ResidueSelectorOP ResidueIndexSelector::clone() const { return ResidueSelectorOP( new ResidueIndexSelector(*this) ); }


ResidueIndexSelector::ResidueIndexSelector( std::string const & index_str )
{
	index_str_ = index_str;
}


ResidueIndexSelector::~ResidueIndexSelector() {}

ResidueSubset
ResidueIndexSelector::apply( core::pose::Pose const & pose ) const
{
	debug_assert( !index_str_.empty() );

	ResidueSubset subset( pose.total_residue(), false );
	std::set< Size > const res_set( get_resnum_list( index_str_, pose ) );

	for ( std::set< Size >::const_iterator it = res_set.begin();
			it != res_set.end(); ++it ) {
		if ( *it == 0 || *it > subset.size() ) {
			std::stringstream err_msg;
			err_msg << "Residue " << *it << " not found in pose!\n";
			throw utility::excn::EXCN_Msg_Exception( err_msg.str() );
		}
		subset[ *it ] = true;
	}
	return subset;
}

void
ResidueIndexSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &)
{
	try {
		set_index( tag->getOption< std::string >( "resnums" ) );
	} catch ( utility::excn::EXCN_Msg_Exception & e ) {
		std::stringstream err_msg;
		err_msg << "Failed to access required option 'resnums' from ResidueIndexSelector::parse_my_tag.\n";
		err_msg << e.msg();
		throw utility::excn::EXCN_Msg_Exception( err_msg.str() );
	}
}

void
ResidueIndexSelector::set_index( std::string const &index_str )
{
	index_str_ = index_str;
}

std::string ResidueIndexSelector::get_name() const {
	return ResidueIndexSelector::class_name();
}

std::string ResidueIndexSelector::class_name() {
	return "Index";
}

void
ResidueIndexSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attributes;
	attributes.push_back( XMLSchemaAttribute::required_attribute( "resnums", xsct_int_cslist ));
	xsd_type_definition_w_attributes( xsd, class_name(), attributes );
}


ResidueSelectorOP
ResidueIndexSelectorCreator::create_residue_selector() const {
	return ResidueSelectorOP( new ResidueIndexSelector );
}

std::string
ResidueIndexSelectorCreator::keyname() const {
	return ResidueIndexSelector::class_name();
}

void
ResidueIndexSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ResidueIndexSelector::provide_xml_schema( xsd );
}

} //namespace residue_selector
} //namespace select
} //namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::select::residue_selector::ResidueIndexSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( index_str_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::select::residue_selector::ResidueIndexSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( index_str_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::ResidueIndexSelector );
CEREAL_REGISTER_TYPE( core::select::residue_selector::ResidueIndexSelector )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_task_residue_selector_ResidueIndexSelector )
#endif // SERIALIZATION
