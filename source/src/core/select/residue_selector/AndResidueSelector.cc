// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/AndResidueSelector.cc
/// @brief  The AndResidueSelector combines logic from multiple ResidueSelectors
/// @author Robert Lindner (rlindner@mpimf-heidelberg.mpg.de)

// Unit headers
#include <core/select/residue_selector/AndResidueSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>

// Package headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh>

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
#include <cereal/types/list.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {


AndResidueSelector::AndResidueSelector() {}

/// @brief Copy constructor
///
AndResidueSelector::AndResidueSelector( AndResidueSelector const &src) :
	selectors_( src.selectors_ )
{}

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
ResidueSelectorOP AndResidueSelector::clone() const { return ResidueSelectorOP( new AndResidueSelector(*this) ); }

AndResidueSelector::~AndResidueSelector() {}

AndResidueSelector::AndResidueSelector(ResidueSelectorCOP selector1)
{
	add_residue_selector( selector1 );
}

AndResidueSelector::AndResidueSelector( ResidueSelectorCOP selector1, ResidueSelectorCOP selector2 )
{
	add_residue_selector( selector1 );
	add_residue_selector( selector2 );
}

ResidueSubset
AndResidueSelector::apply( core::pose::Pose const & pose ) const
{
	debug_assert( num_selectors() > 0 );

	// make subset neutral for AND operations
	ResidueSubset subset( pose.total_residue(), true );
	for ( std::list< ResidueSelectorCOP >::const_iterator
			rs = selectors_.begin();
			rs != selectors_.end();
			++rs ) {
		ResidueSubset tmp = (*rs)->apply( pose );
		apply_and_to_subset(tmp, subset);
	}
	return subset;
}

void AndResidueSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
)
{
	// grab the comma-separated list of residue selectors that should be ANDed together
	// from the tag, and then grab each of the indicated residue selectors from the datamap.

	std::list< ResidueSelectorCOP > local_selectors;
	if ( tag->hasOption("selectors") ) {
		std::string selectors_str;
		try {
			selectors_str = tag->getOption< std::string >( "selectors" );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::stringstream error_msg;
			error_msg << "Failed to access option 'selectors' from AndResidueSelector::parse_my_tag.\n";
			throw utility::excn::EXCN_Msg_Exception( error_msg.str() );
		}
		utility::vector1< std::string > selector_names = utility::string_split( selectors_str, ',' );

		for ( core::Size ii = 1; ii <= selector_names.size(); ++ii ) {
			try {
				ResidueSelectorCOP selector = datamap.get_ptr< ResidueSelector const >( "ResidueSelector", selector_names[ ii ] );
				local_selectors.push_back( selector );
			} catch ( utility::excn::EXCN_Msg_Exception & e ) {
				std::stringstream error_msg;
				error_msg << "Failed to find ResidueSelector named '" << selector_names[ ii ] << "' from the Datamap from AndResidueSelector::parse_my_tag.\n";
				error_msg << e.msg();
				throw utility::excn::EXCN_Msg_Exception( error_msg.str() );
			}
		}
	} // hasOption selectors

	// add selectors from tags
	for ( utility::vector0< utility::tag::TagCOP >::const_iterator itag = tag->getTags().begin();
			itag != tag->getTags().end(); ++itag ) {
		ResidueSelectorCOP rs = ResidueSelectorFactory::get_instance()->new_residue_selector(
			(*itag)->getName(),
			(*itag),
			datamap
		);
		local_selectors.push_back( rs );
	}

	if ( local_selectors.empty() ) { //size() == 0 ) {
		std::stringstream error_msg;
		error_msg << "No ResidueSelectors given to the AndResidueSelector; AndResidueSelector requires at least one ResidueSelector as input\n";
		throw utility::excn::EXCN_Msg_Exception( error_msg.str() );
	}

	for ( std::list< ResidueSelectorCOP >::const_iterator
			iter = local_selectors.begin(), iter_end = local_selectors.end();
			iter != iter_end; ++iter ) {
		add_residue_selector( *iter );
	}
}

void AndResidueSelector::add_residue_selector( ResidueSelectorCOP selector )
{
	selectors_.push_back(selector);
}

Size AndResidueSelector::num_selectors() const
{
	return selectors_.size();
}

void
AndResidueSelector::apply_and_to_subset(ResidueSubset const & newSubset, ResidueSubset & existingSubset) const
{
	debug_assert( existingSubset.size() == newSubset.size() );
	for ( Size ii = 1; ii <= existingSubset.size(); ++ii ) {
		existingSubset[ ii ] = existingSubset[ ii ] && newSubset[ ii ];
	}
}

std::string AndResidueSelector::get_name() const {
	return AndResidueSelector::class_name();
}

std::string AndResidueSelector::class_name() {
	return "And";
}

void
AndResidueSelector::provide_selector_xsd( utility::tag::XMLSchemaDefinition & xsd ) {
	AttributeList no_attributes;
	xsd_type_definition_w_attributes_and_subselectors( xsd, class_name(), 2, no_attributes );
}


ResidueSelectorOP
AndResidueSelectorCreator::create_residue_selector() const {
	return ResidueSelectorOP( new AndResidueSelector );
}

std::string
AndResidueSelectorCreator::keyname() const {
	return AndResidueSelector::class_name();
}

void
AndResidueSelectorCreator::provide_selector_xsd( utility::tag::XMLSchemaDefinition & xsd ) const {
	AndResidueSelector::provide_selector_xsd( xsd );
}


} //namespace residue_selector
} //namespace select
} //namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::select::residue_selector::AndResidueSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( selectors_ ) ); // std::list<ResidueSelectorCOP>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::select::residue_selector::AndResidueSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	std::list< std::shared_ptr< core::select::residue_selector::ResidueSelector > > local_selectors;
	arc( local_selectors ); // std::list<ResidueSelectorCOP>
	for ( auto iter = local_selectors.begin(); iter != local_selectors.end(); ++iter ) {
		selectors_.push_back( *iter );
	}

}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::AndResidueSelector );
CEREAL_REGISTER_TYPE( core::select::residue_selector::AndResidueSelector )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_task_residue_selector_AndResidueSelector )
#endif // SERIALIZATION
