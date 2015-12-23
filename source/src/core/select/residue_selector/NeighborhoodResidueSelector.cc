// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/NeighborhoodResidueSelector.hh
/// @brief  The NeighborhoodResidueSelector selects residues in a given proximity of set focus residues
/// @author Robert Lindner (rlindner@mpimf-heidelberg.mpg.de)

// Unit headers
#include <core/select/residue_selector/NeighborhoodResidueSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/select/residue_selector/util.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>

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
#include <cereal/types/set.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {

NeighborhoodResidueSelector::NeighborhoodResidueSelector():
	focus_(),
	focus_str_(""),
	distance_(10.0),
	focus_selector_(),
	focus_set_(false),
	use_focus_selector_(false)
{}

/// @brief Copy constructor
///
NeighborhoodResidueSelector::NeighborhoodResidueSelector( NeighborhoodResidueSelector const &src) :
	focus_( src.focus_ ),
	focus_str_( src.focus_str_ ),
	distance_( src.distance_ ),
	focus_selector_( src.focus_selector_ ),
	focus_set_( src.focus_set_ ),
	use_focus_selector_( src.use_focus_selector_ )
{}

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
ResidueSelectorOP NeighborhoodResidueSelector::clone() const { return ResidueSelectorOP( new NeighborhoodResidueSelector(*this) ); }

NeighborhoodResidueSelector::NeighborhoodResidueSelector( std::set<core::Size> const & focus, Real distance )
{
	set_focus(focus);
	set_distance(distance);
}


NeighborhoodResidueSelector::~NeighborhoodResidueSelector() {}

ResidueSubset
NeighborhoodResidueSelector::apply( core::pose::Pose const & pose ) const
{
	debug_assert( focus_set_ );

	ResidueSubset subset( pose.total_residue(), false );
	std::set< Size > focus_tmp;

	// set subset to focus.
	get_focus(pose, subset, focus_tmp);

	Real const dst_squared = distance_ * distance_;
	// go through each residue of the pose and check if it's near anything in the focus set
	for ( Size ii = 1; ii < subset.size() ; ++ii ) {
		if ( subset[ ii ] ) continue;
		conformation::Residue const & r1( pose.residue( ii ) );
		for ( std::set< Size >::const_iterator it = focus_tmp.begin();
				it != focus_tmp.end(); ++it ) {
			conformation::Residue const & r2( pose.residue( *it ) );
			Real const d_sq( r1.xyz( r1.nbr_atom() ).distance_squared( r2.xyz( r2.nbr_atom() ) ) );
			if ( d_sq <= dst_squared ) {
				subset[ ii ] = true;
			}
		} // focus set
	} // subset

	return subset;
}

void
NeighborhoodResidueSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap)
{
	if ( tag->hasOption("selector") ) {
		if ( tag->hasOption("resnums") ) {
			throw utility::excn::EXCN_Msg_Exception( "NeighborhoodResidueSelector takes EITHER 'selector' OR 'resnum' options, not both!\n" );
		}
		if ( tag->size() > 1 ) { // 1 if no subtags exist
			throw utility::excn::EXCN_Msg_Exception( "NeighborhoodResidueSelector can only have one ResidueSelector loaded!\n" );
		}
		// grab the ResidueSelector from the selector option
		// and then grab each of the indicated residue selectors from the datamap.
		std::string selector_str;
		try {
			selector_str = tag->getOption< std::string >( "selector" );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::stringstream error_msg;
			error_msg << "Failed to access option 'selector' from NeighborhoodResidueSelector::parse_my_tag.\n";
			error_msg << e.msg();
			throw utility::excn::EXCN_Msg_Exception( error_msg.str() );
		}

		try {
			ResidueSelectorCOP selector = datamap.get_ptr< ResidueSelector const >( "ResidueSelector", selector_str );
			set_focus_selector( selector );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::stringstream error_msg;
			error_msg << "Failed to find ResidueSelector named '" << selector_str << "' from the Datamap from NeighborhoodResidueSelector::parse_my_tag.\n";
			error_msg << e.msg();
			throw utility::excn::EXCN_Msg_Exception( error_msg.str() );
		}
	} else if ( tag->size() > 1 ) { // get focus selector from tag
		if ( tag->hasOption("resnums") ) {
			throw utility::excn::EXCN_Msg_Exception( "NeighborhoodResidueSelector takes EITHER a 'resnums' tag or a selector subtag, not both!\n" );
		}

		utility::vector0< utility::tag::TagCOP > const & tags = tag->getTags();
		if ( tags.size() > 1 ) {
			throw utility::excn::EXCN_Msg_Exception( "NeighborhoodResidueSelector takes at most one ResidueSelector to determine the focus!\n" );
		}
		ResidueSelectorCOP rs = ResidueSelectorFactory::get_instance()->new_residue_selector(
			tags.front()->getName(),
			tags.front(),
			datamap
		);
		set_focus_selector( rs );

	} else { // do not get focus from ResidueSelectors but load resnums string instead
		try {
			set_focus ( tag->getOption< std::string >( "resnums" ) );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::stringstream err_msg;
			err_msg << "Failed to access option 'resnums' from NeighborhoodResidueSelector::parse_my_tag.\n";
			err_msg << e.msg();
			throw utility::excn::EXCN_Msg_Exception( err_msg.str() );
		}
	}

	// finally grab distance
	try {
		set_distance( tag->getOption< Real >( "distance", 10.0 ) );
	}  catch ( utility::excn::EXCN_Msg_Exception & e ) {
		std::stringstream error_msg;
		error_msg << "Failed to access option 'distance' from NeighborhoodResidueSelector::parse_my_tag.\n";
		error_msg << e.msg();
		throw utility::excn::EXCN_Msg_Exception( error_msg.str() );
	}

}

void
NeighborhoodResidueSelector::get_focus(
	core::pose::Pose const & pose,
	ResidueSubset & subset,
	std::set< Size > & focus
) const
{
	if ( focus_selector_ && use_focus_selector_ ) {
		subset = focus_selector_->apply( pose );

		for ( Size ii = 1; ii <= subset.size(); ++ii ) {
			if ( subset[ ii ] ) {
				focus.insert( ii );
			}
		}
	} else { // grab from set and string
		std::set< Size > const res_vec( get_resnum_list( focus_str_, pose ) );
		focus.insert( focus_.begin(), focus_.end() );
		focus.insert( res_vec.begin(), res_vec.end() );
		for ( std::set< Size >::const_iterator it = focus.begin();
				it != focus.end(); ++it ) {
			if ( *it == 0 || *it > subset.size() ) {
				std::stringstream err_msg;
				err_msg << "Residue " << *it << " not found in pose!\n";
				throw utility::excn::EXCN_Msg_Exception( err_msg.str() );
			}
			subset[ *it ] = true; // may want to use a tmp subset so we don't wind up with a half-set subset
		}
	}
}

void
NeighborhoodResidueSelector::set_focus( std::set<Size> const &focus )
{
	focus_ = focus;
	focus_set_ = true;
	use_focus_selector_ = false;
}

void
NeighborhoodResidueSelector::set_focus( std::string const &focus_str )
{
	focus_str_ = focus_str;
	focus_set_ = true;
	use_focus_selector_ = false;
}

void NeighborhoodResidueSelector::set_focus_selector( ResidueSelectorCOP rs )
{
	focus_selector_ = rs;
	focus_set_ = true;
	use_focus_selector_ = true;
}

void
NeighborhoodResidueSelector::set_distance( Real distance )
{
	distance_ = distance;
}


std::string NeighborhoodResidueSelector::get_name() const {
	return NeighborhoodResidueSelector::class_name();
}

std::string NeighborhoodResidueSelector::class_name() {
	return "Neighborhood";
}

void
NeighborhoodResidueSelector::provide_selector_xsd( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attributes;
	attributes.push_back( XMLSchemaAttribute( "selector", xs_string        ));
	attributes.push_back( XMLSchemaAttribute( "resnums",  "int_cslist"     ));
	attributes.push_back( XMLSchemaAttribute( "distance", xs_decimal, true ));
	xsd_type_definition_w_attributes_and_subselector( xsd, class_name(), attributes );
}


ResidueSelectorOP
NeighborhoodResidueSelectorCreator::create_residue_selector() const {
	return ResidueSelectorOP( new NeighborhoodResidueSelector );
}

std::string
NeighborhoodResidueSelectorCreator::keyname() const {
	return NeighborhoodResidueSelector::class_name();
}

void
NeighborhoodResidueSelectorCreator::provide_selector_xsd( utility::tag::XMLSchemaDefinition & xsd ) const {
	NeighborhoodResidueSelector::provide_selector_xsd( xsd );
}

} //namespace residue_selector
} //namespace select
} //namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::select::residue_selector::NeighborhoodResidueSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( focus_ ) ); // std::set<Size>
	arc( CEREAL_NVP( focus_str_ ) ); // std::string
	arc( CEREAL_NVP( distance_ ) ); // Real
	arc( CEREAL_NVP( focus_selector_ ) ); // ResidueSelectorCOP
	arc( CEREAL_NVP( focus_set_ ) ); // _Bool
	arc( CEREAL_NVP( use_focus_selector_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::select::residue_selector::NeighborhoodResidueSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( focus_ ); // std::set<Size>
	arc( focus_str_ ); // std::string
	arc( distance_ ); // Real
	std::shared_ptr< core::select::residue_selector::ResidueSelector > local_focus_selector;
	arc( local_focus_selector ); // ResidueSelectorCOP
	focus_selector_ = local_focus_selector; // copy the non-const pointer(s) into the const pointer(s)
	arc( focus_set_ ); // _Bool
	arc( use_focus_selector_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::NeighborhoodResidueSelector );
CEREAL_REGISTER_TYPE( core::select::residue_selector::NeighborhoodResidueSelector )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_task_residue_selector_NeighborhoodResidueSelector )
#endif // SERIALIZATION
