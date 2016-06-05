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
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com) - 10 A neighbor graph, simplification, ResidueSubset as focus, clean up, etc.

// Unit headers
#include <core/select/residue_selector/NeighborhoodResidueSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/select/residue_selector/util.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/util.hh>

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
	ResidueSelector(),
	focus_selector_(/* NULL */)
{
	set_defaults();
}





NeighborhoodResidueSelector::NeighborhoodResidueSelector( std::set<core::Size> const & focus, Real distance, bool include_focus ):
	ResidueSelector(),
	focus_selector_(/* NULL */)
{
	set_defaults();
	set_focus(focus);
	set_distance(distance);
	set_include_focus_in_subset( include_focus );
}

NeighborhoodResidueSelector::NeighborhoodResidueSelector( ResidueSubset const & focus, Real distance, bool include_focus ):
	ResidueSelector(),
	focus_selector_(/* NULL */)
{
	set_defaults();
	set_focus( focus );
	set_distance( distance );
	set_include_focus_in_subset( include_focus );

}




/// @brief Copy constructor
///
NeighborhoodResidueSelector::NeighborhoodResidueSelector( NeighborhoodResidueSelector const &src) :
	focus_( src.focus_ ),
	focus_str_( src.focus_str_ ),
	distance_( src.distance_ ),
	focus_selector_( src.focus_selector_ ),
	include_focus_in_subset_(src.include_focus_in_subset_)
{}

NeighborhoodResidueSelector::~NeighborhoodResidueSelector() {}

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
ResidueSelectorOP NeighborhoodResidueSelector::clone() const { return ResidueSelectorOP( new NeighborhoodResidueSelector(*this) ); }

void
NeighborhoodResidueSelector::set_defaults() {
	distance_ = 10.0;
	include_focus_in_subset_ = true;
	focus_str_ = "";

}


void
NeighborhoodResidueSelector::set_focus( std::set<Size> const &focus )
{
	clear_focus();
	focus_ = focus;
}

void
NeighborhoodResidueSelector::set_focus( utility::vector1< bool > const & focus)
{
	clear_focus();
	for ( core::Size i = 1; i <= focus.size(); ++i ) {
		if ( focus[ i ] ) focus_.insert( i );
	}

}

void
NeighborhoodResidueSelector::set_focus( std::string const &focus_str )
{
	clear_focus();
	focus_str_ = focus_str;
}

void NeighborhoodResidueSelector::set_focus_selector( ResidueSelectorCOP rs )
{
	focus_selector_ = rs;
}

void
NeighborhoodResidueSelector::set_distance( Real distance )
{
	distance_ = distance;
}

void
NeighborhoodResidueSelector::set_include_focus_in_subset(bool include_focus){
	include_focus_in_subset_ = include_focus;
}

void
NeighborhoodResidueSelector::clear_focus(){
	focus_str_ = "";
	focus_.clear();
	focus_selector_ = NULL;
}

void
NeighborhoodResidueSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap)
{
	if ( tag->hasOption("selector") ) {

		// grab the ResidueSelector from the selector option
		// and then grab each of the indicated residue selectors from the datamap.
		std::string selector_str;
		selector_str = tag->getOption< std::string >( "selector" );


		ResidueSelectorCOP selector = datamap.get_ptr< ResidueSelector const >( "ResidueSelector", selector_str );
		set_focus_selector( selector );

	} else if ( tag->hasOption("resnums") ) {
		set_focus ( tag->getOption< std::string >( "resnums" ) );
	} else {
		throw utility::excn::EXCN_Msg_Exception("You must provide either resnums OR selector to give the focus residues.");
	}

	set_distance( tag->getOption< Real >( "distance", distance_ ) );
	set_include_focus_in_subset( tag->getOption< bool >( "include_focus_in_subset", include_focus_in_subset_));

}

void
NeighborhoodResidueSelector::get_focus(
	core::pose::Pose const & pose,
	ResidueSubset & subset,
	ResidueSubset & focus
) const
{
	debug_assert(pose.total_residue() == subset.size());
	debug_assert(pose.total_residue() == focus.size());

	bool focus_set = false;
	if ( focus_selector_ ) {
		focus = focus_selector_->apply( pose );
		focus_set = true;
	} else if ( focus_.size() > 0 ) {
		for ( std::set< Size >::const_iterator it = focus_.begin();
				it != focus_.end(); ++it ) {

			focus[ *it ] = true;
			focus_set = true;

		}
	} else {

		std::set< Size > const res_vec( get_resnum_list( focus_str_, pose ) );
		for ( std::set< Size >::const_iterator it = res_vec.begin();
				it != res_vec.end(); ++it ) {

			focus[ *it ] = true;
			focus_set = true;
		}
	}

	if ( include_focus_in_subset_ ) {
		subset = focus;

	}
	if ( ! focus_set ) {
		throw utility::excn::EXCN_Msg_Exception("Focus not set for NeighborhoodResidueSelector.  A focus must be set!");
	}
}

ResidueSubset
NeighborhoodResidueSelector::apply( core::pose::Pose const & pose ) const
{


	ResidueSubset subset( pose.total_residue(), false );
	ResidueSubset focus_subset( pose.total_residue(), false);

	// set subset to focus if option is true.  Parse focus from string, or obtain from residue selector.
	get_focus(pose, subset, focus_subset);

	debug_assert( focus_subset.size() > 0 );

	utility::vector1< Size > focus_residues = get_residues_from_subset(focus_subset);
	if ( distance_ > 10.0 ) {
		Real const dst_squared = distance_ * distance_;
		// go through each residue of the pose and check if it's near anything in the focus set
		for ( Size ii = 1; ii < pose.total_residue() ; ++ii ) {
			if ( subset[ ii ] ) continue;
			conformation::Residue const & r1( pose.residue( ii ) );

			for ( core::Size focus_res = 1; focus_res <= focus_residues.size(); ++focus_res ) {
				conformation::Residue const & r2( pose.residue( focus_res ) );
				Real const d_sq( r1.xyz( r1.nbr_atom() ).distance_squared( r2.xyz( r2.nbr_atom() ) ) );
				if ( d_sq <= dst_squared ) {
					subset[ ii ] = true;
				}
			} // focus set
		} // subset
	} else {
		if ( include_focus_in_subset_ ) {
			fill_neighbor_residues( pose, subset, distance_);
		} else {
			subset = get_neighbor_residues(pose, focus_subset, distance_);
		}
	}

	return subset;
}




std::string NeighborhoodResidueSelector::get_name() const {
	return NeighborhoodResidueSelector::class_name();
}

std::string NeighborhoodResidueSelector::class_name() {
	return "Neighborhood";
}

void
NeighborhoodResidueSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attributes;
	attributes + XMLSchemaAttribute( "selector", xs_string        )
		+ XMLSchemaAttribute( "resnums",  "int_cslist"     )
		+ XMLSchemaAttribute::required_attribute( "distance", xs_decimal );
	xsd_type_definition_w_attributes_and_optional_subselector( xsd, class_name(), attributes );
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
NeighborhoodResidueSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	NeighborhoodResidueSelector::provide_xml_schema( xsd );
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
	arc( CEREAL_NVP( include_focus_in_subset_ ) ); //bool
	arc( CEREAL_NVP( focus_selector_ ) ); // ResidueSelectorCOP
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
	arc( include_focus_in_subset_ ); //bool
	focus_selector_ = local_focus_selector; // copy the non-const pointer(s) into the const pointer(s)
}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::NeighborhoodResidueSelector );
CEREAL_REGISTER_TYPE( core::select::residue_selector::NeighborhoodResidueSelector )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_task_residue_selector_NeighborhoodResidueSelector )
#endif // SERIALIZATION
