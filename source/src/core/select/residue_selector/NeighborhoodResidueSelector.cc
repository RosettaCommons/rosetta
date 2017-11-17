// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/NeighborhoodResidueSelector.hh
/// @brief  The NeighborhoodResidueSelector selects residues in a given proximity of set focus residues
/// @author Robert Lindner (rlindner@mpimf-heidelberg.mpg.de)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com) - 10 A neighbor graph, simplification, ResidueSubset as focus, clean up, etc.
/// @author Gerard Daniel (gerardda@uw.edu) added support for using custom atom names for distance measure as an alternative to using neighbor atom. This should be especially useful when selecting around a ligand.

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
#include <core/scoring/Energies.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// C++ headers
#include <utility/assert.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>


// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/string.hpp>



#endif // SERIALIZATION


static basic::Tracer TR("core.select.residue_selector.NeighborhoodResidueSelector");

namespace core {
namespace select {
namespace residue_selector {

NeighborhoodResidueSelector::NeighborhoodResidueSelector():
	ResidueSelector(),
	focus_selector_(nullptr)
{
	set_defaults();
}

NeighborhoodResidueSelector::NeighborhoodResidueSelector( ResidueSubset const & focus, Real distance, bool include_focus ):
	ResidueSelector(),
	focus_selector_(nullptr)
{
	set_defaults();
	set_focus( focus );
	set_distance( distance );
	set_include_focus_in_subset( include_focus );

}

NeighborhoodResidueSelector::NeighborhoodResidueSelector( ResidueSelectorCOP selector, Real distance, bool include_focus ):
	ResidueSelector(),
	focus_selector_(selector)
{
	set_defaults();
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
{
	if ( src.focus_selector_ ) focus_selector_ = src.focus_selector_->clone();
}

NeighborhoodResidueSelector::~NeighborhoodResidueSelector() {}

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
ResidueSelectorOP NeighborhoodResidueSelector::clone() const { return ResidueSelectorOP( new NeighborhoodResidueSelector(*this) ); }

void
NeighborhoodResidueSelector::set_defaults() {
	distance_ = 10.0;
	include_focus_in_subset_ = true;

}


void
NeighborhoodResidueSelector::set_focus( utility::vector1< bool > const & focus)
{
	clear_focus();
	focus_ = focus;
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

/// @brief setter for custom atom names
/// @details A list of atom names, for each focus residue, the positions of which will be used for measuring distance to find neighbors. Since focus residues will be known only during the apply time, a check to see if the number of given atom names is equal to the number of focus residues. If these are not equal, an error is thrown.
void
NeighborhoodResidueSelector::set_atom_names_for_distance_measure( utility::vector1< std::string > const & atom_names ) {
	atom_names_for_distance_measure_ = atom_names;
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
	} else if ( tag->size() > 1 ) {
		//Check for embedded ResidueSelector as focus.
		utility::vector0< utility::tag::TagCOP > const & tags = tag->getTags();
		if ( tags.size() > 1 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "NeighborhoodResidueSelector takes at most one ResidueSelector to determine the focus!\n" );
		}
		ResidueSelectorCOP rs = get_embedded_residue_selector( tag, datamap );
		set_focus_selector( rs );
	} else {
		throw CREATE_EXCEPTION(utility::excn::Exception, "You must provide either resnums or selector tag, or an embedded selector to give the focus residues.");
	}

	set_distance( tag->getOption< Real >( "distance", distance_ ) );
	set_include_focus_in_subset( tag->getOption< bool >( "include_focus_in_subset", include_focus_in_subset_));

	if ( tag->hasOption("atom_names_for_distance_measure") ) {
		set_atom_names_for_distance_measure( utility::string_split( tag->getOption< std::string >( "atom_names_for_distance_measure" ), ',') );
		TR.Warning << "Will use given atom names instead of residue neighbor atoms. The number of specified atoms should be equal to the number of focus residues else an error will be thrown." << std::endl;
	}

}

void
NeighborhoodResidueSelector::get_focus(
	core::pose::Pose const & pose,
	ResidueSubset & subset,
	ResidueSubset & focus
) const
{

	debug_assert(pose.size() == subset.size());
	debug_assert(pose.size() == focus.size());

	bool focus_set = false;
	if ( focus_selector_ ) {
		focus = focus_selector_->apply( pose );
		focus_set = true;
	} else if ( focus_.size() > 0 ) {
		focus = focus_;
		focus_set = true;
	} else if ( focus_str_ != "" ) {
		std::set< Size > const res_vec( get_resnum_list( focus_str_, pose ) );
		for ( std::set< Size >::const_iterator it = res_vec.begin();
				it != res_vec.end(); ++it ) {

			focus[ *it ] = true;
			focus_set = true;
		}
	}

	if ( ! focus_set ) {
		throw CREATE_EXCEPTION(utility::excn::Exception, "Focus not set for NeighborhoodResidueSelector.  A focus must be set!");
	}


	if ( include_focus_in_subset_ ) {
		subset = focus;

	}

}

ResidueSubset
NeighborhoodResidueSelector::apply( core::pose::Pose const & pose ) const
{


	ResidueSubset subset( pose.size(), false );
	ResidueSubset focus_subset( pose.size(), false);

	// set subset to focus if option is true.  Parse focus from string, or obtain from residue selector.
	get_focus(pose, subset, focus_subset);

	debug_assert( focus_subset.size() > 0 );

	utility::vector1< Size > focus_residues = get_residues_from_subset(focus_subset);
	if ( focus_residues.size() == pose.size() ) {
		return subset;
	}

	if ( distance_ > 10.0 || atom_names_for_distance_measure_.size() ) {

		// if custom atoms are given, check to see whether their count is equal to the count of focus residues
		if ( atom_names_for_distance_measure_.size() && atom_names_for_distance_measure_.size() != focus_residues.size() ) {
			utility_exit_with_message( "The number of atom names specified is not equal to the number focus residues!" );
		}

		Real const dst_squared = distance_ * distance_;

		// go through each residue of the pose and check if it's near anything in the focus set
		for ( Size ii = 1; ii <= pose.size() ; ++ii ) {
			if ( subset[ ii ] ) continue;
			if ( ! include_focus_in_subset_ && focus_subset[ ii ] ) continue;
			conformation::Residue const & r1( pose.residue( ii ) );

			core::Size i_atom_names = 0;
			for ( core::Size focus_res : focus_residues ) {
				i_atom_names++;
				if ( focus_res == ii ) {
					if ( include_focus_in_subset_ && ! subset[focus_res] ) {
						subset[focus_res] = true;
					} else {
						continue;
					}
				}

				conformation::Residue const & r2( pose.residue( focus_res ) );

				// if atom names for focus residues are given, use those instead of neighbor atoms
				core::Vector focus_residue_atom_xyz;
				if ( atom_names_for_distance_measure_.size() ) {
					// exits with an error, if atom name is not found in the residue
					core::Size atom_index = r2.atom_index( atom_names_for_distance_measure_[i_atom_names] );
					focus_residue_atom_xyz = r2.xyz( atom_index );
					TR << "Using atom " << atom_names_for_distance_measure_[i_atom_names] << " for residue "
						<< r2.name3() << " to find neighbors." << std::endl;
				} else {
					focus_residue_atom_xyz = r2.xyz( r2.nbr_atom() );
				}

				Real const d_sq( r1.xyz( r1.nbr_atom() ).distance_squared( focus_residue_atom_xyz ) );

				if ( d_sq <= dst_squared ) {
					subset[ ii ] = true;
				}
			} // focus set
		} // subset
	} else {
		// The neighbor graph must be updated in order for these methods to work.
		// Since pose is const, we must clone it if the neighbor graph is bad.
		// Initializing references is a bit tricky which is why this comment is here. -bcov
		bool using_clone = ! pose.energies().residue_neighbors_updated();
		pose::Pose pose_clone;
		if ( using_clone ) {
			pose_clone = pose;
			pose_clone.update_residue_neighbors();
			if ( TR.Warning.visible() ) {
				TR.Warning << "################ Cloning pose and building neighbor graph ################" << std::endl;
				TR.Warning << "Ensure that pose is either scored or has update_residue_neighbors() called" << std::endl;
				TR.Warning << "before using NeighborhoodResidueSelector for maximum performance!" << std::endl;
				TR.Warning << "##########################################################################" << std::endl;
			}
		}
		const pose::Pose & pose_to_use = using_clone ? pose_clone : pose;

		if ( include_focus_in_subset_ ) {
			fill_neighbor_residues( pose_to_use, subset, distance_);
		} else {
			subset = get_neighbor_residues(pose_to_use, focus_subset, distance_);
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
	attributes + XMLSchemaAttribute( "selector", xs_string        , "XRW TO DO" )
		+ XMLSchemaAttribute( "resnums",  xsct_int_cslist     , "XRW TO DO" )
		+ XMLSchemaAttribute( "distance", xsct_real , "XRW TO DO" )
		+ XMLSchemaAttribute("atom_names_for_distance_measure", xs_string, "A list of comma separated atom names, for each focus residue, the positions of which will be used for measuring distance to find neighbors." )
		+ XMLSchemaAttribute::attribute_w_default("include_focus_in_subset", xsct_rosetta_bool, "XRW TO DO", "true");
	xsd_type_definition_w_attributes_and_optional_subselector( xsd, class_name(),"XRW TO DO", attributes );
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
	arc( CEREAL_NVP( focus_ ) ); // utility::vector1<bool>
	arc( CEREAL_NVP( focus_str_ ) ); // std::string
	arc( CEREAL_NVP( distance_ ) ); // Real
	arc( CEREAL_NVP( include_focus_in_subset_ ) ); //bool
	arc( CEREAL_NVP( focus_selector_ ) ); // ResidueSelectorCOP
	arc( CEREAL_NVP( atom_names_for_distance_measure_ ) ); // utility::vector1< std::string >
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::select::residue_selector::NeighborhoodResidueSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( focus_ ); // std::set<Size>
	arc( focus_str_ ); // std::string
	arc( distance_ ); // Real
	arc( include_focus_in_subset_ ); //bool
	std::shared_ptr< core::select::residue_selector::ResidueSelector > local_focus_selector;
	arc( local_focus_selector ); // ResidueSelectorCOP
	focus_selector_ = local_focus_selector; // copy the non-const pointer(s) into the const pointer(s)
	arc( atom_names_for_distance_measure_ ); // utility::vector1< std::string >
}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::NeighborhoodResidueSelector );
CEREAL_REGISTER_TYPE( core::select::residue_selector::NeighborhoodResidueSelector )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_task_residue_selector_NeighborhoodResidueSelector )
#endif // SERIALIZATION
