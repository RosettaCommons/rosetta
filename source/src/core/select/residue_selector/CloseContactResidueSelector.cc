// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/residue_selector/CloseContactResidueSelector.hh
/// @brief  A class that finds the neighboring residues for a particular residue by looking at atom-atom distances
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/select/residue_selector/CloseContactResidueSelector.hh>
#include <core/select/residue_selector/CloseContactResidueSelectorCreator.hh>

#include <core/select/residue_selector/util.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh> // for xml schema utility functions

// Project headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <utility/assert.hh>

static THREAD_LOCAL basic::Tracer TR( "core.select.residue_selector.CloseContactResidueSelector" );

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {

/// @brief Constructor.
///
CloseContactResidueSelector::CloseContactResidueSelector():
	ResidueSelector()
{
}

/// @brief Destructor.
///
CloseContactResidueSelector::~CloseContactResidueSelector() {}

/// @brief Copy Constructor.  Usually not necessary unless you need deep copying (e.g. OPs)
CloseContactResidueSelector::CloseContactResidueSelector(CloseContactResidueSelector const & ) = default;

/// @brief Clone function.
/// @details Copy this object and return owning pointer to the copy (created on the heap).
CloseContactResidueSelector::ResidueSelectorOP
CloseContactResidueSelector::clone() const {
	return ResidueSelectorOP( new CloseContactResidueSelector(*this) );
}

/// @brief "Apply" function.
/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
/// indicating whether each residue is selected ("true") or not ("false").
CloseContactResidueSelector::ResidueSubset
CloseContactResidueSelector::apply(
	core::pose::Pose const & pose
) const {
	ResidueSubset subset( pose.total_residue(), false );
	ResidueSubset central_residues = central_residues_selector_->apply( pose );

	Size count_central = 0;
	for ( auto res : central_residues ) if ( res ) ++count_central;

	core::Real cut2 = close_contact_threshold_ * close_contact_threshold_;

	if ( count_central == 1 ) {
		core::Size central_residue( 0 );
		for ( core::Size ii = 1; ii <= central_residues.size(); ++ii ) {
			if ( central_residues[ ii ] ) {
				central_residue = ii;
				break;
			}
		}
		debug_assert( central_residue != 0 );

		subset[ central_residue ] = true;

		core::conformation::Residue const & central_rsd( pose.residue( central_residue  ) );
		core::Vector central_nbr_atom = central_rsd.xyz( central_rsd.nbr_atom() );

		TR << "Selecting residues around central seqpos: " << central_residue << " with a distance cutoff of " << close_contact_threshold_ << std::endl;

		for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			if ( ii == central_residue ) continue;

			core::conformation::Residue const & iirsd( pose.residue(ii) );
			Real cutoff = central_rsd.nbr_radius() + close_contact_threshold_ + iirsd.nbr_radius();
			if ( iirsd.xyz( iirsd.nbr_atom() ).distance_squared( central_nbr_atom ) > cutoff * cutoff ) continue; // too far

			for ( core::Size jj = 1; jj <= iirsd.natoms(); ++jj ) {
				for ( core::Size kk = 1; kk <= central_rsd.natoms(); ++kk ) {
					if ( iirsd.xyz( jj ).distance_squared( central_rsd.xyz( kk ) ) < cut2 ) {
						subset[ ii ] = true;
						break;
					}
				}
				if ( subset[ ii ] ) break;
			}
		}
	} else {
		subset = central_residues;
		utility::vector1< core::Size > central_res_ids;
		central_res_ids.reserve( count_central );
		for ( core::Size ii = 1; ii <= central_residues.size(); ++ii ) {
			if ( central_residues[ ii ] ) {
				central_res_ids.push_back( ii );
			}
		}

		for ( Size ii = 1; ii <= subset.size(); ++ii ) {
			if ( subset[ ii ] ) continue;
			conformation::Residue const & ii_rsd( pose.residue( ii ) );
			Vector ii_nbr_xyz = ii_rsd.xyz( ii_rsd.nbr_atom() );
			for ( Size jj : central_res_ids ) {
				conformation::Residue const & jj_rsd( pose.residue( jj ) );
				Vector jj_nbr_xyz = jj_rsd.xyz( jj_rsd.nbr_atom() );
				Real cutoff = ii_rsd.nbr_radius() + jj_rsd.nbr_radius() + close_contact_threshold_;
				if ( ii_nbr_xyz.distance_squared( jj_nbr_xyz ) > cutoff * cutoff  ) continue;
				for ( Size kk = 1; kk <= ii_rsd.natoms(); ++kk ) {
					for ( Size ll = 1; ll <= jj_rsd.natoms(); ++ll ) {
						if ( ii_rsd.xyz( kk ).distance_squared( jj_rsd.xyz( ll )) <= cut2 ) {
							subset[ ii ] = true;
							break;
						}
					}
					if ( subset[ ii ] ) break;
				}
				if ( subset[ ii ] ) break;
			}
		}
	}
	return subset;
}

void CloseContactResidueSelector::central_residue_group_selector( ResidueSelectorCOP selector )
{
	central_residues_selector_ = selector;
}

void CloseContactResidueSelector::threshold( core::Real contact_threshold )
{
	close_contact_threshold_ = contact_threshold;
}

core::Real CloseContactResidueSelector::threshold() const
{
	return close_contact_threshold_;
}

/// @brief XML parse.
/// @details Parse RosettaScripts tags and set up this mover.
void
CloseContactResidueSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
)
{
	ResidueSelectorCOP selector;
	if ( tag->hasOption("residue_selector") ) {
		std::string selector_name = tag->getOption< std::string >( "residue_selector" );
		try {
			selector = datamap.get_ptr< ResidueSelector const >( "ResidueSelector", selector_name );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::stringstream error_msg;
			error_msg << "Failed to find ResidueSelector named '" << selector_name << "' from the Datamap from AndResidueSelector::parse_my_tag.\n";
			error_msg << e.msg();
			throw utility::excn::EXCN_Msg_Exception( error_msg.str() );
		}
	} else {
		// add selectors from tags
		auto subtags = tag->getTags();
		if ( subtags.size() != 1 ) {
			std::ostringstream error_msg;
			if ( subtags.size() == 0 ) {
				error_msg << "If you do not provide a \"residue_selector\" attribute, then you must provide a ResidueSelector as a sub-element\n";
			} else {
				error_msg << "You may only provide a single sub-element of the " << class_name() << " residue selector\n";
			}
			throw utility::excn::EXCN_Msg_Exception( error_msg.str() );
		}
		selector = ResidueSelectorFactory::get_instance()->new_residue_selector(
			subtags[0]->getName(),
			subtags[0],
			datamap
		);
	}

	central_residue_group_selector( selector );
	threshold( tag->getOption< core::Real >( "contact_threshold", 4.5 ) );

}

std::string CloseContactResidueSelector::get_name() const
{
	return CloseContactResidueSelector::class_name();
}

std::string CloseContactResidueSelector::class_name()
{
	return "CloseContact";
}

void CloseContactResidueSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute( "residue_selector", xs_string, "The name of the residue selector that's used to identify"
		" the central residue(s); other residues are selected if they are in close contact with this group of central residues"  )
		+ XMLSchemaAttribute::attribute_w_default( "contact_threshold", xsct_real, "The distance, in Angstroms, around the"
		" residues in the central-residues set which defines the atomic-contact cutoff", "4.5" );

	xsd_type_definition_w_attributes_and_optional_subselector( xsd, class_name(), "This class selects all residues in close contact"
		" with a central residue, or a set of central residues, where this decision is based on atom-to-atom distances for the current"
		" set of coordinates in the Pose", attributes );

}

CloseContactResidueSelector::ResidueSelectorOP
CloseContactResidueSelectorCreator::create_residue_selector() const {
	return CloseContactResidueSelector::ResidueSelectorOP( new CloseContactResidueSelector );
}

std::string
CloseContactResidueSelectorCreator::keyname() const {
	return CloseContactResidueSelector::class_name();
}

/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
///
void
CloseContactResidueSelectorCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	CloseContactResidueSelector::provide_xml_schema( xsd );
}





} //core
} //select
} //residue_selector

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::select::residue_selector::CloseContactResidueSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( close_contact_threshold_ ) ); // core::Real
	arc( CEREAL_NVP( central_residues_selector_ ) ); // ResidueSelectorCOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::select::residue_selector::CloseContactResidueSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( close_contact_threshold_ ); // core::Real
	std::shared_ptr< core::select::residue_selector::ResidueSelector > local_central_residues_selector;
	arc( local_central_residues_selector ); // ResidueSelectorCOP
	central_residues_selector_ = local_central_residues_selector; // copy the non-const pointer(s) into the const pointer(s)
}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::CloseContactResidueSelector );
CEREAL_REGISTER_TYPE( core::select::residue_selector::CloseContactResidueSelector )

CEREAL_REGISTER_DYNAMIC_INIT( core_select_residue_selector_CloseContactResidueSelector )
#endif // SERIALIZATION
