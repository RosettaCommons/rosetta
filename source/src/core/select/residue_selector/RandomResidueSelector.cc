// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/RandomResidueSelector.hh
/// @brief  The RandomResidueSelector selects residues using a string containing residue names
/// @author Tom Linsky (tlinsky@uw.edu))

// Unit headers
#include <core/select/residue_selector/RandomResidueSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>

// Core headers
#include <core/conformation/Residue.hh>
#include <core/pose/selection.hh>
#include <core/select/residue_selector/ResidueVector.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Utility Headers
#include <numeric/random/random_permutation.hh>
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

static THREAD_LOCAL basic::Tracer TR( "core.select.residue_selector.RandomResidueSelector" );

namespace core {
namespace select {
namespace residue_selector {

RandomResidueSelector::RandomResidueSelector():
	ResidueSelector(),
	selector_( new TrueResidueSelector ),
	num_residues_( 1 )
{
}

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
ResidueSelectorOP
RandomResidueSelector::clone() const
{
	return ResidueSelectorOP( new RandomResidueSelector( *this ) );
}

RandomResidueSelector::RandomResidueSelector( ResidueSelectorCOP selector, Size const num_residues ):
	ResidueSelector(),
	selector_( selector ),
	num_residues_( num_residues )
{
}

RandomResidueSelector::~RandomResidueSelector() {}

ResidueSubset
RandomResidueSelector::apply( core::pose::Pose const & pose ) const
{
	ResidueVector residue_set( selector_->apply( pose ) );
	numeric::random::random_permutation( residue_set, numeric::random::rg() );
	return subset_from_randomized_vector( pose, residue_set );
}

ResidueSubset
RandomResidueSelector::subset_from_randomized_vector(
	core::pose::Pose const & pose,
	ResidueVector const & random_order_residue_set ) const
{
	TR << "Selected residues: ";
	std::set< core::Size > selected;
	for ( auto const & r : random_order_residue_set ) {
		TR << r << " ";
		selected.insert( r );
		if ( selected.size() >= num_residues_ ) break;
	}
	TR << std::endl;

	ResidueSubset subset( pose.size(), false );
	for ( Size const  r : selected ) {
		debug_assert( r <= subset.size() );
		subset[ r ] = true;
	}
	return subset;
}

void
RandomResidueSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data )
{
	num_residues_ = tag->getOption< Size >( "num_residues", num_residues_ );

	std::string const selectorname = tag->getOption< std::string >( "selector", "" );
	if ( !selectorname.empty() ) {
		ResidueSelectorCOP selector;
		try {
			selector = data.get_ptr< core::select::residue_selector::ResidueSelector const >( "ResidueSelector", selectorname );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::stringstream error_msg;
			error_msg << "Failed to find ResidueSelector named '" << selectorname << "' in the DataMap.\n";
			error_msg << e.msg();
			throw utility::excn::EXCN_Msg_Exception( error_msg.str() );
		}
		debug_assert( selector );
		TR << "Using residue selector " << selectorname << std::endl;
		selector_ = selector;
	}
}

std::string
RandomResidueSelector::get_name() const
{
	return RandomResidueSelector::class_name();
}

std::string
RandomResidueSelector::class_name()
{
	return "RandomResidue";
}

void
RandomResidueSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attributes;
	attributes + XMLSchemaAttribute(
		"num_residues", xsct_non_negative_integer,
		"The number of residues to be randomly selected")
		+ XMLSchemaAttribute(
		"selector", xs_string,
		"Defines the subset from which random residues are chosen.");
	xsd_type_definition_w_attributes_and_optional_subselector(
		xsd, class_name(),
		"Selects residues in the pose at random. Note that this residue "
		"selector is stochastic. This is, it will return a different set "
		"of residues every time it is called. However, the randomly selected "
		"residues can be saved using the StoreResidueSubsetMover and "
		"retrieved using the StoredResidueSubset selector.",
		attributes );
}

ResidueSelectorOP
RandomResidueSelectorCreator::create_residue_selector() const
{
	return ResidueSelectorOP( new RandomResidueSelector );
}

std::string
RandomResidueSelectorCreator::keyname() const
{
	return RandomResidueSelector::class_name();
}

void
RandomResidueSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RandomResidueSelector::provide_xml_schema( xsd );
}

} //namespace residue_selector
} //namespace select
} //namespace core


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::select::residue_selector::RandomResidueSelector::save( Archive & arc ) const
{
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( selector_ ) ); // ResidueSelectorCOP
	arc( CEREAL_NVP( num_residues_ ) ); // Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::select::residue_selector::RandomResidueSelector::load( Archive & arc )
{
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );

	std::shared_ptr< core::select::residue_selector::ResidueSelector > local_selector;
	arc( local_selector ); // ResidueSelectorCOP
	selector_ = local_selector; // copy the non-const pointer(s) into the const pointer(s)

	arc( num_residues_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::RandomResidueSelector );
CEREAL_REGISTER_TYPE( core::select::residue_selector::RandomResidueSelector )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_task_residue_selector_RandomResidueSelector )
#endif // SERIALIZATION
