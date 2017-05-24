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
#include <core/select/residue_selector/NeighborhoodResidueSelector.hh>

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
	num_residues_( 1 ),
	select_res_cluster_( false ),
	distance_cutoff_( 8.0 )
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
	num_residues_( num_residues ),
	select_res_cluster_( false ),
	distance_cutoff_( 8.0 )
{
}

RandomResidueSelector::RandomResidueSelector( ResidueSelectorCOP selector, Size const num_residues, bool const select_res_cluster, Real const distance_cutoff ):
	ResidueSelector(),
	selector_( selector ),
	num_residues_( num_residues ),
	select_res_cluster_( select_res_cluster ),
	distance_cutoff_( distance_cutoff )
{
}

RandomResidueSelector::~RandomResidueSelector() {}

ResidueSubset
RandomResidueSelector::apply( core::pose::Pose const & pose ) const
{
	ResidueVector residue_set( selector_->apply( pose ) );
	// shuffle order of residue set
	numeric::random::random_permutation( residue_set, numeric::random::rg() );
	// return random subset
	return subset_from_randomized_vector( pose, residue_set );
}

ResidueSubset
RandomResidueSelector::subset_from_randomized_vector(
	core::pose::Pose const & pose,
	ResidueVector const & random_order_residue_set ) const
{
	// calc neighbors of random residue?
	ResidueSubset nbrs( pose.size(), false );
	if ( select_res_cluster_ && random_order_residue_set.size() > 0 ) {
		// create subset that is just focus
		ResidueSubset focus( pose.size(), false );
		// focus residue is the first one in the selected set
		focus[ random_order_residue_set[ 1 ] ] = true;
		core::select::residue_selector::NeighborhoodResidueSelector nbrs_selector( focus, distance_cutoff_, true );
		// get residue subset bool vector of neighbors
		nbrs = nbrs_selector.apply( pose );
		TR << "Neigbor residues: ";
		for ( Size ires = 1; ires <= nbrs.size(); ++ires ) {
			if ( nbrs[ ires ] ) TR << ires << " ";
		}
		TR << std::endl;
	}

	TR << "Selected residues: ";

	std::set< core::Size > selected;
	// populate selected with residues from random_order_residue_set until >= num_residues_
	for ( auto const & r : random_order_residue_set ) {
		// if first res in selected is init, only add next residue r if is within cutoff distance of first
		if ( select_res_cluster_ && selected.size() >= 1 && nbrs[ r ] == false ) continue;
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
	select_res_cluster_ = tag->getOption< bool >( "select_res_cluster", select_res_cluster_ );
	distance_cutoff_ = tag->getOption< Real >( "distance_cutoff", distance_cutoff_ );

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
		"select_res_cluster", xsct_rosetta_bool,
		"option to only select multiple residues near each other, only applies to case where num_residues greater than 1, multiple random residues are required to be within distance_cutoff angstroms of one another")
		+ XMLSchemaAttribute(
		"distance_cutoff", xsct_real,
		"only active when select_res_cluster set to true, distance that defines whether two residues are neighbors or not")
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
