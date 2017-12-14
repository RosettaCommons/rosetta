// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/residue_selector/RandomGlycanFoliageSelector.hh
/// @brief  Selects a random carbohydrate residue from a subset or selector, then selects the rest of the glycan foliage.  Used for sampling.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <core/select/residue_selector/RandomGlycanFoliageSelector.hh>
#include <core/select/residue_selector/GlycanResidueSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>

#include <core/select/residue_selector/util.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh> // for xml schema utility functions
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <numeric/random/random.hh>
#include <utility/assert.hh>

static basic::Tracer TR( "core.select.residue_selector.RandomGlycanFoliageSelector" );

namespace core {
namespace select {
namespace residue_selector {

/// @brief Constructor.
///
RandomGlycanFoliageSelector::RandomGlycanFoliageSelector():
	ResidueSelector(),
	selector_(nullptr)
{
}

RandomGlycanFoliageSelector::RandomGlycanFoliageSelector( ResidueSubset const & subset):
	ResidueSelector(),
	selector_(nullptr)
{
	subset_ = subset;
}

RandomGlycanFoliageSelector::RandomGlycanFoliageSelector( ResidueSelectorOP selector):
	ResidueSelector()
{
	selector_ = selector;
}



/// @brief Destructor.
///
RandomGlycanFoliageSelector::~RandomGlycanFoliageSelector() = default;

/// @brief Copy Constructor.  Usually not necessary unless you need deep copying (e.g. OPs)
RandomGlycanFoliageSelector::RandomGlycanFoliageSelector(RandomGlycanFoliageSelector const & src):
	ResidueSelector( src ),
	subset_(src.subset_)
{
	if ( src.selector_ ) selector_ = src.selector_->clone();
}

/// @brief Set a subset to select the glycan root and subsequent foliage on.
void
RandomGlycanFoliageSelector::set_subset(ResidueSubset const & subset){
	selector_ = nullptr;
	subset_ = subset;
}

void
/// @brief Set a selector to set the glycan root and subsequent foliage on.
RandomGlycanFoliageSelector::set_selector( ResidueSelectorCOP selector){
	selector_ = selector;
}


/// @brief Clone function.
/// @details Copy this object and return owning pointer to the copy (created on the heap).
RandomGlycanFoliageSelector::ResidueSelectorOP
RandomGlycanFoliageSelector::clone() const {
	return ResidueSelectorOP( new RandomGlycanFoliageSelector(*this) );
}



RandomGlycanFoliageSelector::ResidueSelectorOP
RandomGlycanFoliageSelectorCreator::create_residue_selector() const {
	return RandomGlycanFoliageSelector::ResidueSelectorOP( new RandomGlycanFoliageSelector );
}

std::string RandomGlycanFoliageSelector::get_name() const
{
	return RandomGlycanFoliageSelector::class_name();
}

std::string RandomGlycanFoliageSelector::class_name()
{
	return "RandomGlycanFoliageSelector";
}

std::string
RandomGlycanFoliageSelectorCreator::keyname() const {
	return RandomGlycanFoliageSelector::class_name();
}

/// @brief XML parse.
/// @details Parse RosettaScripts tags and set up this mover.
void
RandomGlycanFoliageSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data)
{
	if ( tag->hasOption("residue_selector") ) {
		selector_ = parse_residue_selector(tag, data);
	}

}

void RandomGlycanFoliageSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd)
{

	using namespace utility::tag;
	AttributeList attributes;
	attributes_for_parse_residue_selector(attributes);
	xsd_type_definition_w_attributes_and_optional_subselector( xsd, class_name(),
		"Selects a random carbohydrate residue from a subset or selector, then selects the rest of the glycan foliage.  Used for sampling.",
		attributes );

}

/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
///
void
RandomGlycanFoliageSelectorCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	RandomGlycanFoliageSelector::provide_xml_schema( xsd );
}

RandomGlycanFoliageSelector::ResidueSubset
RandomGlycanFoliageSelector::apply(
	core::pose::Pose const & pose
) const {

	TR.Debug << "applying RandomGlycanFoliageSelector" << std::endl;
	GlycanResidueSelector glycan_selector = GlycanResidueSelector();
	glycan_selector.set_include_root( true ); //We want to be able to select the residue and the rest, and we won't really usually have the ASN here.

	ResidueSubset select_on_subset;
	utility::vector1< core::Size > select_on_vector;

	if ( selector_ ) {
		select_on_subset = selector_->apply( pose );
	} else if ( subset_.empty() ) {
		//Use ALL glycan residues
		select_on_subset = glycan_selector.apply(pose);
	} else {
		select_on_subset = subset_;
	}

	//Group active residues, make sure we are carbohydrate.
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( select_on_subset[i] && pose.residue( i ).is_carbohydrate() ) {
			select_on_vector.push_back( i );
		}
	}

	core::Size index = numeric::random::rg().random_range( 1, select_on_vector.size() );
	core::Size resnum = select_on_vector[ index ];
	glycan_selector.set_select_from_branch_residue(resnum);
	utility::vector1< bool > subset = glycan_selector.apply(pose);
	return subset;

}


} //core
} //select
} //residue_selector





