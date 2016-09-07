// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/residue_selectors/StoredResidueSubsetSelector.hh
/// @brief  The StoredResidueSubsetSelector selects residues using a previously stored residue subset
/// @author Tom Linsky (tlinsky@uw.edu))

// Unit headers
#include <protocols/residue_selectors/StoredResidueSubsetSelector.hh>
#include <protocols/residue_selectors/StoredResidueSubsetSelectorCreator.hh>

// Protocol Headers

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/select/residue_selector/CachedResidueSubset.hh>
#include <core/select/residue_selector/util.hh>

// Basic Headers
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// C++ headers

static THREAD_LOCAL basic::Tracer TR( "protocols.residue_selectors.StoredResidueSubsetSelector" );

namespace protocols {
namespace residue_selectors {

core::select::residue_selector::ResidueSelectorOP
StoredResidueSubsetSelectorCreator::create_residue_selector() const
{
	return core::select::residue_selector::ResidueSelectorOP( new StoredResidueSubsetSelector );
}

std::string
StoredResidueSubsetSelectorCreator::keyname() const
{
	return StoredResidueSubsetSelector::class_name();
}

void
StoredResidueSubsetSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	StoredResidueSubsetSelector::provide_xml_schema( xsd );
}

StoredResidueSubsetSelector::StoredResidueSubsetSelector() :
	ResidueSelector(),
	subset_name_( "" )
{}

StoredResidueSubsetSelector::StoredResidueSubsetSelector( std::string  subset_name ):
	ResidueSelector(),
	subset_name_(std::move( subset_name ))
{}

StoredResidueSubsetSelector::~StoredResidueSubsetSelector()
= default;

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
core::select::residue_selector::ResidueSelectorOP
StoredResidueSubsetSelector::clone() const
{
	return core::select::residue_selector::ResidueSelectorOP( new StoredResidueSubsetSelector(*this) );
}

void
quit_no_subset()
{
	std::stringstream msg;
	msg << "StoredResidueSubsetSelector: no subset name was specified!" << std::endl;
	throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
}

core::select::residue_selector::ResidueSubset
StoredResidueSubsetSelector::apply( core::pose::Pose const & pose ) const
{
	if ( subset_name_.empty() ) quit_no_subset();

	// check for presence of cached pose data
	if ( ! pose.data().has( core::pose::datacache::CacheableDataType::STORED_RESIDUE_SUBSET ) ) {
		utility_exit_with_message( "Your pose does not have CacheableData of type STORED_RESIDUE_SUBSET" );
	}

	// get cached data from pose
	debug_assert( utility::pointer::dynamic_pointer_cast< core::select::residue_selector::CachedResidueSubset const >( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::STORED_RESIDUE_SUBSET ) ) );
	core::select::residue_selector::CachedResidueSubset const & stored_subsets =
		*( utility::pointer::static_pointer_cast< core::select::residue_selector::CachedResidueSubset const >( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::STORED_RESIDUE_SUBSET ) ) );

	// check for cached subset in data
	if ( !stored_subsets.has_subset( subset_name_ ) ) {
		utility_exit_with_message( "No stored residue subset exists in the pose with the name " + subset_name_ );
	}

	core::select::residue_selector::ResidueSubsetCOP subset = stored_subsets.get_subset( subset_name_ );
	if ( !subset ) {
		utility_exit_with_message( "Stored residue subset with name " + subset_name_ + " is NULL" );
	}

	if ( pose.total_residue() != subset->size() ) {
		std::stringstream msg;
		msg << "StoredResidueSubsetSelector: Size of stored residue subset \"" << subset_name_ << "\" (" << subset->size()
			<< ") does not match pose size (" << pose.total_residue() << ")" << std::endl;
		utility_exit_with_message( msg.str() );
	}

	return *subset;
}

void
StoredResidueSubsetSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & )
{
	subset_name_ = tag->getOption< std::string >( "subset_name", subset_name_ );
	if ( subset_name_.empty() ) quit_no_subset();
}

void
StoredResidueSubsetSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attributes;
	attributes.emplace_back( "subset_name", xs_string );
	core::select::residue_selector::xsd_type_definition_w_attributes( xsd, class_name(), attributes );
}

std::string
StoredResidueSubsetSelector::get_name() const
{
	return StoredResidueSubsetSelector::class_name();
}

std::string
StoredResidueSubsetSelector::class_name()
{
	return "StoredResidueSubset";
}

} //namespace residue_selectors
} //namespace protocols

