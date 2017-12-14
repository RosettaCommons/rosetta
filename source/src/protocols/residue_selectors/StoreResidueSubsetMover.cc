// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/residue_selectors/StoreResidueSubsetMover.cc
/// @brief The StoreResidueSubsetMover allows you to create a at some point during a
/// RosettaScripts run and save it for access later during the same run. Can be useful
/// for mutating/analyzing a particular set of residues using many different movers/filters.
/// @author Tom Linsky ( tlinsky at uw dot edu )

// Unit Headers
#include <protocols/residue_selectors/StoreResidueSubsetMover.hh>
#include <protocols/residue_selectors/StoreResidueSubsetMoverCreator.hh>

// Protocol Headers
#include <core/select/residue_selector/util.hh>
#include <protocols/rosetta_scripts/util.hh>
/*#include <protocols/toolbox/task_operations/STMStoredTask.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/pose/symmetry/util.hh>
#include <basic/datacache/CacheableData.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/make_symmetric_task.hh>
*/
// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/select/residue_selector/CachedResidueSubset.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>

// Utility Headers
#include <basic/datacache/BasicDataCache.hh>
#include <basic/Tracer.hh>
#include <utility>
#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// C++ Headers

// ObjexxFCL Headers

static basic::Tracer TR( "protocols.residue_selectors.StoreResidueSubsetMover" );

namespace protocols {
namespace residue_selectors {

// @brief default constructor
StoreResidueSubsetMover::StoreResidueSubsetMover():
	protocols::moves::Mover(),
	selector_( new core::select::residue_selector::TrueResidueSelector ),
	subset_name_( "" ),
	overwrite_( false )
{}

StoreResidueSubsetMover::StoreResidueSubsetMover(
	core::select::residue_selector::ResidueSelectorCOP selector,
	std::string const & subset_name,
	bool const overwrite ):
	protocols::moves::Mover(),
	selector_(std::move( selector )),
	subset_name_( subset_name ),
	overwrite_( overwrite )
{}

// @brief destructor
StoreResidueSubsetMover::~StoreResidueSubsetMover() = default;

void
StoreResidueSubsetMover::apply( core::pose::Pose & pose )
{
	runtime_assert( selector_ || subset_ );
	if ( subset_name_.empty() ) {
		utility_exit_with_message( "No subset_name specified to StoreResidueSubset mover.  You must specify one." );
	}
	if ( selector_ ) {
		subset_.reset( new core::select::residue_selector::ResidueSubset( selector_->apply( pose ) ) );
	}
	//core::select::residue_selector::ResidueSubsetCOP const subset( new core::select::residue_selector::ResidueSubset( selector_->apply( pose ) ) );
	/*
	if ( core::pose::symmetry::is_symmetric(pose) ) {
	core::pack::make_symmetric_PackerTask_by_truncation(pose, task); // Does this need to be fixed or omitted?
	}
	*/
	// If the pose doesn't have cached subset data, add blank data
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::STORED_RESIDUE_SUBSET ) ) {
		core::select::residue_selector::CachedResidueSubsetOP blank_stored_subsets( new core::select::residue_selector::CachedResidueSubset );
		pose.data().set( core::pose::datacache::CacheableDataType::STORED_RESIDUE_SUBSET, blank_stored_subsets );
	}

	// Grab a reference to the data
	debug_assert( utility::pointer::dynamic_pointer_cast< core::select::residue_selector::CachedResidueSubset >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::STORED_RESIDUE_SUBSET ) ) );
	core::select::residue_selector::CachedResidueSubset & stored_subsets =
		*( utility::pointer::static_pointer_cast< core::select::residue_selector::CachedResidueSubset >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::STORED_RESIDUE_SUBSET ) ) );

	if ( !overwrite_ && stored_subsets.has_subset( subset_name_ ) ) {
		utility_exit_with_message( "A stored residue subset with the name " + subset_name_ + " already exists; you must set overwrite flag to true to overwrite." );
	}
	stored_subsets.set_subset( subset_, subset_name_ );
}

void
StoreResidueSubsetMover::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data_map, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{
	selector_ = protocols::rosetta_scripts::parse_residue_selector( tag, data_map );
	subset_name_ = tag->getOption< std::string >( "subset_name" );
	overwrite_ = tag->getOption< bool >( "overwrite", overwrite_ );
}

// @brief Identification
// XRW TEMP std::string
// XRW TEMP StoreResidueSubsetMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return StoreResidueSubsetMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP StoreResidueSubsetMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "StoreResidueSubset";
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP StoreResidueSubsetMover::get_name() const
// XRW TEMP {
// XRW TEMP  return "StoreResidueSubset";
// XRW TEMP }


// XRW TEMP protocols::moves::MoverOP
// XRW TEMP StoreResidueSubsetMoverCreator::create_mover() const
// XRW TEMP {
// XRW TEMP  return protocols::moves::MoverOP( new StoreResidueSubsetMover );
// XRW TEMP }

protocols::moves::MoverOP
StoreResidueSubsetMover::clone() const
{
	return protocols::moves::MoverOP( new StoreResidueSubsetMover( *this ) );
}

protocols::moves::MoverOP
StoreResidueSubsetMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new StoreResidueSubsetMover );
}

std::string StoreResidueSubsetMover::get_name() const {
	return mover_name();
}

std::string StoreResidueSubsetMover::mover_name() {
	return "StoreResidueSubset";
}

void StoreResidueSubsetMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	AttributeList attlist;
	core::select::residue_selector::attributes_for_parse_residue_selector(attlist);
	attlist + XMLSchemaAttribute::required_attribute( "subset_name", xs_string, "XRW TO DO")
		+ XMLSchemaAttribute( "overwrite", xsct_rosetta_bool, "XRW TO DO");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string StoreResidueSubsetMoverCreator::keyname() const {
	return StoreResidueSubsetMover::mover_name();
}

protocols::moves::MoverOP
StoreResidueSubsetMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new StoreResidueSubsetMover );
}

void StoreResidueSubsetMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	StoreResidueSubsetMover::provide_xml_schema( xsd );
}


} // residue_selectors
} // protocols
