// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/ByMembraneDepthSelector.hh
/// @brief  Select residues within a certain range of membrane depths
/// @author Rebecca Alford (rfalford12@gmail.com)

// Unit headers
#include <core/select/residue_selector/ByMembraneDepthSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>

// Package headers
#include <core/pose/selection.hh>
#include <core/conformation/Residue.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh>

// Project Headers
#include <core/conformation/Conformation.hh> 
#include <core/conformation/membrane/MembraneInfo.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>

#include <utility/string_util.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

#include <basic/datacache/DataMap.hh>

// C++ headers
#include <utility/assert.hh>
#include <string>

static THREAD_LOCAL basic::Tracer TR( "core.select.residue_selector.ByMembraneDepthSelector" );

namespace core {
namespace select {
namespace residue_selector {

using namespace core::select::residue_selector;

ByMembraneDepthSelector::ByMembraneDepthSelector() :
	ResidueSelector(),
	min_depth_( 0 ),
	max_depth_( 100000 )
{
	check_valid_bounds();
}

ByMembraneDepthSelector::~ByMembraneDepthSelector() {}

/// @brief Select a subset of residues based on the user-specified range of membrane depths (from the center)
/// @note This means you will always extrapolate two subsets of residues - above and below the membrane plane, because
/// the chemistry there is the same
ResidueSubset
ByMembraneDepthSelector::apply( core::pose::Pose const & pose ) const
{

	using namespace core::conformation::membrane;
	
	// Check that the pose is a membrane protein
	if ( !pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Cannot select TM spans on a non-membrane pose or pose without any transmembrane spanning segments" );
	}

	// Loop through every residue in the pose: If the residue falls within
	// the specified range, set true
	ResidueSubset selected_residues( pose.total_residue(), false );
	for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		if ( residue_within_depth( pose, ii ) ) selected_residues[ii] = true;
	}
	return selected_residues;
}

/// @brief Retrieve the minimum membrane depth criteria used to select residues
Real
ByMembraneDepthSelector::minimum_depth() const {
	return min_depth_;
}

/// @brief Modify the minimum membrane depth criteria used to select residues
void
ByMembraneDepthSelector::minimum_depth( core::Real minimum_depth ) {
	min_depth_ = minimum_depth;
}

/// @brief Retrieve the maximum memrbane depth criteria used to select residues
Real
ByMembraneDepthSelector::maximum_depth() const {
	return max_depth_;
}

/// @brief Modify the maximum membrane depth criteria used to select residues
void
ByMembraneDepthSelector::maximum_depth( core::Real maximum_depth ) {
	max_depth_ = maximum_depth;
}

bool
ByMembraneDepthSelector::residue_within_depth( core::pose::Pose const & pose, core::Size rsdnum ) const {
	
	using namespace core::conformation::membrane;
	
	// Bounding is done by [mim, max)
	core::Real rsd_z( pose.conformation().membrane_info()->residue_z_position( pose.conformation(), rsdnum ) );
	if ( std::abs( rsd_z ) >= minimum_depth() && std::abs( rsd_z ) < maximum_depth() ) {
		return true;
	}
	return false;

}

/// @brief Check that the user specified bounds are within a reasonable range
void
ByMembraneDepthSelector::check_valid_bounds() const {

	// Check that user specified bounding is valid!
	if ( min_depth_ > max_depth_ ) {
		utility_exit_with_message( "ByMembraneResidueSelector: Minimum depth selection is greater than maximum depth selection. Invalid and exiting!" );
	}
}

void
ByMembraneDepthSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & ) {

	// Read in minimum and maximum bounds from the user
	if ( tag->hasOption( "mim_depth" ) ) {
		minimum_depth( tag->getOption< core::Real >( "min_depth" ) );
	}
	
	if ( tag->hasOption( "min_depth" ) ) {
		maximum_depth( tag->getOption< core::Real >( "max_depth" ) );
	}

	// Check valid!
	check_valid_bounds();
}

std::string ByMembraneDepthSelector::get_name() const {
	return ByMembraneDepthSelector::class_name();
}

std::string ByMembraneDepthSelector::class_name() {
	return "ByMembraneDepthSelector";
}

ResidueSelectorOP
ByMembraneDepthSelector::clone() const {
	return ResidueSelectorOP( new ByMembraneDepthSelector(*this) );
}

ResidueSelectorOP
ByMembraneDepthSelectorCreator::create_residue_selector() const {
	return ResidueSelectorOP( new ByMembraneDepthSelector );
}

std::string
ByMembraneDepthSelectorCreator::keyname() const {
	return ByMembraneDepthSelector::class_name();
}

void
ByMembraneDepthSelectorCreator::provide_selector_xsd( utility::tag::XMLSchemaDefinition & xsd ) const {
	ByMembraneDepthSelector::provide_selector_xsd( xsd );
}


void
ByMembraneDepthSelector::provide_selector_xsd( utility::tag::XMLSchemaDefinition & xsd ) {

	using namespace utility::tag;
	AttributeList attributes;
	attributes.push_back( XMLSchemaAttribute( "min_depth",  xs_decimal ));
	attributes.push_back( XMLSchemaAttribute( "max_depth", xs_decimal ));
	xsd_type_definition_w_attributes( xsd, class_name(), attributes );
}


} //core
} //select
} //residue_selector
