
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/MembraneResidueSelector.hh
/// @brief  Select the memrbane residue in a membrane pose
/// @author Rebecca Alford (rfalford12@gmail.com)

// Unit headers
#include <core/select/residue_selector/MembraneResidueSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>

// Package headers
#include <core/pose/selection.hh>
#include <core/conformation/Residue.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>

// Project Headers
#include <core/conformation/Conformation.hh> 
#include <core/conformation/membrane/MembraneInfo.hh> 

// Utility Headers
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>
#include <utility/string_util.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// C++ headers
#include <utility/assert.hh>
#include <boost/foreach.hpp>

static THREAD_LOCAL basic::Tracer TR( "core.select.residue_selector.MembraneResidueSelector" );

namespace core {
namespace select {
namespace residue_selector {

using namespace core::select::residue_selector;

MembraneResidueSelector::MembraneResidueSelector() {}

MembraneResidueSelector::~MembraneResidueSelector() {}

ResidueSubset
MembraneResidueSelector::apply( core::pose::Pose const & pose ) const
{

	using namespace core::conformation::membrane;
	
	// Check that the pose is a membrane protein
	if ( !pose.conformation().is_membrane()) {
		utility_exit_with_message( "Cannot select membrane residue in a non-membrane pose" );
	}

	// Add the residues contained in each selected TM span to the residue subset
	ResidueSubset membrane_selection( pose.total_residue(), false); 
	core::Size rsdnum( pose.conformation().membrane_info()->membrane_rsd_num() ); 
	membrane_selection[ rsdnum ] = true; 
	return membrane_selection;
}

void
MembraneResidueSelector::parse_my_tag(
	utility::tag::TagCOP,
	basic::datacache::DataMap & ) 
{}

std::string MembraneResidueSelector::get_name() const {
	return MembraneResidueSelector::class_name();
}

ResidueSelectorOP
MembraneResidueSelector::clone() const {
	return ResidueSelectorOP( new MembraneResidueSelector(*this) );
}

std::string MembraneResidueSelector::class_name() {
	return "MembraneResidueSelector";
}

void
MembraneResidueSelector::provide_selector_xsd( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attributes;
	xsd_type_definition_w_attributes( xsd, class_name(), attributes );
}
	
ResidueSelectorOP
MembraneResidueSelectorCreator::create_residue_selector() const {
	return ResidueSelectorOP( new MembraneResidueSelector );
}

std::string
MembraneResidueSelectorCreator::keyname() const {
	return MembraneResidueSelector::class_name();
}

void
MembraneResidueSelectorCreator::provide_selector_xsd( utility::tag::XMLSchemaDefinition & xsd ) const {
	MembraneResidueSelector::provide_selector_xsd( xsd );
}

} // residue_selector
} // select
} // core
