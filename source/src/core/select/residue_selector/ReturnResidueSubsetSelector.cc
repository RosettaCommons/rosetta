// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/residue_selector/ReturnResidueSubsetSelector.hh
/// @brief  A simple selector that returns the set subset.  This to enable simplification of code-based interfaces to residue selectors,
/// so that one may accept only selectors, but using this selector, we can set subsets.  This greatly reduces the interface complexity
/// and code-complexity arising from accepting BOTH ResidueSubsets and ResidueSelectors (Which I'm terribly sick of doing at this point).
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <core/select/residue_selector/ReturnResidueSubsetSelector.hh>
//#include <core/select/residue_selector/ResidueSelectorCreators.hh>

#include <core/select/residue_selector/util.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh> // for xml schema utility functions
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <utility/assert.hh>

static THREAD_LOCAL basic::Tracer TR( "core.select.residue_selector.ReturnResidueSubsetSelector" );

namespace core {
namespace select {
namespace residue_selector {

/// @brief Constructor.
///
ReturnResidueSubsetSelector::ReturnResidueSubsetSelector():
	ResidueSelector()
{
}

ReturnResidueSubsetSelector::ReturnResidueSubsetSelector( ResidueSubset const & subset):
	ResidueSelector(),
	subset_(subset)
{

}

void
ReturnResidueSubsetSelector::set_residue_subset( ResidueSubset const & subset ){
	subset_ = subset;
}


/// @brief Destructor.
///
ReturnResidueSubsetSelector::~ReturnResidueSubsetSelector() {}

 ///@brief Copy Constructor.  Usually not necessary unless you need deep copying (e.g. OPs)
 ReturnResidueSubsetSelector::ReturnResidueSubsetSelector(ReturnResidueSubsetSelector const & src):
 ResidueSelector( src ),
 subset_( src.subset_ )
{

}

/// @brief Clone function.
/// @details Copy this object and return owning pointer to the copy (created on the heap).
ReturnResidueSubsetSelector::ResidueSelectorOP
ReturnResidueSubsetSelector::clone() const {
	return ResidueSelectorOP( new ReturnResidueSubsetSelector(*this) );
}

std::string ReturnResidueSubsetSelector::get_name() const
{
	return ReturnResidueSubsetSelector::class_name();
}

std::string ReturnResidueSubsetSelector::class_name()
{
	return "ReturnResidueSubsetSelector";
}


/// @brief XML parse.
/// @details Parse RosettaScripts tags and set up this mover.
void
ReturnResidueSubsetSelector::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap & )
{
}

void ReturnResidueSubsetSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & )
{
	//Syntax Example:
	//using namespace utility::tag;
	//AttributeList attributes;
	//attributes
	// + XMLSchemaAttribute::attribute_w_default(  "select_positive_phi",      xsct_rosetta_bool, "true" )
	// + XMLSchemaAttribute::attribute_w_default(  "ignore_unconnected_upper", xsct_rosetta_bool, "true" );
	//xsd_type_definition_w_attributes( xsd, class_name(), attributes );

}


/// @brief "Apply" function.
/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
/// indicating whether each residue is selected ("true") or not ("false").
ReturnResidueSubsetSelector::ResidueSubset
ReturnResidueSubsetSelector::apply(
	core::pose::Pose const & pose//pose
) const {

	if ( subset_.empty() ) {
		utility_exit_with_message("A subset must be set for the ReturnResidueSubsetSelector!");
	}

	debug_assert(subset_.size() == pose.size()); //Assert that we are returning a subset FOR THIS POSE //
	return subset_;
}


} //core
} //select
} //residue_selector



