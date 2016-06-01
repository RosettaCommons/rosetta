// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/GlycanResidueSelector.hh
/// @brief  A ResidueSelector for carbohydrates and individual carbohydrate trees. 
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <core/select/residue_selector/GlycanResidueSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>

#include <core/select/residue_selector/util.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh> // for xml schema utility functions
#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>
#include <core/pose/carbohydrates/util.hh>
#include <core/conformation/Residue.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <utility/assert.hh>
#include <utility/string_util.hh>

static THREAD_LOCAL basic::Tracer TR( "core.select.residue_selector.GlycanResidueSelector" );


namespace core {
namespace select {
namespace residue_selector {

using namespace core::select::residue_selector;
using namespace core::pose::carbohydrates;

/// @brief Constructor.
///
GlycanResidueSelector::GlycanResidueSelector():
	ResidueSelector(),
	ref_pose_name_("")
{
	
}

GlycanResidueSelector::GlycanResidueSelector( utility::vector1< bool > root_residues ):
	ResidueSelector(),
	ref_pose_name_("")
{
	set_select_from_branch_residues( root_residues );
}

/// @brief Destructor.
///
GlycanResidueSelector::~GlycanResidueSelector() {}

void
GlycanResidueSelector::set_select_from_branch_residues(utility::vector1<bool> root_residues){
	
	root_residues_ = root_residues;

}


/// @brief Clone function.
/// @details Copy this object and return owning pointer to the copy (created on the heap).
core::select::residue_selector::ResidueSelectorOP
GlycanResidueSelector::clone() const {
	return core::select::residue_selector::ResidueSelectorOP(
		utility::pointer::dynamic_pointer_cast<core::select::residue_selector::ResidueSelector>(
			GlycanResidueSelectorOP( new GlycanResidueSelector(*this) )
		)
	);
}

/// @brief XML parse.
/// @details Parse RosettaScripts tags and set up this mover.
void
GlycanResidueSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & )
{

	if ( tag->hasOption("branch") ) {
		parsed_positions_.push_back( tag->getOption< std::string >("branch"));
	} else if ( tag->hasOption("branches") ) {
		parsed_positions_ = utility::string_split_multi_delim( tag->getOption< std::string >("branches"), ",'`~+*&|;. ");
	}
	
	ref_pose_name_ = tag->getOption< std::string >( "ref_pose_name", ref_pose_name_ );
	
}

std::string GlycanResidueSelector::get_name() const
{
	return GlycanResidueSelector::class_name();
}

std::string GlycanResidueSelector::class_name()
{
	return "GlycanResidueSelector";
}

void GlycanResidueSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute(  "branch",      xs_string )
		+ XMLSchemaAttribute(  "branches",    xs_string )
		+ XMLSchemaAttribute(  "ref_pose_name", xs_string );
	
	xsd_type_definition_w_attributes( xsd, class_name(), attributes );
	
	
}

ResidueSelectorOP
GlycanResidueSelectorCreator::create_residue_selector() const {
	return core::select::residue_selector::ResidueSelectorOP(
		utility::pointer::dynamic_pointer_cast< core::select::residue_selector::ResidueSelector > (
			GlycanResidueSelectorOP( new GlycanResidueSelector )
		)
	);
}

std::string
GlycanResidueSelectorCreator::keyname() const {
	return GlycanResidueSelector::class_name();
}

/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
///
void
GlycanResidueSelectorCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	GlycanResidueSelector::provide_xml_schema( xsd );
}


/// @brief "Apply" function.
/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
/// indicating whether each residue is selected ("true") or not ("false").
ResidueSubset
GlycanResidueSelector::apply(
	core::pose::Pose const & pose
) const {
	
	utility::vector1< bool > root_residues = root_residues_;
	
	utility::vector1< bool > subset (pose.total_residue(), false);
	if (parsed_positions_.size() > 0){
		root_residues.resize(pose.total_residue(), false);
		for ( core::Size i = 1; i <= parsed_positions_.size(); ++ i ) {
			core::Size resnum = core::pose::parse_resnum( parsed_positions_[ i ], pose);
			root_residues[ resnum ] = true;
		}
	}
	
	if (root_residues_.size() > 0){
		for ( core::Size i = 1; i <= root_residues_.size(); ++ i ) {
		
			if (! root_residues[ i ] ) continue;
			
			core::Size resnum = i;
			
			if ( ref_pose_name_ != "" ) {
				resnum = pose.corresponding_residue_in_current( resnum, ref_pose_name_ );
			}
			
			std::pair< utility::vector1< core::Size >, utility::vector1< core::Size > > res_and_tips;

			res_and_tips = get_carbohydrate_residues_upstream( pose, resnum );
			utility::vector1< core::Size > branching_residues = res_and_tips.first;
			//branching_residues.push_back( resnum ); This selects the root and includes it.  Probably not what we want, especially if we select on
			
			for ( core::Size x = 1; x <= branching_residues.size(); ++x ) {
				core::Size branching_resnum = branching_residues[ x ];
				subset[ branching_resnum ] = true;
			
			}
		}
	}
	else{
		for (core::Size i = 1; i <= pose.total_residue(); ++i){
			if (pose.residue( i ).is_carbohydrate()){
				subset[ i ] = true;
			}
		}
	}
	
	return subset;
	
}




} //core
} //select
} //residue_selector
