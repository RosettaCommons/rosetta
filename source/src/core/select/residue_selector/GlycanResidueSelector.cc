// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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

static basic::Tracer TR( "core.select.residue_selector.GlycanResidueSelector" );


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
	include_root_ = false;
	root_residue_ = 0;
}

GlycanResidueSelector::GlycanResidueSelector( utility::vector1< bool > root_residues, bool include_root /*false*/ ):
	ResidueSelector(),
	ref_pose_name_(""),
	include_root_(include_root)
{
	root_residue_ = 0;
	set_select_from_branch_residues( root_residues );
}

GlycanResidueSelector::GlycanResidueSelector( core::Size root_residue, bool include_root /*false*/ ):
	ResidueSelector(),
	ref_pose_name_(""),
	include_root_(include_root)
{
	set_select_from_branch_residue( root_residue );
}

/// @brief Destructor.
///
GlycanResidueSelector::~GlycanResidueSelector() = default;

void
GlycanResidueSelector::set_select_from_branch_residues(utility::vector1<bool> root_residues){
	root_residue_ = 0;
	root_residues_ = root_residues;

}

void
GlycanResidueSelector::set_select_from_branch_residue(core::Size root_residue){
	root_residues_.clear();
	root_residue_ = root_residue;

}

void
GlycanResidueSelector::set_include_root(bool include_root){
	include_root_ = include_root;
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

ResidueSelectorOP
GlycanResidueSelectorCreator::create_residue_selector() const {
	return core::select::residue_selector::ResidueSelectorOP(
		utility::pointer::dynamic_pointer_cast< core::select::residue_selector::ResidueSelector > (
		GlycanResidueSelectorOP( new GlycanResidueSelector )
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
	include_root_ = tag->getOption< bool >("include_root", include_root_);
}

std::string GlycanResidueSelector::get_name() const
{
	return GlycanResidueSelector::class_name();
}

std::string GlycanResidueSelector::class_name()
{
	return "Glycan";
}

void GlycanResidueSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attributes;
	attributes

		+ XMLSchemaAttribute(  "branch",      xs_string ,
		"Set the residue to select from.  "
		"These can be the branch points of the glycans or carbohydrate residues from which to select the downstream branch from. "
		"Like the rest of a tree from a particular position.  "
		"That position could be the trunk or individual branches, which keep branching out. "
		"Pose (23) or PDB numbering (24A)" )

		+ XMLSchemaAttribute(  "branches",    xs_string ,
		"Set the residues to select from. Separated by spaces or commas. "
		"These can be the branch points of the glycans or carbohydrate residues from which to select the downstream branch from. "
		"Like the rest of a tree from a particular position.  "
		"That position could be the trunk or individual branches, which keep branching out. Pose (23) or PDB numbering (24A)" )

		+ XMLSchemaAttribute(  "ref_pose_name", xs_string ,
		"If using a Reference Pose, set the name" )

		+ XMLSchemaAttribute::attribute_w_default(  "include_root", xsct_rosetta_bool ,
		"Should we include the root (branch/branches) residues in our selection or just go from that?", "false");

	xsd_type_definition_w_attributes( xsd, class_name(),
		"A ResidueSelector for carbohydrates and individual carbohydrate trees. "
		"Selects all Glycan residues if no option is given or the branch going out from the root residue. "
		"Selecting from root residues allows you to choose the whole glycan branch or only tips, etc.", attributes );



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
	utility::vector1< bool > subset (pose.size(), false);

	if ( root_residue_ != 0 ) {
		root_residues.clear();
		root_residues.resize(pose.total_residue(), false);
		root_residues[ root_residue_ ] = true;

	}

	if ( ! parsed_positions_.empty() ) {
		root_residues.resize(pose.size(), false);
		for ( core::Size i = 1; i <= parsed_positions_.size(); ++ i ) {
			core::Size resnum = core::pose::parse_resnum( parsed_positions_[ i ], pose);
			root_residues[ resnum ] = true;
		}
	}
	if ( ! root_residues.empty() ) {
		for ( core::Size i = 1; i <= root_residues.size(); ++i ) {
			if ( ! root_residues[ i ] ) continue;

			core::Size resnum = i ;

			if ( ref_pose_name_ != "" ) {
				resnum = pose.corresponding_residue_in_current( resnum, ref_pose_name_ );
			}

			utility::vector1< core::Size > branching_residues = get_carbohydrate_residues_of_branch( pose, resnum );

			if ( include_root_ ) {
				branching_residues.push_back( resnum );
			}

			for ( core::Size x = 1; x <= branching_residues.size(); ++x ) {
				core::Size branching_resnum = branching_residues[ x ];
				subset[ branching_resnum ] = true;

			}
		}
	} else {
		for ( core::Size i = 1; i <= pose.size(); ++i ) {
			if ( pose.residue( i ).is_carbohydrate() ) {
				subset[ i ] = true;
			}
		}
	}
	return subset;

}




} //core
} //select
} //residue_selector
