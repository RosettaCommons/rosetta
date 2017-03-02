// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/movemap/MoveMapFactory.cc
/// @brief  The MoveMapFactory class builds a MoveMap given a Pose using instructions
///         from ResidueSelectors, and JumpSelectors (and eventually, TorsionIDSelectors)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/select/movemap/MoveMapFactory.hh>

// Package headers
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/jump_selector/JumpSelector.hh>
#include <core/select/jump_selector/util.hh>

// Project headers
#include <core/kinematics/MoveMap.hh>

// Basic headers
#include <basic/datacache/DataMap.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/vector1.hh>


namespace core {
namespace select {
namespace movemap {


MoveMapFactory::MMResAction::MMResAction() : action( mm_enable ) {}
MoveMapFactory::MMResAction::MMResAction( move_map_action act, residue_selector::ResidueSelectorCOP sel ) :
	action( act ),
	selector( sel )
{}
MoveMapFactory::MMResAction::~MMResAction() {}

MoveMapFactory::MMResIndexAction::MMResIndexAction() : index(0), action( mm_enable ) {}
MoveMapFactory::MMResIndexAction::MMResIndexAction( Size ind, move_map_action act, residue_selector::ResidueSelectorCOP sel ) :
	index( ind ),
	action( act ),
	selector( sel )
{}
MoveMapFactory::MMResIndexAction::~MMResIndexAction() {}

MoveMapFactory::MMJumpAction::MMJumpAction() : action( mm_enable ) {}
MoveMapFactory::MMJumpAction::MMJumpAction( move_map_action act, jump_selector::JumpSelectorCOP sel) :
	action( act ),
	selector( sel )
{}
MoveMapFactory::MMJumpAction::~MMJumpAction() {}

MoveMapFactory::MoveMapFactory() :
	use_all_bb_( false ),
	all_bb_setting_( false ),
	use_all_chi_( false ),
	all_chi_setting_( false ),
	use_all_nu_( false ),
	all_nu_setting_( false ),
	use_all_branches_( false ),
	all_branches_setting_( false ),
	use_all_jumps_( false ),
	all_jumps_setting_( false )
{}

MoveMapFactory::MoveMapFactory( MoveMapFactory const & ) = default;

MoveMapFactory & MoveMapFactory::operator = ( MoveMapFactory const & ) = default;

MoveMapFactory::~MoveMapFactory() = default;

void MoveMapFactory::all_bb( bool setting ) { use_all_bb_ = true; all_bb_setting_ = setting; }
void MoveMapFactory::add_bb_action(
	move_map_action action,
	residue_selector::ResidueSelectorCOP selector
) {
	MMResAction res_action{ action, selector };
	bb_actions_.push_back( res_action );
}

void MoveMapFactory::add_bb_index_action(
	Size index,
	move_map_action action,
	residue_selector::ResidueSelectorCOP selector
) {
	MMResIndexAction res_ind_action{ index, action, selector };
	bb_index_actions_.push_back( res_ind_action );
}

void MoveMapFactory::all_chi( bool setting ) { use_all_chi_ = true; all_chi_setting_ = setting; }

void MoveMapFactory::add_chi_action(
	move_map_action action,
	residue_selector::ResidueSelectorCOP selector
) {
  MMResAction res_action{ action, selector };
	chi_actions_.push_back( res_action );
}

void MoveMapFactory::all_nu( bool setting ) { use_all_nu_ = true; all_nu_setting_ = setting; }

void MoveMapFactory::add_nu_action(
	move_map_action action,
	residue_selector::ResidueSelectorCOP selector
) {
	MMResAction res_action{ action, selector };
	nu_actions_.push_back( res_action );
}

void MoveMapFactory::all_branches( bool setting ) { use_all_branches_ = true; all_branches_setting_ = setting; }

void MoveMapFactory::add_branches_action(
	move_map_action action,
	residue_selector::ResidueSelectorCOP selector
) {
	MMResAction res_action{ action, selector };
	branches_actions_.push_back( res_action );
}

void MoveMapFactory::all_jumps( bool setting ) { use_all_jumps_ = true; all_jumps_setting_ = setting; }

void MoveMapFactory::add_jump_action(
	move_map_action action,
	jump_selector::JumpSelectorCOP selector
) {
	MMJumpAction jump_action{ action, selector };
	jump_actions_.push_back( jump_action );
}

/// @brief Construct a MoveMap from the internally held ResidueSelectors and JumpSelectors
kinematics::MoveMapOP
MoveMapFactory::create_movemap_from_pose(
	core::pose::Pose const & pose
) const
{
	kinematics::MoveMapOP mm( new kinematics::MoveMap );
	edit_movemap_given_pose( pose, *mm );
	return mm;
}

void
MoveMapFactory::edit_movemap_given_pose(
	core::pose::Pose const & pose,
	kinematics::MoveMap & mm
) const
{
	using residue_selector::ResidueSelectorCOP;
	using residue_selector::ResidueSubset;
	using jump_selector::JumpSelectorCOP;
	using jump_selector::JumpSubset;

	std::map< ResidueSelectorCOP, ResidueSubset > res_selections;
	std::map< JumpSelectorCOP, JumpSubset > jump_selections;

	///////// Backbone /////////
	if ( use_all_bb_ ) mm.set_bb( all_bb_setting_ );
	for ( auto const & bb_act : bb_actions_ ) {
		ResidueSubset selection;
		if ( res_selections.count( bb_act.selector ) ) {
			selection = res_selections[ bb_act.selector ];
		} else {
			selection = bb_act.selector->apply( pose );
			res_selections[ bb_act.selector ] = selection;
		}

		for ( Size ii = 1; ii <= selection.size(); ++ii ) {
			if ( selection[ ii ] ) {
				mm.set_bb( ii, bb_act.action );
			}
		}
	}

	for ( auto const & bb_indact : bb_index_actions_ ) {
		ResidueSubset selection;
		if ( res_selections.count( bb_indact.selector ) ) {
			selection = res_selections[ bb_indact.selector ];
		} else {
			selection = bb_indact.selector->apply( pose );
			res_selections[ bb_indact.selector ] = selection;
		}

		for ( Size ii = 1; ii <= selection.size(); ++ii ) {
			if ( selection[ ii ] ) {
				mm.set_bb( ii, bb_indact.index, bb_indact.action );
			}
		}
	}

	///////// Chi /////////
	if ( use_all_chi_ ) mm.set_chi( all_chi_setting_ );
	for ( auto const & chi_act : chi_actions_ ) {
		ResidueSubset selection;
		if ( res_selections.count( chi_act.selector ) ) {
			selection = res_selections[ chi_act.selector ];
		} else {
			selection = chi_act.selector->apply( pose );
			res_selections[ chi_act.selector ] = selection;
		}

		for ( Size ii = 1; ii <= selection.size(); ++ii ) {
			if ( selection[ ii ] ) {
				mm.set_chi( ii, chi_act.action );
			}
		}
	}

	///////// Nu /////////
	if ( use_all_nu_ ) mm.set_nu( all_nu_setting_ );
	for ( auto const & nu_act : nu_actions_ ) {
		ResidueSubset selection;
		if ( res_selections.count( nu_act.selector ) ) {
			selection = res_selections[ nu_act.selector ];
		} else {
			selection = nu_act.selector->apply( pose );
			res_selections[ nu_act.selector ] = selection;
		}

		for ( Size ii = 1; ii <= selection.size(); ++ii ) {
			if ( selection[ ii ] ) {
				mm.set_nu( ii, nu_act.action );
			}
		}
	}

	///////// Branches /////////
	if ( use_all_branches_ ) mm.set_branches( all_branches_setting_ );
	for ( auto const & branches_act : branches_actions_ ) {
		ResidueSubset selection;
		if ( res_selections.count( branches_act.selector ) ) {
			selection = res_selections[ branches_act.selector ];
		} else {
			selection = branches_act.selector->apply( pose );
			res_selections[ branches_act.selector ] = selection;
		}

		for ( Size ii = 1; ii <= selection.size(); ++ii ) {
			if ( selection[ ii ] ) {
				mm.set_branches( ii, branches_act.action );
			}
		}
	}


	///////// Jumps /////////
	if ( use_all_jumps_ ) mm.set_jump( all_jumps_setting_ );
	for ( auto const & jump_act : jump_actions_ ) {
		JumpSubset selection;
		if ( jump_selections.count( jump_act.selector ) ) {
			selection = jump_selections[ jump_act.selector ];
		} else {
			selection = jump_act.selector->apply( pose );
			jump_selections[ jump_act.selector ] = selection;
		}

		for ( Size ii = 1; ii <= selection.size(); ++ii ) {
			if ( selection[ ii ] ) {
				mm.set_jump( ii, jump_act.action );
			}
		}
	}

	///////// TO DO: Individual Torsions ////////////
}

/// @brief Initialize any data members of this instance from an input tag
/// and a DataMap object
void MoveMapFactory::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datacache
)
{
	debug_assert( tag->getName() == "MoveMapFactory" );
	if ( tag->hasOption( "bb" ) ) {
		all_bb( tag->getOption< bool >( "bb" ));
	}
	if ( tag->hasOption( "chi" ) ) {
		all_chi( tag->getOption< bool >( "chi" ));
	}
	if ( tag->hasOption( "nu" ) ) {
		all_nu( tag->getOption< bool >( "nu" ));
	}
	if ( tag->hasOption( "branches" ) ) {
		all_branches( tag->getOption< bool >( "branches" ));
	}
	if ( tag->hasOption( "jumps" ) ) {
		all_jumps( tag->getOption< bool >( "jumps" ));
	}

	for ( auto const & subtag : tag->getTags() ) {
		move_map_action enable = move_map_action( subtag->getOption< bool >( "enable", true ) );
		if ( subtag->getName() == "Backbone" ) {
			auto selector = residue_selector::parse_residue_selector( subtag, datacache );
			if ( ! selector ) {
				throw utility::excn::EXCN_Msg_Exception( "Failed to find the (required) residue_selector for a Backbone sub-element of the MoveMapFactory" );
			}
			if ( subtag->hasOption( "bb_tor_index" ) ) {
				Size index = subtag->getOption< Size >( "bb_tor_index" );
				MMResIndexAction action;
				action.index = index;
				action.selector = selector;
				action.action = move_map_action( enable );
				add_bb_index_action( index, enable, selector );
			} else {
				add_bb_action( enable, selector );
			}
		} else if ( subtag->getName() == "Chi" ) {
			auto selector = residue_selector::parse_residue_selector( subtag, datacache );
			if ( ! selector ) {
				throw utility::excn::EXCN_Msg_Exception( "Failed to find the (required) residue_selector for a Chi sub-element of the MoveMapFactory" );
			}
			add_chi_action( enable, selector );
		} else if ( subtag->getName() == "Nu" ) {
			auto selector = residue_selector::parse_residue_selector( subtag, datacache );
			if ( ! selector ) {
				throw utility::excn::EXCN_Msg_Exception( "Failed to find the (required) residue_selector for a Nu sub-element of the MoveMapFactory" );
			}
			add_nu_action( enable, selector );
		} else if ( subtag->getName() == "Branches" ) {
			auto selector = residue_selector::parse_residue_selector( subtag, datacache );
			if ( ! selector ) {
				throw utility::excn::EXCN_Msg_Exception( "Failed to find the (required) residue_selector for a Branches sub-element of the MoveMapFactory" );
			}
			add_branches_action( enable, selector );
		} else if ( subtag->getName() == "Jumps" ) {
			auto selector = jump_selector::parse_jump_selector( subtag, datacache );
			if ( ! selector ) {
				throw utility::excn::EXCN_Msg_Exception( "Failed to find the (required) jump_selector for a Jumps sub-element of the MoveMapFactory" );
			}
			add_jump_action( enable, selector );
		}
		// TO DO: TorsionSelection

	}

}

std::string
MoveMapFactory::element_name()
{
	return "MoveMapFactory";
}

std::string
MoveMapFactory::mmf_ct_namer( std::string const & element )
{
	return "mm_factory_" + element + "_type";
}

void MoveMapFactory::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	XMLSchemaAttribute enable_attr = XMLSchemaAttribute::attribute_w_default(
		"enable", xsct_rosetta_bool, "Enable this DOF? (false for disable)", "true" );

	XMLSchemaSimpleSubelementList subelements;

	AttributeList backbone_attrs;
	backbone_attrs.push_back( enable_attr );
	residue_selector::attributes_for_parse_residue_selector_when_required( backbone_attrs );
	backbone_attrs + XMLSchemaAttribute( "bb_tor_index", xsct_non_negative_integer,
		"If this option is set, then this element will apply to only a particular backbone index; e.g. 3."
		" Indices that are out-of-bounds (e.g.'10' for a residue with 3 backbone dihedrals) have no effect" );

	AttributeList chi_nu_branch_attrs;
	chi_nu_branch_attrs.push_back( enable_attr );
	residue_selector::attributes_for_parse_residue_selector_when_required( chi_nu_branch_attrs );

	AttributeList jump_attrs;
	jump_attrs.push_back( enable_attr );
	jump_selector::attributes_for_parse_jump_selector_when_required( jump_attrs );

	subelements.add_simple_subelement( "Backbone", backbone_attrs, "Tag for enabling or disabling backbone motion for particular residues or, optionally, for enabling particular backbone torsions on particular residues by using a previously-defined ResidueSelector" );
	subelements.add_simple_subelement( "Chi", chi_nu_branch_attrs, "Tag for enabling or disabling sidechain (chi) motion for particular residues by using a previously-defined ResidueSelector" );
	subelements.add_simple_subelement( "Nu", chi_nu_branch_attrs, "Tag for enabling or disabling nu-torsion motion for particular residues by using a previously-defined ResidueSelector" );
	subelements.add_simple_subelement( "Branches", chi_nu_branch_attrs, "Tag for enabling or disabling branch-torsion motion for particular residues by using a previously-defined ResidueSelector" );
	subelements.add_simple_subelement( "Jumps", jump_attrs, "Tag for enabling or disabling particular Jumps using a previously-defined JumpSelector" );

	AttributeList mmfact_attributes;
	mmfact_attributes
		+ XMLSchemaAttribute( "bb", xsct_rosetta_bool, "Enable or disable movement for all backbone torsions (default: disabled)" )
		+ XMLSchemaAttribute( "chi", xsct_rosetta_bool, "Enable or disable movement for all chi torsions (default: disabled)" )
		+ XMLSchemaAttribute( "nu", xsct_rosetta_bool, "Enable or disable movement for all nu torsions (default: disabled)" )
		+ XMLSchemaAttribute( "branches", xsct_rosetta_bool, "Enable or disable movement for all branch torsions (default: disabled)" )
		+ XMLSchemaAttribute( "jump", xsct_rosetta_bool, "Enable or disable movement for all jump DOFs (default: disabled)" );

	XMLSchemaComplexTypeGenerator xsct;
	xsct.element_name( element_name() )
		.complex_type_naming_func( & mmf_ct_namer )
		.description( "A MoveMapFactory can be used to restrict conformational flexibility to a specific set of DOFs taking into account"
			" the conformation of the input Pose. It will construct a MoveMap, (or edit an existing MoveMap), making a series of modifications"
			" to it based on the instructions given. First it will make the highest-level modifications given by the attributes"
			" 'bb', 'chi', 'nu', 'branches', and 'jumps'; it will make these modifications to 'true' (enable) or 'false' (disable) if"
			" they are given in the Tag, but it will not make a modification if it is not given. (It is worth knowing that a default"
			" MoveMap says 'nothing is free to move' but that many Movers will use the MoveMapFactory to edit an existing MoveMap they have"
			" initialized with behaviors such as 'all backbones and sidechain dihedrals are free' so make sure to consult the documentation"
			" for the MoveMapFactory consumer you intend to use.) After it makes the highest-level modifications, it makes intermediate-level"
			" modifications specified by the sub-elements in the order they are provided. These are the Backbone, Chi, Nu, Branches, and Jumps"
			" tags (but the Backbone tag only when the 'bb_tor_index' attribute is absent). For each of these tags, a ResidueSelector is used"
			" to define a group of residues to operate on, and then those residues selected (i.e. whose values are marked 'true') will be operated on."
			" The operation performed is specified by the 'enable' attribute of these sub-elements; if this attribute is not given, then"
			" the operation will be to enable that DOF. Residues that are not selected by the ResidueSelector are not operated on; if your intention"
			" is to disable flexibility for a set of residues, and these residues are already marked as flexible, it is not good enough to leave them"
			" out of a selection that enables a different set of residues. You will have to explicitly select those residues and then perform a disabling"
			" action on them. Finally, the lowest-level operations are performed. For now, this is the particular-backbone-dihedrals enable/disable"
			" actions for Backbone tag (i.e. when the 'bb_tor_index' attribute is provided). These lowest level operations happen last, but within the set"
			" of all lowest-level operations, they are performed in the order that they are provided; this is perhaps confusing." )
		.add_optional_name_attribute()
		.set_subelements_repeatable( subelements )
		.add_attributes( mmfact_attributes )
		.write_complex_type_to_schema( xsd );
}




} //namespace movemap
} //namespace select
} //namespace core

