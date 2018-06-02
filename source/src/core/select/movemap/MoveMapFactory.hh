// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/movemap/MoveMapFactory.hh
/// @brief  The MoveMapFactory class builds a MoveMap given a Pose using instructions
///         from ResidueSelectors, and JumpSelectors (and eventually, TorsionIDSelectors)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_select_movemap_MoveMapFactory_HH
#define INCLUDED_core_select_movemap_MoveMapFactory_HH

// Unit headers
#include <core/select/movemap/MoveMapFactory.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/select/jump_selector/JumpSelector.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>

// Basic headers
#include <basic/datacache/DataMap.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <list>
#include <string>

namespace core {
namespace select {
namespace movemap {

/// @brief Are we turning on or turning off these torsions?
enum move_map_action {
	mm_disable = 0, // To be convertable with "false"
	mm_enable = 1 // To be convertable with "true"
};


class MoveMapFactory : public utility::pointer::ReferenceCount {
public:

	/// @brief Constructor.
	MoveMapFactory();

	MoveMapFactory( MoveMapFactory const & src );
	MoveMapFactory & operator = ( MoveMapFactory const & rhs );

	MoveMapFactoryOP clone() const;

	/// @brief Destructor.
	///
	virtual ~MoveMapFactory();

	void all_bb( bool setting );
	void add_bb_action( move_map_action action, residue_selector::ResidueSelectorCOP selector );
	void add_bb_index_action( Size index, move_map_action action, residue_selector::ResidueSelectorCOP selector );

	void all_chi( bool setting );
	void add_chi_action( move_map_action action, residue_selector::ResidueSelectorCOP selector );

	void all_nu( bool setting );
	void add_nu_action( move_map_action action, residue_selector::ResidueSelectorCOP selector );

	void all_branches( bool setting );
	void add_branches_action( move_map_action action, residue_selector::ResidueSelectorCOP selector );

	void all_jumps( bool setting );
	void add_jump_action( move_map_action action, jump_selector::JumpSelectorCOP selector );

	void all_bondangles( bool setting );
	void add_bondangles_action( move_map_action action, residue_selector::ResidueSelectorCOP selector );

	void all_bondlengths( bool setting );
	void add_bondlengths_action( move_map_action action, residue_selector::ResidueSelectorCOP selector );


	///@brief Option to use specific cartesian settings for the movemap.
	/// Used for glycans as dihedral vs cart MM settings are completely different to move the same
	/// underlying coordinates in IUPAC definitions
	void set_cartesian( bool setting );
	bool cartesian() const;


	/// @brief Construct a MoveMap from the internally held ResidueSelectors and JumpSelectors
	kinematics::MoveMapOP
	create_movemap_from_pose(
		core::pose::Pose const & pose
	) const;

	void
	edit_movemap_given_pose(
		core::pose::Pose const & pose,
		kinematics::MoveMap & mm
	) const;

	/// @brief Initialize any data members of this instance from an input tag
	/// and a DataMap object
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datacache
	);

	static std::string element_name();
	static std::string mmf_ct_namer( std::string const & element );
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	/// @brief An action to be applied to a residue all members of a particular part
	/// (i.e. all backbone torsions, or all chi torsions)
	struct MMResAction {
		MMResAction();
		MMResAction( move_map_action, residue_selector::ResidueSelectorCOP );
		~MMResAction();
		move_map_action action;
		residue_selector::ResidueSelectorCOP selector;
	};

	/// @brief An action for a particular indexed thing (i.e. a particular backbone, or a particular chi torsion).
	struct MMResIndexAction
	{
		MMResIndexAction();
		MMResIndexAction( Size, move_map_action, residue_selector::ResidueSelectorCOP );
		~MMResIndexAction();
		Size index;
		move_map_action action;
		residue_selector::ResidueSelectorCOP selector;
	};

	/// @brief An action for a set of Jumps
	struct MMJumpAction
	{
		MMJumpAction();
		MMJumpAction( move_map_action, jump_selector::JumpSelectorCOP );
		~MMJumpAction();
		move_map_action action;
		jump_selector::JumpSelectorCOP selector;
	};


private:
	bool use_all_bb_ = false; // true if the all_bb_setting_ should be used
	bool all_bb_setting_ = false;
	std::list< MMResAction > bb_actions_;
	std::list< MMResIndexAction > bb_index_actions_;

	bool use_all_chi_ = false; // true if the all_chi_setting_ should be used
	bool all_chi_setting_ = false;
	std::list< MMResAction > chi_actions_;
	//std::list< MMResIndexAction > chi_index_actions_;

	bool use_all_nu_ = false; // true if the all_nu_setting_ should be used
	bool all_nu_setting_ = false;
	std::list< MMResAction > nu_actions_;

	bool use_all_branches_ = false; // true if the all_branches_setting_ should be used
	bool all_branches_setting_ = false;
	std::list< MMResAction > branches_actions_;

	bool use_all_jumps_ = false; // true if the all_jumps_setting_ should be used
	bool all_jumps_setting_ = false;
	std::list< MMJumpAction > jump_actions_;

	bool use_all_bondangles_ = false; // true if the all_bondangles_setting_ should be used
	bool all_bondangles_setting_ = false;
	std::list< MMResAction > bondangles_actions_;

	bool use_all_bondlengths_ = false; // true if the all_bondlengths_setting_ should be used
	bool all_bondlengths_setting_ = false;
	std::list< MMResAction > bondlengths_actions_;

	bool cartesian_specific_ = false; //Option to use specific cartesian settings for the movemap.  Used for glycans as dihedral vs cart MM settings are completely different to move the same underlying coordinates in IUPAC definitions
	// TO DO:
	// Specific torsion selectors!

};


} //namespace movemap
} //namespace select
} //namespace core


#endif
