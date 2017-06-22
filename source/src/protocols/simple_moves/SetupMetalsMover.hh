// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_moves/SetupMetalsMover.hh
/// @brief  class definition for SetupMetalsMover
/// @author Sharon Guffy (guffy@email.unc.edu)

#ifndef INCLUDED_protocols_simple_moves_SetupMetalsMover_HH
#define INCLUDED_protocols_simple_moves_SetupMetalsMover_HH

// Unit headers
#include <protocols/simple_moves/SetupMetalsMover.fwd.hh>

// Package headers
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

namespace protocols {
namespace simple_moves {

///////////////////////////////////////////////////////////////////////////////
/// @brief A protocols::moves::Mover that
class SetupMetalsMover : public protocols::moves::Mover {

public:
	// default constructor
	/// @brief Constructs a SetupMetalsMover
	SetupMetalsMover();

	/// @brief Copy constructor
	SetupMetalsMover( SetupMetalsMover const & src );

	~SetupMetalsMover() override;


	//General mover methods
	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override { return mover_name(); }
	void show(std::ostream & output=std::cout) const override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief Called by protocols::moves::MoverFactory when constructing new protocols::moves::Movers. Takes care of the specific mover's parsing.
	void parse_my_tag(
		TagCOP,
		basic::datacache::DataMap &,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & ) override;

	//Static methods
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string mover_name();

	//Getters
	//bool
	//get_prevent_setup_metal_bb_variants() const;

	bool
	get_constraints_only() const;

	core::Real
	get_metals_detection_LJ_multiplier() const;

	core::Real
	get_metals_distance_constraint_multiplier() const;

	core::Real
	get_metals_angle_constraint_multiplier() const;

	core::select::residue_selector::ResidueSelectorCOP
	get_metal_selector() const;

	std::string
	get_metal_resnums_string() const;

	bool
	get_remove_hydrogens() const;

	//Setters
	//void
	//set_prevent_setup_metal_bb_variants( bool );

	void
	set_constraints_only( bool );

	void
	set_metals_detection_LJ_multiplier( core::Real );

	void
	set_metals_distance_constraint_multiplier( core::Real );

	void
	set_metals_angle_constraint_multiplier( core::Real );

	void
	set_metal_selector( core::select::residue_selector::ResidueSelectorCOP );

	void
	set_metal_resnums_string( std::string );

	void
	set_remove_hydrogens( bool );

protected:
	///@brief If a residue selector or resnum string is provided, returns any metal ions within the selection.
	///Otherwise returns an empty vector
	utility::vector1< core::Size >
	find_metal_resnums( core::pose::Pose const & );

	void
	set_defaults_from_command_line();
private:
	// data
	//These data members are based on the options currently available on the command line
	//bool prevent_setup_metal_bb_variants_ = false; //default false
	bool remove_hydrogens_ = true;
	bool constraints_only_ = false;
	core::Real metals_detection_LJ_multiplier_ = 1.0; //default 1.0
	core::Real metals_distance_constraint_multiplier_ = 1.0; //default 1.0
	core::Real metals_angle_constraint_multiplier_ = 1.0; //default 1.0

	//TODO: Use residue selector/metal resnums option to only set up specified metal ions
	core::select::residue_selector::ResidueSelectorCOP metal_selector_ = nullptr;
	std::string metal_resnums_string_ = "";
};

std::ostream &operator<< (std::ostream &os, SetupMetalsMover const &mover);

}  // simple_moves
}  // protocols

#endif  // INCLUDED_protocols_simple_moves_SetupMetalsMover_HH
