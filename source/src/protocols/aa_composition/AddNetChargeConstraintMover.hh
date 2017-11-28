// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/aa_composition/AddNetChargeConstraintMover.hh
/// @brief Headers for the AddNetChargeConstraintMover.
/// @details Assigns an NetChargeConstraint to a pose, initializing it from a file.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_aa_composition_AddNetChargeConstraintMover_hh
#define INCLUDED_protocols_aa_composition_AddNetChargeConstraintMover_hh

#include <protocols/aa_composition/AddNetChargeConstraintMover.fwd.hh>
#include <core/scoring/netcharge_energy/NetChargeConstraint.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <core/scoring/constraints/ConstraintSet.fwd.hh>

#include <protocols/filters/Filter.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace aa_composition {

class AddNetChargeConstraintMover : public protocols::moves::Mover {

public:
	typedef core::scoring::constraints::ConstraintSet ConstraintSet;
	typedef core::scoring::constraints::ConstraintSetOP ConstraintSetOP;
	typedef core::scoring::constraints::ConstraintSetCOP ConstraintSetCOP;

public:
	/// @brief Constructor.
	///
	AddNetChargeConstraintMover();

	/// @brief Copy Constructor.
	///
	AddNetChargeConstraintMover( AddNetChargeConstraintMover const &src );

	/// @brief Destructor.
	///
	~AddNetChargeConstraintMover() override;

	/// @brief Copy this object and return a pointer to the copy.
	///
	protocols::moves::MoverOP clone() const override;

	/// @brief Create a new object of this type and return a pointer to it.
	///
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief Actually apply the mover to a pose.
	///
	void apply( Pose & ) override;

	/// @brief Parse RosettaScripts XML tag to set up the mover.
	///
	void
	parse_my_tag( TagCOP tag, basic::datacache::DataMap &data_map, Filters_map const &filters_map, protocols::moves::Movers_map const &movers_map, Pose const &pose ) override;

	/// @brief Get an owning pointer to the NetChargeConstraint object created by this mover.
	/// @details NULL unless the create_constraint_from_file() function has been called.
	inline core::scoring::netcharge_energy::NetChargeConstraintCOP constraint() const { return constraint_; }

	/// @brief Create the NetChargeConstraint object and initialize it from a .charge file.
	///
	void create_constraint_from_file( std::string const &filename );

	/// @brief Create the NetChargeConstraint object from the data from a .charge file.
	/// @details Allows external code to create the constraint object without having it read directly from disk.
	void create_constraint_from_file_contents( std::string const &filecontents );

	/// @brief Add a ResidueSelector to the constraint to use as a mask.
	/// @details The constraint must already have been created with the create_constraint_from_file() function before this function is called.
	void add_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector_in );

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	//Private member variables:

	/// @brief Owning pointer for the NetChargeConstraint object that will be added to the pose.
	/// @details NULL by default; only created when the create_constraint_from_file() function is called.  Must be called before apply time.
	core::scoring::netcharge_energy::NetChargeConstraintOP constraint_;

};

} // aa_composition
} // protocols

#endif
