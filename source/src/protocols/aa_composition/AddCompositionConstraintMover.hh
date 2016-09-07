// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/aa_composition/AddCompositionConstraintMover.hh
/// @brief Headers for the AddCompositionConstraintMover.
/// @details Assigns an AACompositionConstraint to a pose, initializing it from a file.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_aa_composition_AddCompositionConstraintMover_hh
#define INCLUDED_protocols_aa_composition_AddCompositionConstraintMover_hh

#include <protocols/aa_composition/AddCompositionConstraintMover.fwd.hh>
#include <core/scoring/aa_composition_energy/AACompositionConstraint.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <core/scoring/constraints/ConstraintSet.fwd.hh>

#include <protocols/filters/Filter.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace aa_composition {

class AddCompositionConstraintMover : public protocols::moves::Mover {

public:
	typedef core::scoring::constraints::ConstraintSet ConstraintSet;
	typedef core::scoring::constraints::ConstraintSetOP ConstraintSetOP;
	typedef core::scoring::constraints::ConstraintSetCOP ConstraintSetCOP;

public:
	/// @brief Constructor.
	///
	AddCompositionConstraintMover();

	/// @brief Copy Constructor.
	///
	AddCompositionConstraintMover( AddCompositionConstraintMover const &src );

	/// @brief Destructor.
	///
	~AddCompositionConstraintMover() override;

	/// @brief Copy this object and return a pointer to the copy.
	///
	protocols::moves::MoverOP clone() const override;

	/// @brief Create a new object of this type and return a pointer to it.
	///
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief Returns the name of this mover ("AddCompositionConstraintMover").
	///
	std::string get_name() const override;

	/// @brief Actually apply the mover to a pose.
	///
	void apply( Pose & ) override;

	/// @brief Parse RosettaScripts XML tag to set up the mover.
	///
	void
	parse_my_tag( TagCOP tag, basic::datacache::DataMap &data_map, Filters_map const &filters_map, protocols::moves::Movers_map const &movers_map, Pose const &pose ) override;

	/// @brief Get an owning pointer to the AACompositionConstraint object created by this mover.
	/// @details NULL unless the create_constraint_from_file() function has been called.
	inline core::scoring::aa_composition_energy::AACompositionConstraintCOP constraint() const { return constraint_; }

	/// @brief Create the AACompositionConstraint object and initialize it from a .comp file.
	///
	void create_constraint_from_file( std::string const &filename );

	/// @brief Create the AACompositionConstraint object from the data from a .comp file.
	/// @details Allows external code to create the constraint object without having it read directly from disk.
	void create_constraint_from_file_contents( std::string const &filecontents );

	/// @brief Add a ResidueSelector to the constraint to use as a mask.
	/// @details The constraint must already have been created with the create_constraint_from_file() function before this function is called.
	void add_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector_in );

private:
	//Private member variables:

	/// @brief Owning pointer for the AACompositionConstraint object that will be added to the pose.
	/// @details NULL by default; only created when the create_constraint_from_file() function is called.  Must be called before apply time.
	core::scoring::aa_composition_energy::AACompositionConstraintOP constraint_;

};

} // aa_composition
} // protocols

#endif
