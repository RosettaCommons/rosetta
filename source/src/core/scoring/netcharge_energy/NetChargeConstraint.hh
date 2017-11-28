// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/netcharge_energy/NetChargeConstraint.hh
/// @brief Headers for a constraint for constraining sequences to have a desired net charge, analogous to a geometric constraint.
/// @details The corresponding energy term for this constraint is the NetChargeEnergy (netcharge in wts files).
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_core_scoring_netcharge_energy_NetChargeConstraint_hh
#define INCLUDED_core_scoring_netcharge_energy_NetChargeConstraint_hh

#include <core/scoring/aa_composition_energy/SequenceConstraint.hh>
#include <core/scoring/netcharge_energy/NetChargeConstraint.fwd.hh>

#include <core/scoring/constraints/Constraint.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/scoring/netcharge_energy/NetChargeEnergySetup.fwd.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/func/XYZ_Func.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>

#include <core/scoring/ScoreType.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>

#include <string>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace netcharge_energy {

class NetChargeConstraint : public core::scoring::aa_composition_energy::SequenceConstraint {

public: //Constructor, destructor, copy, clone:

	/// @brief Constructor
	NetChargeConstraint();

	/// @brief Copy constructor
	NetChargeConstraint( NetChargeConstraint const & src );

	/// @brief Destructor
	~NetChargeConstraint();

	/// @brief Clone operator
	virtual constraints::ConstraintOP clone() const;

	virtual
	bool operator == ( constraints::Constraint const & /*other*/ ) const;

	virtual
	bool
	same_type_as_me( constraints::Constraint const & other ) const;


public: //Functions that actually do stuff:

	/// @brief Set the selector to be used by this constraint.
	/// @details Clones the input.
	void set_selector( select::residue_selector::ResidueSelectorCOP selector_in );

	/// @brief Const access to the selector.
	/// @details Returns a null pointer if no selector has been specified.
	select::residue_selector::ResidueSelectorCOP selector() const;

	/// @brief Const access to the NetChargeEnergySetup object.
	NetChargeEnergySetupCOP
	netcharge_energy_setup() const;

	/// @brief Initialize the NetChargeEnergySetup object from a file.
	void initialize_from_file( std::string const &filename );

	/// @brief Initialize the NetChargeEnergySetup object from the contents of a file.
	/// @details Allows external code to initialize a constriant object without having the
	/// object read directly from disk.
	void initialize_from_file_contents( std::string const &filecontents );

	/// @brief Print info on the constraint
	void show_def (std::ostream &TO, pose::Pose const &pose) const;

private:
	// Member variables

	/// @brief Owning pointer to a ResidueSelector.
	/// @details Optional; will serve as a mask for this NetChargeConstraint if provided.
	select::residue_selector::ResidueSelectorOP selector_;

	/// @brief NetChargeEnergySetup object stored by this object.
	/// @details Created on object construction
	NetChargeEnergySetupOP netcharge_setup_;

};

} // netcharge_energy
} // scoring
} // core

#endif
