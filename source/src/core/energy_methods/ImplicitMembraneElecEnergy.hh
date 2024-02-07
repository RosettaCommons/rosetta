// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/energy_methods/ImplicitMembraneElecEnergy.hh
/// @brief Energy method for computing the depth- and membrane-dependent electrostatics energy
/// @author rfalford12 (rfalford12@gmail.com)

#ifndef INCLUDED_core_energy_methods_ImplicitMembraneElecEnergy_hh
#define INCLUDED_core_energy_methods_ImplicitMembraneElecEnergy_hh

// Unit headers
#include <core/energy_methods/ImplicitMembraneElecEnergy.fwd.hh>
#include <core/scoring/methods/TwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethod.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/elec/FA_ElecEnergy.hh>
#include <core/energy_methods/ImplicitMembraneCoulomb.fwd.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.hh>


// Package headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/MinimizationData.fwd.hh>

#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>

#include <core/scoring/DerivVectorPair.fwd.hh>
#include <utility/vector1.hh>

namespace core {
namespace energy_methods {
// namespace elec {

///@brief Energy method for computing the depth- and membrane-dependent electrostatics energy
/// The Two Body Energy Method specifies an interface for all two body methods: both
/// long and short range, both context dependent and independent.  Any two body method
/// must implement this interface as well as the EnergyMethod interface.
class ImplicitMembraneElecEnergy : public core::scoring::elec::FA_ElecEnergy {

public:

	//Simply to keep the interface clean

	typedef core::pose::Pose Pose;
	typedef core::conformation::Residue Residue;
	typedef conformation::RotamerSetBase RotamerSetBase;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::scoring::EnergyMap EnergyMap;
	typedef core::scoring::ResPairMinimizationData ResPairMinimizationData;

	typedef FA_ElecEnergy parent;
	// typedef ContextDependentTwoBodyEnergy grandparent;

public:


	ImplicitMembraneElecEnergy( core::scoring::methods::EnergyMethodOptions const & options );

	// copy constructor for deep copies
	ImplicitMembraneElecEnergy( ImplicitMembraneElecEnergy const & src );

	core::scoring::methods::EnergyMethodOP clone() const override;

	// destructor (important for properly forward-declaring smart-pointer members)
	~ImplicitMembraneElecEnergy() override;

	void
	setup_for_derivatives( pose::Pose & pose, core::scoring::ScoreFunction const & ) const override;

	// void
	// setup_for_packing( pose::Pose & pose, utility::vector1< bool > const &, utility::vector1< bool > const & ) const override;

	void
	setup_for_scoring( pose::Pose & pose, core::scoring::ScoreFunction const & ) const override;

	/// @brief overrides parent class implementation which would have
	/// created several tries
	// void
	// prepare_rotamers_for_packing(
	//  pose::Pose const & pose,
	//  conformation::RotamerSetBase & set ) const override;

	/// @brief overrides parent class implementation which would have
	/// updated a trie
	// void
	// update_residue_for_packing( pose::Pose & pose, Size resid ) const override;

	/// @brief Evaluate the interaction between a given residue pair
	/// accumulating the unweighted energies in an EnergyMap
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::scoring::EnergyMap & emap
	) const override;

	Real
	score_atom_pair(
		core::conformation::Residue const & rsd1,
		core::conformation::Residue const & rsd2,
		core::Size const at1,
		core::Size const at2,
		core::Real const at1_hyd,
		core::Real const at2_hyd,
		core::scoring::EnergyMap & emap,
		core::Real const cpweight,
		core::Real & d2
	) const;


	void
	eval_atom_derivative(
		id::AtomID const & atom_id,
		pose::Pose const & pose,
		kinematics::DomainMap const &,// domain_map,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const override;

	/// @brief Returns "true" because this energy method has not been updated to
	/// use the new derivative evaluation machinery.  Note that this class requires
	/// the definition of this method because it's parent class, FA_ElecEnergy,
	/// HAS been updated to use the new derivative evaluation machinery, and,
	/// if this class did not return "true", it would be asked to evaluate derivatives
	/// in ways it cannot yet evaluate them in.
	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const override { return true; }


	void
	indicate_required_context_graphs( utility::vector1< bool > & /* context_graphs_required */ ) const override;

private:

	core::Size version() const override;
	core::energy_methods::ImplicitMembraneCoulombOP membrane_coulomb_;

};

// } // elec
} // scoring
} // core

#endif //core/energy_methods_ImplicitMembraneElecEnergy_hh
