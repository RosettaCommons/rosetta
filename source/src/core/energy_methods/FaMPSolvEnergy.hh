// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/energy_methods/FaMPSolvEnergy.hh
///
/// @brief  LK-Type Membrane Solvation Energy
///
/// @author  Patrick Barth (Original)
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_FaMPSolvEnergy_hh
#define INCLUDED_core_scoring_membrane_FaMPSolvEnergy_hh

// Unit Headers
#include <core/energy_methods/FaMPSolvEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.hh>

// Project Headers

// Package headers
#include <core/conformation/Atom.fwd.hh>
#include <core/scoring/etable/Etable.fwd.hh>
#include <core/scoring/memb_etable/MembEtable.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

//#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray3D.fwd.hh>

// C++ Headers

namespace core {
namespace energy_methods {

/// @brief Energy Method: Membrane Fullaotm Solvation Energy (LK)
class FaMPSolvEnergy : public core::scoring::methods::ContextDependentTwoBodyEnergy {

public:

	typedef ContextDependentTwoBodyEnergy  parent;

public:

	/// @brief Construct MP Solv energy from standard and membrane etable
	FaMPSolvEnergy(
		core::scoring::etable::EtableCAP etable_in,
		core::scoring::etable::MembEtableCAP memb_etable_in,
		bool const analytic_membetable_evaluation
	);


	/// @brief Clone Energy Method
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/// @brief Setup Energy Method for Derivatives
	void
	setup_for_derivatives(
		pose::Pose & pose,
		core::scoring::ScoreFunction const & scfxn
	) const override;

	/// @brief Evaluate Derivatives
	/// @details Called during graident-based minimization inside dfunc
	/// note: f1 and f2 are not zeroed - contributions are summed
	void
	eval_atom_derivative(
		id::AtomID const & id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		core::scoring::ScoreFunction const & sfxn,
		core::scoring::EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const override;


	/// @brief Compute Residue Pair Energy
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & emap
	) const override;

	/// @brief Define Use of Intraresidue Energies
	bool
	defines_intrares_energy( core::scoring::EnergyMap const & /*weights*/ ) const override { return false; }

	/// @brief Evaluate Intra-Residue Energies
	void
	eval_intrares_energy(
		conformation::Residue const &,
		pose::Pose const &,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & emap
	) const override;

	/// @brief Specify Interaction Cutoff for computing pair energies
	Distance
	atomic_interaction_cutoff() const override;

	/// @brief Provide context graphs
	void
	indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const override;

	/// @brief Setup Energy Method for Scoring
	void
	setup_for_scoring(
		pose::Pose & pose, core::scoring::ScoreFunction const &
	) const override;


	/// @brief Finalize method after computing totals
	void
	finalize_total_energy(
		pose::Pose & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & emap // totals
	) const override;

private: // helper methods

	/// @brief Compute Residue Pair Energies
	void
	get_residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		Real & fa_mbsolv_score
	) const;

	/// @brief Evaluate LK Energy
	Real
	eval_lk(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		Real const & d2,
		Real const & f1,
		Real const & f2,
		bool & debug
	) const;

	/// @brief Compute Change in Energy over distance (for minimization)
	Real
	eval_dE_dR_over_r(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		core::scoring::EnergyMap const & weights,
		Vector & F1,
		Vector & F2,
		Real const & f1,
		Real const & f2
	) const;

	/// @brief Versioning
	core::Size version() const override;

	/// @brief Initialize Energy Method data for derivatives
	void
	init( pose::Pose & pose ) const;

	//solvation component of full atom membrane solvation energy of atom i and j in water and chex
	void
	solvationE(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		Real dis2,
		Real & solvE1,
		Real & solvE2,
		Real & membsolvE1,
		Real & membsolvE2
	) const;

	//solvation function used in solvationE that is independent of membrane depth of atom
	Real
	solv(
		int atom1type,
		int atom2type,
		Real dis
	) const;


	//A portion of the solvation calculated in solv that is only dependent on one atom
	Real
	solv_piece(
		int atom1type,
		Real d
	) const;

	//solvation partial derivative wrt distance from atom i and j
	void
	dsolvationE(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		Real dis2,
		Real & dsolvE1,
		Real & dsolvE2,
		Real & dmembsolvE1,
		Real & dmembsolvE2
	) const;

	Real
	solv_deriv(
		conformation::Atom const & atom,
		Real dis
	) const;

	/// @brief Allocate memory for derivatives
	void setup_for_fullatom( pose::Pose & pose ) const;

private: // data

	// Store Etables at construction
	core::scoring::etable::EtableCAP etable_;
	core::scoring::etable::MembEtableCAP memb_etable_;

	// Store Standard Energies from Etable
	ObjexxFCL::FArray3D< Real > const & solv1_;
	ObjexxFCL::FArray3D< Real > const & solv2_;
	ObjexxFCL::FArray3D< Real > const & dsolv1_;
	ObjexxFCL::FArray3D< Real > const & dsolv2_;
	ObjexxFCL::FArray3D< Real > const & dsolv_;

	// Store Membrane Energies from the etable
	ObjexxFCL::FArray3D< Real > const & memb_solv1_;
	ObjexxFCL::FArray3D< Real > const & memb_solv2_;
	ObjexxFCL::FArray3D< Real > const & memb_dsolv1_;
	ObjexxFCL::FArray3D< Real > const & memb_dsolv2_;

	//Store energies and parameters from etable
	//from Lazaridis 2003 https://doi.org/10.1002/prot.10410
	utility::vector1< Real > const & lk_dgfree_;
	utility::vector1< Real > const & memb_lk_dgfree_;
	utility::vector1< Real > const & lj_radius_;
	utility::vector1< Real > const & lk_volume_;
	utility::vector1< Real > const & lk_lambda_;


	Real const safe_max_dis2_;
	Real const get_bins_per_A2_;

	Real max_dis_;
	Real max_normal_dis_;
	bool const analytic_etable_evaluation_;
	//bool const verbose_;

	// Used only when computing derivatives.
	mutable Real fa_weight_;

	// Arrays used for computing derivatives
	mutable utility::vector1 < utility::vector1 < Real > > fa_proj_;

	mutable utility::vector1 < utility::vector1 < Vector > > fa_f1_;
	mutable utility::vector1 < utility::vector1 < Vector > > fa_f2_;
};

} // scoring
} // core

#endif // INCLUDED_core_energy_methods_FaMPSolvEnergy_hh
