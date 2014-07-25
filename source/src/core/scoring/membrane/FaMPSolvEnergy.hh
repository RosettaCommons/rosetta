// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		core/scoring/membrane/FaMPSolvEnergy.hh
///
/// @brief		LK-Type Membrane Solvation Energy
/// @details	Last Modified: 5/13/14
///
/// @author		Patrick Barth (Original)
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_FaMPSolvEnergy_hh
#define INCLUDED_core_scoring_membrane_FaMPSolvEnergy_hh

// Unit Headers
#include <core/scoring/membrane/FaMPSolvEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.hh>

// Project Headers
#include <core/conformation/Conformation.hh> 
#include <core/conformation/membrane/SpanningTopology.fwd.hh>

// Package headers
#include <core/conformation/Atom.fwd.hh>
#include <core/scoring/etable/Etable.fwd.hh>
#include <core/scoring/memb_etable/MembEtable.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

#include <ObjexxFCL/FArray3D.fwd.hh>

// C++ Headers
#include <cstdlib>

namespace core {
namespace scoring {
namespace membrane {

/// @brief Energy Method: Membrane Fullaotm Solvation Energy (LK)
class FaMPSolvEnergy : public methods::ContextDependentTwoBodyEnergy {

public:

	typedef ContextDependentTwoBodyEnergy  parent;

public:
	
	/// @brief Construct MP Solv energy from standard and membrane etable
	FaMPSolvEnergy(
		etable::EtableCAP etable_in,
		etable::MembEtableCAP memb_etable_in
		);
	
	
	/// @brief Clone Energy Method
	virtual
	methods::EnergyMethodOP
	clone() const;
	
	/// @brief Setup Energy Method for Derivatives
	virtual
	void
	setup_for_derivatives(
	  pose::Pose & pose,
	  ScoreFunction const & scfxn
	  ) const;
	 
	/// @brief Evaluate Derivatives
	/// @details Called during graident-based minimization inside dfunc
	///	note: f1 and f2 are not zeroed - contributions are summed
	virtual
	void
	eval_atom_derivative(
	 id::AtomID const & id,
	 pose::Pose const & pose,
	 kinematics::DomainMap const & domain_map,
	 ScoreFunction const & sfxn,
	 EnergyMap const & weights,
	 Vector & F1,
	 Vector & F2
	 ) const;
	
	
	/// @brief Compute Residue Pair Energy
	virtual
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
		) const;
	
	/// @brief Define Use of Intraresidue Energies
	virtual
	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const { return false; }
	
	/// @brief Evaluate Intra-Residue Energies
	virtual
	void
	eval_intrares_energy(
						 conformation::Residue const &,
						 pose::Pose const &,
						 ScoreFunction const &,
						 EnergyMap & emap
						 ) const;
	
	/// @brief Specify Interaction Cutoff for computing pair energies
	virtual
	Distance
	atomic_interaction_cutoff() const;
	
	/// @brief Provide context graphs
	void
	indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const;
	
	/// @brief Setup Energy Method for Scoring
	void
	setup_for_scoring(
		pose::Pose & pose, ScoreFunction const &
		) const;
	
	
	/// @brief Finalize method after computing totals
	virtual
	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & emap // totals
		) const;
	
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
		Real & deriv,
		Real const & f1,
		Real const & f2,
		bool & debug
		) const;
	
	/// @brief Compute Change in Energy over distance (for minimization)
	Real
	eval_dE_dR_over_r(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2,
		Real const & f1,
		Real const & f2
		) const;
	
	/// @brief Versioning
	virtual
	core::Size version() const;
	
	/// @brief Initialize Energy Method data for derivatives
	void
	init( pose::Pose & pose ) const;
	
	/// @brief Helper Method - Compute Fa Proj
	core::Real
	compute_fa_proj(
		core::Real z_position,
		core::Real thickness,
		core::Real steepness
		) const;

		
	/// @brief Allocate memory for derivatives
	void setup_for_fullatom( pose::Pose & pose ) const;
	
private: // data

	// Store Etables at construction
	etable::EtableCAP etable_;
	etable::MembEtableCAP memb_etable_;
	
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
	
	Real const safe_max_dis2_;
	Real const get_bins_per_A2_;
	
	bool const verbose_;
	
	// Used only when cmputing derivatives.
	mutable Real fa_weight_;
	
	// Arrays used for computing derivatives
	mutable utility::vector1 < utility::vector1 < Real > > fa_proj_;
	mutable utility::vector1 < utility::vector1 < Real > > fa_z_position_;

};

} // membrane
} // scoring
} // core

#endif // INCLUDED_core_scoring_methods_FaMPSolvEnergy_hh
