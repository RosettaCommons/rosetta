// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		core/scoring/membrane/FaMPEnvEnergy.hh
///
/// @brief		LK-Type Membrane Environment Energy
/// @details	Last Modified: 5/13/14
///
/// @author		Patrick Barth (Original)
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_FaMPEnvEnergy_hh
#define INCLUDED_core_scoring_membrane_FaMPEnvEnergy_hh

// Unit headers
#include <core/scoring/membrane/FaMPEnvEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.hh>

// Project Headers
#include <core/conformation/membrane/SpanningTopology.fwd.hh>

// Package headers
#include <core/scoring/memb_etable/MembEtable.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Conformation.fwd.hh> 

#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <ObjexxFCL/FArray1.fwd.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

// C++ Headers
#include <cstdlib>

namespace core {
namespace scoring {
namespace membrane {

/// @brief Fullatom Membrane Environment Energy
class FaMPEnvEnergy : public methods::ContextDependentOneBodyEnergy {
	
public:

	typedef ContextDependentOneBodyEnergy parent;
	
	/// @brief Construct Energy Method from Etable
	FaMPEnvEnergy( etable::MembEtableCAP memb_etable_in );
	
	/// @brief Clone Energy Method
	virtual
	methods::EnergyMethodOP
	clone() const;
	
	/// @brief Compute Per-Residue Energies
	virtual
	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		EnergyMap & emap
		) const;
	
	/// @brief Fianlzie Total Per-Residue Energies
	virtual
	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & emap
		) const;

	/// @brief Setup for Computing Derivatives
	virtual
	void
	setup_for_derivatives(
		pose::Pose & pose,
		ScoreFunction const & scfxn
		) const;

	/// @brief Evaluate Per-Atom Derivatives
	virtual
	void
	eval_atom_derivative(
		id::AtomID const & id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		ScoreFunction const & sfxn,
		EnergyMap const & emap,
		Vector & F1,
		Vector & F2
		) const;
	
	/// @brief Fa_MbenvEnergy is context independent
	virtual
	void
	indicate_required_context_graphs( utility::vector1< bool > & ) const;
	
	/// @brief Setup Method for initial scoring
	void
	setup_for_scoring(
		pose::Pose & pose,
		ScoreFunction const &
		) const;
	
private: // helper methods
	
	/// @brief Evaluate Per-Atom Env term
	Real
	eval_fa_mbenv(
				  conformation::Atom const & atom1,
				  Real const & f1
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
	
	/// @brief Helper Method - Compute Fa Derivatives
	core::Real
	compute_fa_deriv(
		core::Real z_position,
		core::Real thickness,
		core::Real steepness
		) const;
	
	/// @brief Helper Method - Compute Fa Proj coordinate
	core::Vector
	compute_fa_proj_coord(
		core::Real z_position,
		core::Vector xyz,
		core::Vector center,
		core::Vector normal
		) const;
		
	/// @brief Allocate memory for derivatives
	void setup_for_fullatom( pose::Pose & pose ) const;
	
private:
	
	// Store a copy of the etable in construction
	etable::MembEtableCAP memb_etable_;
	
	// Make copies from the etable
	ObjexxFCL::FArray1< Real > const & lk_dgrefce_;
	ObjexxFCL::FArray1< Real > const & memb_lk_dgrefce_;
	
	// Store mbenv weight when computing derivatives
	mutable Real fa_mbenv_weight_;
	
	// Arrays used for computing derivatives
	mutable utility::vector1 < utility::vector1 < Real > > fa_proj_;
	mutable utility::vector1 < utility::vector1 < Real > > fa_z_position_;
	mutable utility::vector1 < utility::vector1 < Vector > > fa_proj_coord_;
	mutable utility::vector1 < utility::vector1 < Real > > fa_proj_deriv_;
	
};
	
} // membrane
} // scoring
} // core

#endif // INCLUDED_core_scoring_methods_FaMPEnvEnergy_hh
