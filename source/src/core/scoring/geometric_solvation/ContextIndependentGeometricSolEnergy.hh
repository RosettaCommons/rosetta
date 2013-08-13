// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/ContextIndependentGeometricSolEnergy.hh
/// @brief  Geometric solvation energy.
/// @author Parin Sripakdeevong (sripakpa@stanford.edu), Rhiju Das (rhiju@stanford.edu)



#ifndef INCLUDED_core_scoring_geometric_solvation_ContextIndependentGeometricSolEnergy_hh
#define INCLUDED_core_scoring_geometric_solvation_ContextIndependentGeometricSolEnergy_hh

// Unit Headers
#include <core/scoring/geometric_solvation/ContextIndependentGeometricSolEnergy.fwd.hh>
#include <core/scoring/geometric_solvation/GeometricSolEnergyEvaluator.fwd.hh>
#include <core/types.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

//Auto Headers
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>


namespace core {
namespace scoring {
namespace geometric_solvation {

///
class ContextIndependentGeometricSolEnergy : public methods::ContextIndependentTwoBodyEnergy  {
public:
	typedef methods::ContextIndependentTwoBodyEnergy  parent;
public:

	///
	ContextIndependentGeometricSolEnergy( methods::EnergyMethodOptions const & options );

	///@brief copy c-tor
	ContextIndependentGeometricSolEnergy( ContextIndependentGeometricSolEnergy const & src );

	/// clone
	virtual
	methods::EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	/// This evaluates everything for now,
	/// but eventually may want to split this
	/// based on backbone/backbone vs. others,
	/// as is carried out in HBondEnergy.cc
	virtual
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const;

	// Undefined, commenting out to fix PyRosetta build
	/* void
	eval_atom_derivative_intra_RNA(
		 id::AtomID const & atom_id,
		 pose::Pose const & pose,
		 EnergyMap const & weights,
		 Vector & F1,
		 Vector & F2
	) const;
	*/

	/// f1 and f2 are zeroed
	virtual
	void
	eval_atom_derivative(
		id::AtomID const & atom_id,
		pose::Pose const & pose,
		kinematics::DomainMap const &,
		ScoreFunction const &,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const;


	virtual
	bool
	defines_intrares_energy( EnergyMap const & weights ) const;

	virtual
	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & ,
		EnergyMap & emap
	) const;

	virtual
	Distance
	atomic_interaction_cutoff() const;

	virtual
	core::Size version() const;

	///@brief GeometricSolEnergy is context sensitive
	virtual
	void indicate_required_context_graphs(
		utility::vector1< bool > & context_graphs_required ) const;

private:

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////
	methods::EnergyMethodOptionsOP options_;

	GeometricSolEnergyEvaluatorOP evaluator_;


};

} // hbonds
} // scoring
} // core

#endif
