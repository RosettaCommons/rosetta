// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/rna/RNP_LowResEnergy.hh
/// @brief  Statistically derived rotamer pair potential class declaration
/// @author Rhiju Das


#ifndef INCLUDED_core_scoring_rna_RNP_LowResEnergy_hh
#define INCLUDED_core_scoring_rna_RNP_LowResEnergy_hh

// Unit Headers
#include <core/scoring/rna/RNP_LowResEnergy.fwd.hh>
#include <core/scoring/rna/RNP_LowResPotential.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/rna/RNA_LowResolutionPotential.fwd.hh>
#include <core/scoring/rna/RNA_RawBaseBaseInfo.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <core/kinematics/Stub.fwd.hh>
#include <utility/vector1.hh>


// Utility headers


namespace core {
namespace scoring {
namespace rna {


class RNP_LowResEnergy : public methods::ContextIndependentTwoBodyEnergy  {
public:
	typedef methods::ContextIndependentTwoBodyEnergy  parent;
public:


	RNP_LowResEnergy();


	/// clone
	virtual
	methods::EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	virtual
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	//virtual
	//void
	//setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const;

	//virtual
	//void
	//setup_for_packing( pose::Pose & pose, utility::vector1< bool > const &, utility::vector1< bool > const & designing_residues ) const {};

	virtual
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const;

	virtual
	void
	eval_intrares_energy(
		conformation::Residue const &,
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap &
	) const {}

//	virtual
//	void
//	eval_atom_derivative(
//		id::AtomID const & atom_id,
//		pose::Pose const & pose,
//		kinematics::DomainMap const & domain_map,
//		ScoreFunction const & scorefxn,
//		EnergyMap const & weights,
//		Vector & F1,
//		Vector & F2
//	) const {};

	virtual
	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const { return false; }

	//virtual
	//void
	//finalize_total_energy(
	//	pose::Pose & pose,
	//	ScoreFunction const &,
	//	EnergyMap &// totals
	//) const;

	virtual
	Distance
	atomic_interaction_cutoff() const;

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const {}


	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	// const-ref to scoring database
	rna::RNP_LowResPotential const & potential_;

	virtual
	core::Size version() const;

};


} //rna
} //scoring
} //core

#endif // INCLUDED_core_scoring_ScoreFunction_HH
