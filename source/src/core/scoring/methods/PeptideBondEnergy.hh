// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/PeptideBondEnergy.hh
/// @brief  scoring method that defines ideal bond lengths between
/// carbonyl carbon of residue N and nitrogen of residue N+1.
/// @author James Thompson

#ifndef INCLUDED_core_scoring_methods_PeptideBondEnergy_hh
#define INCLUDED_core_scoring_methods_PeptideBondEnergy_hh

// Unit headers
#include <core/scoring/methods/PeptideBondEnergy.fwd.hh>

// Package headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/id/AtomID.fwd.hh>

#include <utility/vector1.hh>



namespace core {
namespace scoring {
namespace methods {

/// @brief PeptideBondEnergy class iterates across all residues in finalize()
/// and determines the penalty between residues i and i+1 by the distance
/// the C-N bond. Evantually I'd also like to add bond angle constraints as
/// well, but that's handled by OmegaTether at the moment.
//class PeptideBondEnergy : public WholeStructureEnergy  {
class PeptideBondEnergy : public ContextIndependentTwoBodyEnergy  {
public:
	typedef ContextIndependentTwoBodyEnergy  parent;

public:

	PeptideBondEnergy();

	/// clone
	virtual
	EnergyMethodOP
	clone() const
	{
		return EnergyMethodOP( new PeptideBondEnergy );
	}

	/// called at the end of energy evaluation
	//virtual
	//void
	//finalize_total_energy(
	//	pose::Pose & pose,
	//	ScoreFunction const &,
	//	EnergyMap & totals
	//) const;

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

	/// called during gradient-based minimization inside dfunc
	/**
		 F1 and F2 are not zeroed -- contributions from this atom are
		 just summed in
	**/
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

	virtual
	void
	indicate_required_context_graphs( utility::vector1< bool > & ) const;

	virtual
	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const {
		return false;
	}

	virtual
	Distance
	atomic_interaction_cutoff() const;
virtual
core::Size version() const;

private:

};

} // methods
} // scoring
} // core


#endif
