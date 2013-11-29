// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/CustomAtomPairEnergy.hh
/// @brief
/// @author James Thompson


#ifndef INCLUDED_core_scoring_methods_CustomAtomPairEnergy_hh
#define INCLUDED_core_scoring_methods_CustomAtomPairEnergy_hh

// Unit Headers
#include <core/scoring/methods/CustomAtomPairEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/func/SOGFunc_Impl.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>

//#include <core/pack/task/PackerTask.fwd.hh>

// Utility headers

namespace core {
namespace scoring {
namespace methods {

///
class CustomAtomPairEnergy : public ContextIndependentTwoBodyEnergy  {
public:
	typedef ContextIndependentTwoBodyEnergy parent;

public:

	///
	CustomAtomPairEnergy( Size const cst_seq_sep );

	/// clone
	virtual
	EnergyMethodOP
	clone() const;

	///
	virtual
	void
	setup_for_packing( pose::Pose & pose, utility::vector1< bool > const &, utility::vector1< bool > const & ) const;

	///
	virtual
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	///
	virtual
	void
	setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const;

	virtual
	void
	prepare_rotamers_for_packing(
		pose::Pose const & pose,
		conformation::RotamerSetBase & set
	) const;

	virtual
	void
	update_residue_for_packing(
		pose::Pose &,
		Size resid ) const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

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
	eval_atom_derivative(
		id::AtomID const & atom_id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		ScoreFunction const &,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const;

	virtual
	Distance
	atomic_interaction_cutoff() const;

	/// @details non-virtual accessor for speed
	Distance
	interaction_cutoff() const;

	virtual
	void indicate_required_context_graphs(
		utility::vector1< bool > & context_graphs_required
	) const;

	virtual
	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const;

	virtual
	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	// constraints with more than this sequence separation are not active
	Size const max_cst_seq_sep_;
	mutable std::string fn_;

	//mutable bool init_;
	mutable utility::vector1< utility::vector1< core::scoring::constraints::SOGFunc_Impl > > funcs_;
	mutable utility::vector1< utility::vector1< bool > > have_cst_;
virtual
core::Size version() const;
};

} // constraints
} // scoring
} // core

#endif // INCLUDED_core_scoring_ScoreFunction_HH
