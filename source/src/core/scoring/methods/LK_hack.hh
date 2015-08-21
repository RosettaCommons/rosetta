// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/LK_hack.hh
/// @brief  LK Solvation using hemisphere culling class declaration
/// @author David Baker
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_core_scoring_methods_LK_hack_hh
#define INCLUDED_core_scoring_methods_LK_hack_hh

// Unit Headers
#include <core/scoring/methods/LK_hack.fwd.hh>

// Package headers
#include <core/conformation/Atom.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/etable/Etable.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

// Utility headers

#include <utility/vector1.hh>
#include <ObjexxFCL/FArray3D.fwd.hh>


namespace core {
namespace scoring {
namespace methods {


class LK_hack : public ContextIndependentTwoBodyEnergy  {
public:
	typedef ContextIndependentTwoBodyEnergy  parent;
public:

	LK_hack( etable::Etable const & etable_in );


	/// clone
	virtual
	EnergyMethodOP
	clone() const;

	LK_hack( LK_hack const & src );

	virtual
	void
	setup_for_derivatives(
		pose::Pose & pose,
		ScoreFunction const & scfxn
	) const;

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
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const;


	virtual
	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const { return false; }

	virtual
	void
	eval_intrares_energy(
		conformation::Residue const &,
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap &
	) const {}

	virtual
	Distance
	atomic_interaction_cutoff() const;


	void indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const;

private:

	void
	allocate_appropriate_memory( pose::Pose const & pose ) const;

	void
	calculate_orientation_vectors_and_pseudo_base_atoms( pose::Pose const & pose ) const;

	void
	calculate_derivatives_for_atoms_and_pseudo_base_atoms( pose::Pose const & pose ) const;

	void
	calculate_derivatives_for_residue_pair
	(
		pose::Pose const & pose,
		Size lower_res_id,
		Size upper_res_id
	) const;

	Real
	eval_dE_dR_over_r(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		Distance const one_over_d,
		DistanceSquared const d2,
		Vector & f1_1,
		Vector & f2_1,
		Vector & f1_2,
		Vector & f2_2
	) const;

	void
	distribute_pseudo_base_atom_derivatives( pose::Pose const & pose ) const;

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:
	etable::Etable const & etable_; // shouldn't this be a pointer? Reference count information is (dangerously) lost when
	//a reference is taken, instead of a smart pointer.  There's the potential for a dangling reference with this.


	/// these guys are taken from the etable
	ObjexxFCL::FArray3D< Real > const & solv1_;
	ObjexxFCL::FArray3D< Real > const & solv2_;

	ObjexxFCL::FArray3D< Real > const & dsolv1_;

	Real const safe_max_dis2_;
	Real const get_bins_per_A2_;

	/// Used soley when calculating derivatives
	/// Could/should be moved into the Pose's cachable data.
	mutable utility::vector1< utility::vector1< Size > >   nneighbs_;
	mutable utility::vector1< utility::vector1< Vector > > orientation_vectors_;
	mutable utility::vector1< utility::vector1< Vector > > base_pseudo_atom_centers_;
	mutable utility::vector1< utility::vector1< std::pair< Vector, Vector > > > atom_f1_f2s_;
	mutable utility::vector1< utility::vector1< std::pair< Vector, Vector > > > base_pseudo_atom_f1_f2s_;
	mutable Real lk_hack_weight_; // hold this while calculating derivatives.
	virtual
	core::Size version() const;
};

}
}
}

#endif // INCLUDED_core_scoring_methods_LK_hack_HH
