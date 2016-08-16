// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/MembraneEnvSmoothEnergy.hh
/// @brief  Statistically derived smooth membrane environment potential class declaration
/// @brief  and based on EnvSmoothEnergy.hh developed by Andrew Leaver-Fay, Mike Tyka and Phil Bradley
/// @author Vladimir Yarov-Yarovoy

#ifndef INCLUDED_core_scoring_methods_MembraneEnvSmoothEnergy_hh
#define INCLUDED_core_scoring_methods_MembraneEnvSmoothEnergy_hh

// Unit Headers
#include <core/scoring/methods/MembraneEnvSmoothEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextDependentOneBodyEnergy.hh>

// Project headers
#include <core/chemical/AA.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility headers
#include <ObjexxFCL/FArray3D.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {


class MembraneEnvSmoothEnergy : public ContextDependentOneBodyEnergy  {
public:
	typedef ContextDependentOneBodyEnergy  parent;

public:


	MembraneEnvSmoothEnergy();

	/// clone
	virtual
	EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	/// @brief computes dScore/dNumNeighbors for all residues for rapid use in later
	/// atom derivate calculations
	virtual
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	/// @brief causes a neighbor graph update
	virtual
	void
	setup_for_derivatives( pose::Pose & pose, ScoreFunction const & sf) const;

	/// @brief evaluates the one-body energy for a residue
	virtual
	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const &,// pose,
		EnergyMap &// emap
	) const;

	/// @brief increments the F1 and F2 derivative vectors for an atom
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

	/// @brief unused by the MembraneEnvSmoothEnergy class, returns 0
	Distance
	atomic_interaction_cutoff() const;

	/// @brief Tells the scoring function to maintain the TwelveANeighborGraph
	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const;

private:

	/// @brief returns the atom name for the atom used to represent the sidechain for
	/// a particular amino acid; this atom was used to derive the statistics this potential
	/// is based on.
	std::string const &
	representative_atom_name( chemical::AA const aa ) const;

	/// @brief convert the neighbor counts for a residue and its aa type into
	/// a score and a score derivative.
	void
	calc_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		Real const neighbor_count,
		chemical::AA const aa,
		Real & score,
		Real & dscore_dneighbor_count
	) const;

	/// @brief given the square distance between a representative atom and a neighbor atom,
	/// return the neighborlyness.  Ramps from 1 down to 0 over a range.
	Real
	sigmoidish_neighbor( DistanceSquared const sqdist ) const;

	/// @brief given a pair of atoms, one of which is a neighbor atom and the other of which
	/// is a representative atom, and given the weighted score derivative, increments
	/// the F1 and F2 derivatives.
	void
	increment_f1_f2_for_atom_pair(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		Real weighted_dScore_dN,
		Vector & F1,
		Vector & F2
	) const;

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	mutable utility::vector1< Real > residue_N_;
	mutable utility::vector1< Real > residue_E_;
	mutable utility::vector1< Real > residue_dEdN_;

	ObjexxFCL::FArray3D< Real > mem_env_log10_;
	virtual
	core::Size version() const;

};


}
}
}

#endif
