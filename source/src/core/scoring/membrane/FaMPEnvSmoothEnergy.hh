// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/scoring/membrane/FaMPEnvSmoothEnergy.hh
///
/// @brief  Fullatom Smoothed Membrane Environment Energy
/// @details Updated residue-environment energy (fullatom) by Vladmir in 2010 - smoothed
///    derivatives based on updated statistics. Adapted for mpframework by Rebecca
///    @GrayLab.
///    Last Modified: 7/6/14
///
/// @author  Vladmir Yarov-Yaravoy
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_FaMPEnvSmoothEnergy_hh
#define INCLUDED_core_scoring_membrane_FaMPEnvSmoothEnergy_hh

// Unit Headers
#include <core/scoring/membrane/FaMPEnvSmoothEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextDependentOneBodyEnergy.hh>

// Project headers
#include <core/chemical/AA.hh>
#include <core/conformation/Atom.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <ObjexxFCL/FArray3D.hh>
#include <utility/vector1.hh>

namespace core {
namespace scoring {
namespace membrane {

/// @brief Fullatom Membrane Solvation Energy - Statistically Derived,
/// and smoothed derivatives
class FaMPEnvSmoothEnergy : public core::scoring::methods::ContextDependentOneBodyEnergy  {

public:

	typedef core::scoring::methods::ContextDependentOneBodyEnergy  parent;

public:

	// Constructors ////////////////////

	/// @brief Default Constructor
	FaMPEnvSmoothEnergy();

	/// @brief Create a clone of this energy method
	virtual
	core::scoring::methods::EnergyMethodOP
	clone() const;

	// Scoring Methods ////////////////

	/// @brief Computes dScore/dNumNeighbors for all residues for rapid use in later
	/// atom derivate calculations
	virtual
	void
	setup_for_scoring( core::pose::Pose & pose, core::scoring::ScoreFunction const & ) const;

	/// @brief Causes a neighbor graph update
	virtual
	void
	setup_for_derivatives( core::pose::Pose & pose, core::scoring::ScoreFunction const & sf) const;

	/// @brief Evaluates the one-body energy for a residue
	virtual
	void
	residue_energy(
		core::conformation::Residue const & rsd,
		core::pose::Pose const &,
		core::scoring::EnergyMap &
	) const;

	/// @brief Increments the F1 and F2 derivative vectors for an atom
	virtual
	void
	eval_atom_derivative(
		core::id::AtomID const & atom_id,
		core::pose::Pose const & pose,
		core::kinematics::DomainMap const &,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap const & weights,
		core::Vector & F1,
		core::Vector & F2
	) const;

	/// @brief Unused by the FaMPEnvSmoothEnergy class, returns 0
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
	representative_atom_name( core::chemical::AA const aa ) const;

	/// @brief convert the neighbor counts for a residue and its aa type into
	/// a score and a score derivative.
	void
	calc_energy(
		core::conformation::Residue const & rsd,
		core::pose::Pose const & pose,
		core::Real const neighbor_count,
		core::chemical::AA const aa,
		core::Real & score,
		core::Real & dscore_dneighbor_count
	) const;

	/// @brief given the square distance between a representative atom and a neighbor atom,
	/// return the neighborlyness.  Ramps from 1 down to 0 over a range.
	core::Real
	sigmoidish_neighbor( DistanceSquared const sqdist ) const;

	/// @brief given a pair of atoms, one of which is a neighbor atom and the other of which
	/// is a representative atom, and given the weighted score derivative, increments
	/// the F1 and F2 derivatives.
	void
	increment_f1_f2_for_atom_pair(
		core::conformation::Atom const & atom1,
		core::conformation::Atom const & atom2,
		core::Real weighted_dScore_dN,
		core::Vector & F1,
		core::Vector & F2
	) const;

private:

	mutable utility::vector1< core::Real > residue_N_;
	mutable utility::vector1< core::Real > residue_E_;
	mutable utility::vector1< core::Real > residue_dEdN_;

	ObjexxFCL::FArray3D< core::Real > mem_env_log10_;

	virtual
	core::Size version() const;

};

} // membrane
} // scoring
} // core

#endif // INCLUDED_core_scoring_membrane_FaMPEnvSmoothEnergy_hh
