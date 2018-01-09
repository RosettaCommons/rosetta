// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/scoring/membrane/MPHelicalityEnergy.hh
///
/// @brief  Fullatom and centroid level smooth membrane non-helicality penalty
/// @details
/// @FlesihmanLab
/// @Last Modified: 20/2/17
///
/// @author Jonathan Weinstein (jonathan.weinstein@weizmann.ac.il)
/// @author Assaf Elazar
/// @author Sarel Fleishman

#ifndef INCLUDED_core_scoring_membrane_MPHelicalityEnergy_hh
#define INCLUDED_core_scoring_membrane_MPHelicalityEnergy_hh

// Unit Headers
#include <core/scoring/membrane/MPHelicalityEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextDependentOneBodyEnergy.hh>

// Project headers
#include <core/chemical/AA.hh>
#include <core/conformation/Atom.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/graph/Graph.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <ObjexxFCL/FArray3D.hh>
#include <utility/vector1.hh>
#include <map>
#include <numeric/interpolation/spline/CubicSpline.hh>

namespace core {
namespace scoring {
namespace membrane {


/// @brief Fullatom Membrane Solvation Energy - Statistically Derived,
/// and smoothed derivatives
class MPHelicalityEnergy : public core::scoring::methods::ContextDependentOneBodyEnergy  {

public:

	typedef ContextDependentOneBodyEnergy  parent;

public:
	// Constructors ////////////////////

	/// @brief Default Constructor
	MPHelicalityEnergy();

	/// @brief Create a clone of this energy method
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	// Scoring Methods ////////////////

	/// @brief Computes dScore/dNumNeighbors for all residues for rapid use in later
	/// atom derivate calculations
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const override;

	/// @brief Causes a neighbor graph update
	void
	setup_for_derivatives( pose::Pose &, ScoreFunction const & ) const override;

	/// @brief Evaluates the one-body energy for a residue
	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const &,
		EnergyMap &
	) const override;

	utility::vector1< core::Size >
	centroid_neighbors(
		pose::Pose const & pose,
		conformation::Residue const & rsd
	) const;

	/// @brief Increments the F1 and F2 derivative vectors for an atom
	void
	eval_atom_derivative(
		id::AtomID const & atom_id,
		pose::Pose const & pose,
		kinematics::DomainMap const &,
		ScoreFunction const &,
		EnergyMap const &,
		Vector &,
		Vector &
	) const override;

	/// @brief Unused by the MPHelicalityEnergy class, returns 0
	Distance
	atomic_interaction_cutoff() const;

	/// @brief Tells the scoring function to maintain the TwelveANeighborGraph
	void indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const override;

	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const override;

	/// @brief returns the atom name for the atom used to represent the sidechain for
	/// a particular amino acid; this atom was used to derive the statistics this potential
	/// is based on.
	// std::string const &
	// representative_atom_name( chemical::AA const aa ) const;

	/// @brief convert the neighbor counts for a residue and its aa type into
	/// a score and a score derivative.
	void
	calc_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		Real & score
	) const;

	// void
	// calc_energy_old(
	//  conformation::Residue const & rsd,
	//  pose::Pose const & pose,
	//  Real const neighbor_count,
	//  chemical::AA const aa,
	//  Real & score
	// ) const;
	/// @brief given the square distance between a representative atom and a neighbor atom,
	/// return the neighborlyness.  Ramps from 1 down to 0 over a range.
	// Real
	// sigmoidish_neighbor( DistanceSquared const sqdist ) const;

	/// @brief given a pair of atoms, one of which is a neighbor atom and the other of which
	/// is a representative atom, and given the weighted score derivative, increments
	/// the F1 and F2 derivatives.
	// void
	// increment_f1_f2_for_atom_pair(
	//  conformation::Atom const & atom1,
	//  conformation::Atom const & atom2,
	//  Real weighted_dScore_dN,
	//  Vector & F1,
	//  Vector & F2
	// ) const;

	utility::vector1 < core::Size >
	neighboring_atoms(
		pose::Pose const & pose,
		conformation::Residue const & rsd,
		Real const & cutoff_1,
		Real const & cutoff_2
	) const;

	Real
	burial_sigmoid(
		core::Size n_atoms,
		Real const slope,
		Real const offset
	) const;

	core::Real
	calc_residue_burial(
		pose::Pose const & pose,
		conformation::Residue const & rsd
	) const;


private:
	// mutable utility::vector1< Real > residue_N_;
	// mutable utility::vector1< Real > residue_E_;
	// mutable utility::vector1< Real > residue_dEdN_;
	// ObjexxFCL::FArray3D< Real > mem_env_log10_;

	core::Size version() const override;

};

} // membrane
} // scoring
} // core

#endif // INCLUDED_core_scoring_membrane_MPHelicalityEnergy_hh
