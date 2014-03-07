// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief  Source file for ClosureSolution.
/// @author Kale Kundert (kale.kundert@ucsf.edu)

// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/internal.hh>
#include <protocols/kinematic_closure/ClosureProblem.hh>
#include <protocols/kinematic_closure/ClosureSolution.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoringManager.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/constants.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <numeric/kinematic_closure/closure.hh>

#include <numeric/angle.functions.hh>

// Utility headers
#include <utility/exit.hh>
#include <boost/foreach.hpp>

// External headers
#include <Eigen/Dense>

// C++ headers
#include <algorithm>
#include <cmath>

using namespace std;

namespace protocols {
namespace kinematic_closure {

using namespace std;
using core::Size;
using core::Real;
using core::pose::Pose;
using numeric::constants::r::pi;

// {{{1
/// @details This constructor is called by ClosureProblem, which is a friend 
/// class.  It is declared as private because it should only be called from 
/// code that has been specifically written to solve a closure problem.

ClosureSolution::ClosureSolution(
		ClosureProblemCOP problem, Size const index,
		ParameterList const & torsion_angles,
		ParameterList const & bond_angles,
		ParameterList const & bond_lengths)

		: problem_(problem), index_(index),
		  bond_lengths_(bond_lengths),
		  bond_angles_(bond_angles),
		  torsion_angles_(torsion_angles),
		  jacobian_(-1) {}

void ClosureSolution::apply(Pose & pose) const { // {{{1
	problem_->apply_internal_coordinates(
			bond_lengths_, bond_angles_, torsion_angles_, pose);
}

bool ClosureSolution::apply_if_reasonable( // {{{1
		Pose & pose, bool rama_on, bool bump_on, bool be_lenient) const {

	// It's very important to perform the rama check before updating the pose.  
	// If the pose is updated before the rama check, the rama check will cause 
	// the conformation to update its coordinates, which becomes a bottleneck.  
	// Since so many solutions fail the rama check, it is much better just to 
	// update the pose afterwards.  Note that this means the pose passed into the 
	// rama check is expected to have out-of-date coordinates!

	Real const temperature = be_lenient ? 2.0 : 1.0;

	if (rama_on and not check_rama(pose, temperature)) {
		num_rama_filter_fails += 1;
		return false;
	}

	// On the other hand, the solution must be applied to the pose before the 
	// bump check is performed, because the bump check uses only the cartesian 
	// coordinates in the given pose.  If the solution is not applied before 
	// calling this method, the check will meaninglessly act on whatever is 
	// already in the pose, even though the check is a method of the current 
	// solution.

	apply(pose);

	Real const scale_factor = be_lenient ? 0.40 : 0.49;

	if (bump_on and not check_overlap(pose, scale_factor)) {
		problem_->restore(pose);
		num_bump_filter_fails += 1;
		return false;
	}

	return true;
}

bool ClosureSolution::check_rama( // {{{1
		Pose const & pose, Real const temperature) const {

	using core::chemical::AA;
	using core::scoring::Ramachandran;
	using core::scoring::ScoringManager;
	using numeric::random::uniform;
	using numeric::conversions::degrees;

	Ramachandran const & rama =
		ScoringManager::get_instance()->get_Ramachandran();

	for (Size i = 1; i <= 3; i++) {
		Size ca = problem_->pivot_atoms()[i];
		AA type = pose.aa(problem_->pivot_residues()[i]);

		Real new_phi = degrees(torsion_angles_[ca - 1]);
		Real new_psi = degrees(torsion_angles_[ca]);
		Real new_score = rama.eval_rama_score_residue(type, new_phi, new_psi);

		// If we get the maximum possible rama score, bail out immediately.
		if (new_score >= 20.0) return false;

		Real old_phi = degrees(problem_->unperturbed_torsions_[ca - 1]);
		Real old_psi = degrees(problem_->unperturbed_torsions_[ca]);
		Real old_score = rama.eval_rama_score_residue(type, old_phi, old_psi);

		// Apply the Metropolis criterion to decide whether or not the new rama 
		// score should be accepted.  This seems like a pretty arbitrary way to 
		// make a decision, but performance is much worse with a fixed cutoff.

		if (new_score > old_score) {
			Real const difference = old_score - new_score;
			Real const probability = exp(difference / temperature);
			if (uniform() >= probability) return false;
		}
	}

	return true;
}

// {{{1
/// @details Apply the solution to the pose before calling this filter.
bool ClosureSolution::check_overlap(
		Pose const & pose, Real const scale_factor) const {

	using core::Vector;
	using core::conformation::Residue;

	Size first_residue = problem_->first_residue();
	Size last_residue = problem_->last_residue();

	// Iterate over loop residues.
	for (Size i = first_residue; i <= last_residue; i++) {
		Residue const & residue_i = pose.residue(i);
		Vector const & vector_i = residue_i.xyz(residue_i.nbr_atom());

		// Iterate over all other residues in the protein.
		for (Size j = 1; j <= pose.total_residue(); j++ ) {

			// Don't do adjacent residues.
			if ((j == i) || (j == i+1) || (j == i-1)) continue;
			
			// Don't do loop residues multiple times.
			if ((j >= first_residue) && (j <= i)) continue;

			Residue const & residue_j = pose.residue(j);
			Vector const & vector_j = residue_j.xyz(residue_j.nbr_atom());

			// Determine the neighbor cutoff from the radii of the neighbor atoms.
			Real const nbr_cutoff = residue_i.nbr_radius() + residue_j.nbr_radius();
			Real const nbr_cutoff_sq = nbr_cutoff * nbr_cutoff;
			Real const nbr_distance_sq = (vector_i - vector_j).length_squared();

			if (nbr_distance_sq > nbr_cutoff_sq) continue;

			// Check for clashes between the N, CA, C, O, and CB (except for glycine) 
			// atoms of the two residues.
			Size num_atoms_i = min<Size>(5, residue_i.nheavyatoms());

			for (Size m = 1; m <= num_atoms_i; m++) {
				Size num_atoms_j = residue_j.is_protein() ?
					min<Size>(5, residue_j.nheavyatoms()) : residue_j.nheavyatoms();

				for (Size n = 1; n <= num_atoms_j; n++) {
					Vector const & atom_i = residue_i.xyz(m);
					Vector const & atom_j = residue_j.xyz(n);

					Real const radius_i = residue_i.atom_type(m).lj_radius();
					Real const radius_j = residue_j.atom_type(n).lj_radius();

					// Compare the squared sum of Lennard-Jones radii.
					Real const bump_cutoff = radius_i + radius_j;
					Real const bump_cutoff_sq = scale_factor * bump_cutoff * bump_cutoff;
					Real const bump_distance_sq = (atom_i - atom_j).length_squared();

					if (bump_distance_sq < bump_cutoff_sq) return false;
				}
			}
		}
	}

	return true;
}


Size ClosureSolution::get_index() const { // {{{1
	return index_;
}

// {{{1
/// @details This quantity indicates how much the dihedral space around the 
/// pivots was warped by the choice of controls and is used as a normalization 
/// factor when picking a move in such a way that obeys detailed balance.  The 
/// return value is cached, so it's cheap to call this function multiple times.

Real ClosureSolution::get_jacobian() const {
	using std::abs;
	using numeric::max;

	if (jacobian_ < 0) {
		Eigen::Matrix<Real, 6, 3> r1, r2, s1, s2, gamma;
		Eigen::Matrix<Real, 3, 1> cross_i, cross_45, delta;
		Eigen::Matrix<Real, 4, 4> J;

		// Convert the solution from internal coordinates to  cartesian ones.  
		// Since the jacobian is invariant under rigid-body transformations, the 
		// origin and reference frame used in this conversion are unimportant.

		Coordinate dummy_origin;
		CoordinateList dummy_frame, atom_xyzs;

		numeric::kinematic_closure::radians::chainXYZ(
				problem_->num_atoms(),
				bond_lengths_,
				bond_angles_,
				torsion_angles_,
				atom_xyzs);

		for (Size i = 1; i <= 3; i++) {
			Size pivot = problem_->pivot_atoms()[i];
			Size j = 2 * (i - 1);

			r1(j, 0) = atom_xyzs[pivot-1][1];
			r1(j, 1) = atom_xyzs[pivot-1][2];
			r1(j, 2) = atom_xyzs[pivot-1][3];

			r1(j+1, 0) = atom_xyzs[pivot][1];
			r1(j+1, 1) = atom_xyzs[pivot][2];
			r1(j+1, 2) = atom_xyzs[pivot][3];

			r2(j, 0) = atom_xyzs[pivot][1];
			r2(j, 1) = atom_xyzs[pivot][2];
			r2(j, 2) = atom_xyzs[pivot][3];

			r2(j+1, 0) = atom_xyzs[pivot+1][1];
			r2(j+1, 1) = atom_xyzs[pivot+1][2];
			r2(j+1, 2) = atom_xyzs[pivot+1][3];
		}

		// Calculate the jacobian following the method outlined by Nilmeier, Hua, 
		// Coutsias, and Jacobson in their 2011 JCTC paper.

		for (Size i = 0; i < 6; i++) {
			delta = r2.row(i) - r1.row(i);
			gamma.row(i) = delta.normalized();
		}

		cross_45 = gamma.row(4).cross(gamma.row(5));

		for (Size i = 0; i < 4; i++) {
			cross_i = gamma.row(i).cross(r1.row(5) - r1.row(i));
			Real dot_i = gamma.row(i).dot(cross_45);

			J(0, i) = cross_i(0);
			J(1, i) = cross_i(1);
			J(2, i) = cross_i(2);
			J(3, i) = dot_i;
		}

		jacobian_ = 1 / max(abs(J.determinant()), 1e-100);
	}

	return jacobian_;
}

// {{{1
/// @details Note that this is not a rigorous distance metric.  It's just meant 
/// to distinguish one solution that's nearly identical to the given problem 
/// from several solutions that aren't.
Real ClosureSolution::get_distance(ClosureProblemCOP problem) const {
	Real distance = 0;

	for (Size i = 1; i <= problem->num_atoms(); i++) {
		distance += pow(
				(bond_lengths_[i] - problem->unperturbed_lengths_[i]) / 1.48, 2);
		distance += 0.5 - 0.5 * cos(
				bond_angles_[i] - problem->unperturbed_angles_[i]);
		distance += 0.5 - 0.5 * cos(
				torsion_angles_[i] - problem->unperturbed_torsions_[i]);
	}

	return distance;
}
// }}}1

} // namespace kinematic_closure
} // namespace protocols

