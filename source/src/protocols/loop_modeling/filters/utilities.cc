// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/internal.hh>
#include <protocols/kinematic_closure/ClosureProblem.hh>
#include <protocols/kinematic_closure/ClosureSolution.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/chemical/AtomType.hh>
#include <core/id/NamedAtomID.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/Ramachandran.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

// Utility headers
#include <boost/foreach.hpp>

// External headers
#include <Eigen/Dense>

// C++ headers
#include <cmath>
#include <algorithm>

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/// This file contains source code for the rama and bump check.  It's just  ///
/// meant to be a reference; the real implementations are in the filter     ///
/// classes.  Once they have been tested, this file can be safely removed.  ///
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

#define foreach BOOST_FOREACH

namespace protocols {
namespace kinematic_closure {

using namespace std;
using core::Size;
using core::Real;
using core::pose::PoseCOP;
static numeric::random::RandomGenerator RG(49297);

// bool ClosureSolution::check_rama(Real temperature) const // {{{1

/// @details This function works by comparing old and new rama scores for each 
/// pivot residue.  But KIC moves by their nature have very little memory, so 
/// it seems rather arbitrary to compare against the previous state.  That 
/// said, the only other way I can imagine doing this is to set a rama 
/// threshold.  That would probably be even worse.

bool ClosureSolution::check_rama(Real temperature) const {
	using core::scoring::Ramachandran;
	using core::scoring::ScoringManager;
	using core::conformation::Residue;

	PoseCOP new_pose = pose_;
	PoseCOP old_pose = problem_->closed_pose_;

	Ramachandran const & rama =
		ScoringManager::get_instance()->get_Ramachandran();

	foreach (Size i, problem_->pivot_residues_) {
		Residue const & new_residue = new_pose->conformation().residue(i);
		Residue const & old_residue = old_pose->conformation().residue(i);

		Real new_score = rama.eval_rama_score_residue(new_residue);
		Real old_score = rama.eval_rama_score_residue(old_residue);

		// Allow bad moves occasionally.
		if ( new_score > old_score ) {
			Real const difference = old_score - new_score;
			Real const probability = std::exp(difference / temperature);

			if (RG.uniform() >= probability) return false;
		}
	}

	return true;
}

bool ClosureSolution::check_overlap() const { // {{{1
	using core::Vector;
	using core::conformation::Residue;

	Size first_residue = problem_->first_residue_;
	Size last_residue = problem_->last_residue_;

	// Iterate over loop residues.
	for (Size i = first_residue; i <= last_residue; i++) {
		Residue const & residue_i = pose_->conformation().residue(i);
		Vector const & vector_i = residue_i.xyz(residue_i.nbr_atom());

		// Iterate over all other residues in the protein.
		for (Size j = 1; j <= pose_->total_residue(); j++ ) {

			// Don't do same or adjacent residues.
			if ((j == i) || (j == i+1) || (j == i-1)) continue;
			
			// Don't do loop residues multiple times.
			if ((j >= first_residue) && (j <= i)) continue;

			Residue const & residue_j = pose_->conformation().residue(j);
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
					Real const bump_cutoff_sq = 0.49 * bump_cutoff * bump_cutoff;
					Real const bump_distance_sq = (atom_i - atom_j).length_squared();

					if (bump_distance_sq < bump_cutoff_sq) return false;
				}
			}
		}
	}

	return true;
}

// }}}1
// }}}1

} // end namespace kinematic_closure
} // end namespace protocols
