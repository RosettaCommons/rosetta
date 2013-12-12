// Unit headers
#include <protocols/loop_modeling/filters/BumpCheck.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunction.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>

// Utility headers
#include <numeric/xyzVector.hh>
#include <boost/foreach.hpp>

// C++ headers
#include <cmath>
#include <algorithm>

#define foreach BOOST_FOREACH

namespace protocols {
namespace loop_modeling {
namespace filters {

using namespace std;

bool BumpCheck::apply(
		Pose & pose, Loop const & loop, ScoreFunctionCOP score_function) {

	using core::Size;
	using core::Real;
	using core::Vector;
	using core::conformation::Residue;

	Size first_residue = loop.start();
	Size last_residue = loop.stop();

	// Iterate over loop residues.
	for (Size i = first_residue; i <= last_residue; i++) {
		Residue const & residue_i = pose.residue(i);
		Vector const & vector_i = residue_i.xyz(residue_i.nbr_atom());

		// Iterate over all other residues in the protein.
		for (Size j = 1; j <= pose.total_residue(); j++ ) {

			// Don't do same or adjacent residues.
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
					Real const bump_cutoff_sq = 0.49 * bump_cutoff * bump_cutoff;
					Real const bump_distance_sq = (atom_i - atom_j).length_squared();

					if (bump_distance_sq < bump_cutoff_sq) return false;
				}
			}
		}
	}

	return true;
}

} // namespace filters
} // namespace kinematic_closure
} // namespace protocols

