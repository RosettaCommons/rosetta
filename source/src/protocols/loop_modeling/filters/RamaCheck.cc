// Unit headers
#include <protocols/loop_modeling/filters/RamaCheck.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/Ramachandran.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>

// Utility headers
#include <boost/foreach.hpp>

#define foreach BOOST_FOREACH

namespace protocols {
namespace loop_modeling {
namespace filters {

using namespace std;

// This function used to work by comparing old and new rama scores for each 
// pivot residue.  But KIC moves by their nature have very little memory, so it 
// seemed rather arbitrary to compare against the previous state.  I think it 
// would be better to compare against a fixed threshold, but it will take some 
// effort to figure out what that threshold should be.

bool RamaCheck::apply(
		Pose & pose, Loop const & loop, ScoreFunctionCOP score_function) {

	using core::scoring::Ramachandran;
	using core::scoring::ScoringManager;
	using core::conformation::Residue;

	/*
	Ramachandran const & rama =
		ScoringManager::get_instance()->get_Ramachandran();

	Size pivot_residues[] = {
		loop.start(), loop.cut(), loop.stop()};

	foreach (Size i, pivot_residues) {
		Residue const & residue = pose.residue(i);
		Real score = rama.eval_rama_score_residue(residue);

		// Return false below some cutoff.
	}
	*/

	return true;
}

} // namespace filters
} // namespace kinematic_closure
} // namespace protocols

