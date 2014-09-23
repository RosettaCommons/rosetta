// Headers {{{1
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/internal.hh>
#include <protocols/kinematic_closure/BalancedKicMover.hh>
#include <protocols/kinematic_closure/ClosureProblem.hh>
#include <protocols/kinematic_closure/ClosureSolution.hh>
#include <protocols/kinematic_closure/perturbers/Perturber.hh>
#include <protocols/kinematic_closure/perturbers/PerturberSet.hh>
#include <protocols/kinematic_closure/perturbers/RamaPerturber.hh>
#include <protocols/kinematic_closure/pivot_pickers/PivotPicker.hh>
#include <protocols/kinematic_closure/pivot_pickers/StandardPivots.hh>

// Core headers
#include <core/id/types.hh>
#include <core/id/TorsionID.hh>
#include <core/id/TorsionID_Range.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/loops_main.hh>

// Utility headers
#include <utility/exit.hh>
#include <numeric/constants.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <boost/foreach.hpp>

// C++ headers
#include <iostream>
#include <limits>

// Global Names {{{1

using namespace std;
using core::id::TorsionID;
using core::id::TorsionID_Range;
using protocols::kinematic_closure::pivot_pickers::PivotPickerOP;
using utility::vector1;
// }}}1

namespace protocols {
namespace kinematic_closure {

// Member functions:

BalancedKicMover::BalancedKicMover() { // {{{1
	perturbers_ = new perturbers::PerturberSet;
	perturbers_->add(perturbers::PerturberOP( new perturbers::RamaPerturber ));
	perturbers_->mark_as_default();

	loop_ = Loop(0, 0);
	pivot_picker_ = new pivot_pickers::StandardPivots;
	is_fold_tree_stale_ = true;
}

BalancedKicMover::~BalancedKicMover() {} // {{{1

void BalancedKicMover::apply(Pose & pose) { // {{{1
	ClosureProblemOP problem = new ClosureProblem();
	ClosureSolutionCOP solution;
	SolutionList unperturbed_solutions, perturbed_solutions;

	if (loop_.start() == 0 && loop_.stop() == 0) {
		utility_exit_with_message(
				"Before calling BalancedKicMover.apply(), you must provide a loop "
				"via BalancedKicMover.set_loop().");
	}
	if (is_fold_tree_stale_) {
		protocols::loops::set_single_loop_fold_tree(pose, loop_);
	}

	// Solve the unperturbed problem.
	problem->frame(pose, loop_, pivot_picker_);
	unperturbed_solutions = problem->solve();

	// Solve the perturbed problem.
	perturbers_->perturb_with_balance(pose, problem);
	perturbed_solutions = problem->solve();

	// Pick a solution to apply.
	solution = pick_solution(unperturbed_solutions, perturbed_solutions);
	solution->apply(pose);

	// Decide if a real move was made.
	bool is_trivial = is_solution_trivial(
			problem, solution, unperturbed_solutions, perturbed_solutions);
	type(is_trivial ? "balanced-kic-no-op" : "balanced-kic");
}

void BalancedKicMover::set_loop(Loop const & loop) { // {{{1
	loop_ = loop;
	is_fold_tree_stale_ = true;
}

void BalancedKicMover::add_perturber(perturbers::PerturberOP perturber) { // {{{1
	perturbers_->add(perturber) ;
}

void BalancedKicMover::set_pivot_picker(PivotPickerOP picker) { // {{{1
	pivot_picker_ = picker;
}

vector1<TorsionID_Range> BalancedKicMover::torsion_id_ranges(Pose &) { // {{{1
	using core::id::BB;

	vector1<TorsionID_Range> results;
	Real static const pi = numeric::constants::r::pi;

	for (Size residue = loop_.start(); residue <= loop_.stop(); residue++) {
		TorsionID phi(residue, BB, 1);
		TorsionID psi(residue, BB, 2);
		TorsionID omega(residue, BB, 3);

		results.push_back(TorsionID_Range(phi, -pi, pi));
	}

	return results;
}
// }}}1

// Static member functions:

// {{{1
/// @details Note that this is a static method, which means that it can be used
/// outside the context of this class.  The inputs are two sets of closure
/// solutions.  The first set should contain the starting conformation, and the
/// second should contain novel conformations.  The first set can be generated
/// by running the KIC algorithm without perturbing the non-pivot torsions.
/// For each solution in these two sets, a jacobian will have to be calculated.
/// The solution that is returned will be free of bias.  See apply() for an
/// example of how this method is used.

ClosureSolutionCOP BalancedKicMover::pick_solution(
		SolutionList const & unperturbed_solutions,
		SolutionList const & perturbed_solutions) {

	Real total_jacobian = 0;
	Real selection_chance = 0;
	Real random_threshold = numeric::random::uniform();
	ChainedSolutionList all_solutions(
			unperturbed_solutions, perturbed_solutions);

	runtime_assert(! all_solutions.empty())

	// Calculate the Jacobian for each solution and keep the sum.  The Jacobian
	// relates to the probability that a certain set of pivot torsions will lead
	// to a closed solution.

	BOOST_FOREACH(ClosureSolutionCOP solution, all_solutions) {
		total_jacobian += solution->get_jacobian();
	}

	// Use the Jacobian weights to pick a balanced solution.  The pick must be
	// made from both the perturbed and unperturbed pools of solutions, otherwise
	// the forward and reverse move probabilities won't be equivalent.

	BOOST_FOREACH(ClosureSolutionCOP solution, all_solutions) {
		selection_chance += solution->get_jacobian() / total_jacobian;
		if (selection_chance >= random_threshold) return solution;
	}

	// Execution will only get this far if random_threshold is very close to one
	// and floating point error causes the sum of all the selection_chance terms
	// to be slightly less than that.  In this case, the right course of action
	// is clearly to return the last solution.

	return all_solutions.back();
}

// {{{1
/// @details In order to obey detailed balance, the set of unperturbed
/// solutions must contain one solution that is identical to the input pose.
/// When this solution is picked and applied, the resulting move will pass the
/// Metropolis criterion and lead to an artificially inflated acceptance rate.
/// This method provides a way to report what really happened in the move, so
/// that an accurate acceptance rate can be conveyed.  Note that this is a
/// static method, so it can be used outside the context of this class.  The
/// input pose conformation is inferred from the given problem.

bool BalancedKicMover::is_solution_trivial(
		ClosureProblemCOP problem,
		ClosureSolutionCOP picked_solution,
		SolutionList const & unperturbed_solutions,
		SolutionList const & perturbed_solutions)
{
	Real distance, closest_distance = numeric_limits<Real>::infinity();
	ClosureSolutionCOP closest_solution = NULL;
	ChainedSolutionList all_solutions(
			unperturbed_solutions, perturbed_solutions);

	BOOST_FOREACH(ClosureSolutionCOP solution, all_solutions) {
		distance = solution->get_distance(problem.get());
		if (distance < closest_distance) {
			closest_distance = distance;
			closest_solution = solution;
		}
	}

	runtime_assert(closest_solution.get() != NULL);
	return picked_solution.get() == closest_solution.get();
}
// }}}1

}
}
