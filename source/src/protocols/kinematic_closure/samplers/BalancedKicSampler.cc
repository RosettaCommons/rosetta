// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/internal.hh>
#include <protocols/kinematic_closure/utilities.hh>
#include <protocols/kinematic_closure/ClosureProblem.hh>
#include <protocols/kinematic_closure/ClosureSolution.hh>
#include <protocols/kinematic_closure/samplers/BalancedKicSampler.hh>
#include <protocols/kinematic_closure/perturbers/Perturber.hh>
#include <protocols/kinematic_closure/perturbers/PerturberSet.hh>
#include <protocols/kinematic_closure/perturbers/RamaPerturber.hh>
#include <protocols/kinematic_closure/pivot_pickers/PivotPicker.hh>
#include <protocols/kinematic_closure/pivot_pickers/StandardPivots.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>
#include <protocols/loop_modeling/loggers/Logger.hh>
#include <protocols/loop_modeling/loggers/NullLogger.hh>

// Utility headers
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <boost/foreach.hpp>

// C++ headers
#include <iostream>

#define foreach BOOST_FOREACH
using namespace std;

namespace protocols {
namespace kinematic_closure {
namespace samplers {

static numeric::random::RandomGenerator RG(48291);
using protocols::kinematic_closure::pivot_pickers::PivotPickerOP;
using protocols::loop_modeling::loggers::LoggerOP;

BalancedKicSampler::BalancedKicSampler() { // {{{1
	perturbers_ = new perturbers::PerturberSet;
	perturbers_->add(new perturbers::RamaPerturber);
	perturbers_->mark_as_default();

	pivot_picker_ = new pivot_pickers::StandardPivots;
	logger_ = new protocols::loop_modeling::loggers::NullLogger;
	setup_called_ = false;
}

BalancedKicSampler::~BalancedKicSampler() {} // {{{1

void BalancedKicSampler::setup(Pose & pose, Loop const & loop) { // {{{1
	setup_fold_tree(pose, loop);
	setup_called_ = true;
}

void BalancedKicSampler::apply(Pose & pose, Loop const & loop) { // {{{1
	ClosureProblemOP problem = new ClosureProblem(logger_);
	ClosureSolutionCOP solution;
	SolutionList unperturbed_solutions, perturbed_solutions;

	if (not setup_called_) setup(pose, loop);

	// Solve the unperturbed problem.
	problem->frame(pose, loop, pivot_picker_);
	problem->solve(unperturbed_solutions);

	// Solve the perturbed problem.
	perturbers_->perturb_with_balance(pose, problem);
	problem->solve(perturbed_solutions);

	// Pick a solution to return.
	solution = pick_solution(unperturbed_solutions, perturbed_solutions);
	solution->apply(pose);
}

// {{{1
/// @details Note that this is a static method, which means that it can be used 
/// outside the context of this class.  The inputs are two sets of closure 
/// solutions.  The first set should contain the starting conformation, and the 
/// second should contain novel conformations.  The first set can be generated 
/// by running the KIC algorithm without perturbing the non-pivot torsions.  
/// For each solution in these two sets, a jacobian will have to be calculated.  
/// The solution that is returned will be free of bias.  See apply() for an 
/// example of how this method is used.

ClosureSolutionCOP BalancedKicSampler::pick_solution (
		SolutionList const & unperturbed_solutions,
		SolutionList const & perturbed_solutions) {

	Real total_jacobian = 0;
	Real selection_chance = 0;
	Real random_threshold = numeric::random::uniform();
	ChainedSolutionList all_solutions(
			unperturbed_solutions, perturbed_solutions);

	if (all_solutions.empty()) return NULL;

	// Calculate the Jacobian for each solution and keep the sum.  The Jacobian 
	// relates to the probability that a certain set of pivot torsions will lead 
	// to a closed solution.

	foreach (ClosureSolutionCOP solution, all_solutions) {
		total_jacobian += solution->get_jacobian();
	}

	// Use the Jacobian weights to pick a balanced solution.  The pick must be 
	// made from both the perturbed and unperturbed pools of solutions, otherwise 
	// the forward and reverse move probabilities won't be equivalent.

	foreach (ClosureSolutionCOP solution, all_solutions) {
		selection_chance += solution->get_jacobian() / total_jacobian;
		if (selection_chance >= random_threshold)
			return solution;
	}

	// Execution will only get this far if random_threshold is very close to one 
	// and floating point error causes the sum of all the selection_chance terms 
	// to be slightly less than that.  In this case, the right course of action 
	// is clearly to return the last solution.

	return all_solutions.back();
}

void BalancedKicSampler::add_perturber(perturbers::PerturberOP perturber) { // {{{1
	perturbers_->add(perturber) ;
}

void BalancedKicSampler::set_pivot_picker(PivotPickerOP picker) { // {{{1
	pivot_picker_ = picker;
}

void BalancedKicSampler::log_filters(LoggerOP logger) { // {{{1
	logger_ = logger;
}
// }}}1

}
}
}
