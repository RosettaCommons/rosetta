// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_kinematic_closure_BalancedKicMover_HH
#define INCLUDED_protocols_kinematic_closure_BalancedKicMover_HH

// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/BalancedKicMover.fwd.hh>
#include <protocols/kinematic_closure/ClosureProblem.fwd.hh>
#include <protocols/kinematic_closure/ClosureSolution.fwd.hh>
#include <protocols/kinematic_closure/perturbers/Perturber.fwd.hh>
#include <protocols/kinematic_closure/perturbers/PerturberSet.fwd.hh>
#include <protocols/kinematic_closure/pivot_pickers/PivotPicker.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/id/TorsionID_Range.hh>

// Protocol headers
#include <protocols/canonical_sampling/ThermodynamicMover.hh>
#include <protocols/loops/Loop.hh>

// Utility headers
#include <utility/vector1.hh>
#include <boost/noncopyable.hpp>

namespace protocols {
namespace kinematic_closure {

/// @brief Make a kinematic closure move that obeys detailed balance.
///
/// @details This class is very conceptually similar to KicMover, so check 
/// out its documentation for a general overview of the kinematic closure 
/// algorithm.  Here I will just highlight some details associated with making 
/// a balanced version of the move.  Detailed balance is a useful property, 
/// because it allows a Monte Carlo simulation to recapitulate ensembles with 
/// correct equilibrium populations (so long as sampling is good, of course).  
/// There are a two reasons why the standard KicMover algorithm does not obey 
/// detailed balance.  The first is that the geometry of the closure move 
/// itself introduces some inherent bias which has to be explicitly canceled 
/// out.  The second is that care needs to be taken to perturb the non-pivot 
/// torsions in a way that also obeys detailed balance, as well.
///
/// The add_perturber() method works much like it does in KicMover.  The only 
/// conceptual difference is that when the added perturbers are used internally 
/// within apply(), perturbers::Perturber::perturb_with_balance() is called 
/// instead of perturbers::Perturber::perturb().  This makes it easy to make 
/// variants of the perturber algorithms which obey detailed balance.  The 
/// set_pivot_picker() method is no different from the KicMover version.

class BalancedKicMover
	: public protocols::canonical_sampling::ThermodynamicMover,
	  private boost::noncopyable {

public:

	/// @brief Default constructor.
	BalancedKicMover();

	/// @brief Default destructor.
	~BalancedKicMover();

public:

	/// @copydoc KicMover::apply
	void apply(Pose & pose);

	/// @copydoc KicMover::get_name
	string get_name() const { return "BalancedKicMover"; }

public:

	/// @copydoc KicMover::set_loop
	void set_loop(Loop const & loop);

	/// @copydoc KicMover::add_perturber
	void add_perturber(perturbers::PerturberOP perturber);

	/// @copydoc KicMover::set_pivot_picker
	void set_pivot_picker(pivot_pickers::PivotPickerOP picker);

public:

	/// @brief Return true, because this mover always obeys detailed balance.
	bool preserve_detailed_balance() const { return true; }

	/// @brief This mover always obeys detailed balance, so this is a no-op.
	void set_preserve_detailed_balance(bool) {}

	/// @details Right now the proposal probabilities are balanced internally, so 
	/// this ratio will always be unity.  This could change eventually, though.
	Real last_proposal_density_ratio() { return 1; }

	/// @brief Indicate that each torsion in the loop may take on any value.
	utility::vector1<core::id::TorsionID_Range> torsion_id_ranges(Pose & pose);

public:

	/// @brief Pick a solution in a way that cancels out the geometrical bias of 
	/// the kinematic closure algorithm.
	static ClosureSolutionCOP pick_solution(
			SolutionList const & unperturbed_solutions,
			SolutionList const & perturbed_solutions);


	/// @brief Return true if the given solution is the same as the input pose.  
	/// This allows for a more accurate reporting of Monte Carlo statistics.
	static bool is_solution_trivial(
			ClosureProblemCOP problem,
			ClosureSolutionCOP solution,
			SolutionList const & unperturbed_solutions,
			SolutionList const & perturbed_solutions);

private:
	bool is_fold_tree_stale_;
	protocols::loops::Loop loop_;
	perturbers::PerturberSetOP perturbers_;
	pivot_pickers::PivotPickerOP pivot_picker_;

};

}
}

#endif
