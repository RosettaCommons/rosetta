// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_kinematic_closure_BalancedKicSampler_HH
#define INCLUDED_protocols_kinematic_closure_BalancedKicSampler_HH

// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/ClosureSolution.fwd.hh>
#include <protocols/kinematic_closure/samplers/BalancedKicSampler.fwd.hh>
#include <protocols/kinematic_closure/perturbers/Perturber.fwd.hh>
#include <protocols/kinematic_closure/perturbers/PerturberSet.fwd.hh>
#include <protocols/kinematic_closure/pivot_pickers/PivotPicker.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>
#include <protocols/loop_modeling/loggers/Logger.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <boost/noncopyable.hpp>

namespace protocols {
namespace kinematic_closure {
namespace samplers {

/// @brief Make a kinematic closure move that obeys detailed balance.
///
/// @details This class is very conceptually similar to KicSampler, so check 
/// out its documentation for a general overview of the kinematic closure 
/// algorithm.  Here I will just highlight some details associated with making 
/// a balanced version of the move.  Detailed balance is a useful property, 
/// because it allows a Monte Carlo simulation to recapitulate ensembles with 
/// correct equilibrium populations (so long as sampling is good, of course).  
/// There are a two reasons why the standard KicSampler algorithm does not obey 
/// detailed balance.  The first is that the geometry of the closure move 
/// itself introduces some inherent bias which has to be explicitly canceled 
/// out.  The second is that care needs to be taken to perturb the non-pivot 
/// torsions in a way that also obeys detailed balance, as well.
///
/// The add_perturber() method works much like it does in KicSampler.  The only 
/// conceptual difference is that when the added perturbers are used internally 
/// within apply(), perturbers::Perturber::perturb_with_balance() is called 
/// instead of perturbers::Perturber::perturb().  This makes it easy to make 
/// variants of the perturber algorithms which obey detailed balance.  The 
/// set_pivot_picker() method is no different from the KicSampler version.

class BalancedKicSampler
	: public utility::pointer::ReferenceCount, private boost::noncopyable {

public:

	/// @brief Default constructor.
	BalancedKicSampler();

	/// @brief Destructor.
	~BalancedKicSampler();

public:

	/// @copydoc KicSampler::setup
	void setup(Pose & pose, Loop const & loop);

	/// @copydoc KicSampler::apply
	void apply(Pose & pose, Loop const & loop);

public:

	/// @copydoc KicSampler::add_perturber
	void add_perturber(perturbers::PerturberOP perturber);

	/// @copydoc KicSampler::set_pivot_picker
	void set_pivot_picker(pivot_pickers::PivotPickerOP picker);

	/// @copydoc KicSampler::log_filters
	void log_filters(protocols::loop_modeling::loggers::LoggerOP logger);

public:

	/// @brief Pick a solution in a way that cancels out the geometrical bias of 
	/// the kinematic closure algorithm.
	static ClosureSolutionCOP pick_solution(
			SolutionList const & unperturbed_solutions,
			SolutionList const & perturbed_solutions);

private:
	bool setup_called_;
	perturbers::PerturberSetOP perturbers_;
	pivot_pickers::PivotPickerOP pivot_picker_;
	protocols::loop_modeling::loggers::LoggerOP logger_;

};

}
}
}

#endif
