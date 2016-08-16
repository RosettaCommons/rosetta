// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file relax_initialization_protocols
/// @brief
/// @details
///   For diagnosis. Run with this TrialMover and get a log-file containing output from the provided PoseEvalutor
///    and if it was accepted or not.
///
///
/// @author Oliver Lange


#ifndef INCLUDED_protocols_simple_filters_EvaluatedTrialMover_hh
#define INCLUDED_protocols_simple_filters_EvaluatedTrialMover_hh


// Unit Headers


// Package Headers
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>

// Project Headers
#include <core/io/silent/silent.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

//// C++ headers
#include <list>

#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace simple_filters {


class EvaluatedTrialMover;
typedef  utility::pointer::shared_ptr< EvaluatedTrialMover >  EvaluatedTrialMoverOP;


class EvaluatedTrialMover : public moves::TrialMover {
public:
	/// c'stor
	EvaluatedTrialMover(
		moves::MoverOP mover_in,
		moves::MonteCarloOP mc_in,
		evaluation::PoseEvaluatorOP evaluator_in,
		std::string tag
	);
	~EvaluatedTrialMover();

	/// apply does a single trial (which is a mover apply and a boltzmann)
	virtual void apply( core::pose::Pose & pose );

	virtual std::string get_name() const;

	/// write all collected output to file (appends if file exist )
	void dump_file( std::string fn );
private:
	evaluation::MetaPoseEvaluatorOP evaluator_;
	typedef utility::vector1< core::io::silent::SilentStructOP > SilentInfoList;
	SilentInfoList evals_;
	std::string tag_;  //the decoy_output tag
};


}
}

#endif
