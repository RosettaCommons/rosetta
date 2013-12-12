// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_loop_modeling_LoopProtocol_HH
#define INCLUDED_protocols_loop_modeling_LoopProtocol_HH

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopProtocol.fwd.hh>
#include <protocols/loop_modeling/LoopMover.fwd.hh>
#include <protocols/loop_modeling/loggers/Logger.fwd.hh>

// Protocols headers
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/loops/Loop.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility headers
#include <boost/utility.hpp>

namespace protocols {
namespace loop_modeling {

class LoopProtocol
	: public protocols::moves::Mover, private boost::noncopyable {

	public:
		LoopProtocol();
		~LoopProtocol();

		void apply(Pose & pose);
		string get_name() const { return "LoopProtocol"; }

	private:
		void start_protocol(Pose & pose);
		void ramp_temperature(Size iteration);
		void ramp_score_function(Size iteration);
		void attempt_loop_move(Pose & pose, Size i, Size j, Size k);
		void finish_protocol(Pose const & pose);

	public:
		void set_mover(LoopMoverOP value);
		void set_loop(Loop const & value);
		void set_score_function(ScoreFunctionOP function);
		void set_iterations(IndexList values);
		void set_iterations(Size i, Size j, Size k);
		void set_temperature_schedule(Real initial, Real final);
		void set_temperature_ramping(bool value);
		void set_score_function_ramping(bool value);
		void add_logger(loggers::LoggerOP logger);

	private:
		Loop loop_;
		LoopMoverOP mover_;
		ScoreFunctionOP score_function_;
		protocols::moves::MonteCarloOP monte_carlo_;
		IndexList iterations_;
		loggers::LoggerList loggers_;

		bool ramp_score_function_;
		bool ramp_temperature_;

		Real initial_temperature_;
		Real final_temperature_;
		Real temperature_scale_factor_;
};

}
}

#endif

