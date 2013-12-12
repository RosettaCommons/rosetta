// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_loop_modeling_LoopMover_HH
#define INCLUDED_protocols_loop_modeling_LoopMover_HH

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopMover.fwd.hh>
#include <protocols/loop_modeling/LoopMoverTask.hh>
#include <protocols/loop_modeling/loggers/Logger.hh>

// Core headers
#include <core/scoring/ScoreFunction.hh>

// Protocols headers
#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/loops/Loop.hh>

// Utility headers
#include <utility/vector1.hh>
#include <boost/utility.hpp>

namespace protocols {
namespace loop_modeling {

class LoopMover : public protocols::moves::Mover, private boost::noncopyable {

public:
	LoopMover();
	LoopMover(Loop const & loop, ScoreFunctionOP score_function);

	void apply(Pose & pose);
	void setup(Pose & pose);
	virtual string get_name() const { return "LoopMover"; }

public:
	void set_loop(Loop const & loop);
	void set_score_function(ScoreFunctionOP score_function);
	void add_task(LoopMoverTaskOP task);
	void add_filter(protocols::filters::FilterOP filter);
	void add_logger(loggers::LoggerOP logger);

	Loop get_loop() const;
	ScoreFunctionOP get_score_function() const;
	bool was_successful() const;

private:
	Loop loop_;
	ScoreFunctionOP score_function_;
	utility::vector1<LoopMoverTaskOP> tasks_;
	loggers::LoggerList loggers_;

	bool setup_needed_;
	bool last_move_successful_;

};

}
}

#endif

