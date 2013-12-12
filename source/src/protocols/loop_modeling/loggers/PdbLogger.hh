// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_loop_modeling_loggers_PdbLogger_HH
#define INCLUDED_protocols_loop_modeling_loggers_PdbLogger_HH

// Unit headers
#include <protocols/loop_modeling/loggers/PdbLogger.fwd.hh>
#include <protocols/loop_modeling/loggers/Logger.hh>

// Protocol headers
#include <protocols/moves/MonteCarlo.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

namespace protocols {
namespace loop_modeling {
namespace loggers {

using core::Size;
using core::pose::Pose;

class PdbLogger : public Logger {

public:

	PdbLogger();
	PdbLogger(bool log_tasks);

protected:

	void log_beginning_(Pose const & pose);
	void log_iteration_(Pose const & pose);
	void log_task_(Pose const & pose, string name, bool successful);
	void log_monte_carlo_(MonteCarlo const & monte_carlo);

private:

	void log_pose(Pose const & pose) const;
	void log_pose(Pose const & pose, string task) const;
	void log_pose(Pose const & pose, Size iteration) const;
	void log_pose(Pose const & pose, Size iteration, string task) const;

private:

	Size counter_;
	bool log_tasks_;

};

}
}
}

#endif

