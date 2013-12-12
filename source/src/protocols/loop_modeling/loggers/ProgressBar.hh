// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_loop_modeling_loggers_ProgressBar_HH
#define INCLUDED_protocols_loop_modeling_loggers_ProgressBar_HH

// Unit headers
#include <protocols/loop_modeling/loggers/ProgressBar.fwd.hh>
#include <protocols/loop_modeling/loggers/Logger.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace loop_modeling {
namespace loggers {

using core::Size;
using core::pose::Pose;

class ProgressBar : public Logger {

public:

	void log_iteration_(Pose const & pose);
	void log_ending_(Pose const & pose);

};

}
}
}

#endif

