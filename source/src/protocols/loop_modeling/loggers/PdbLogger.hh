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
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/loggers/Logger.hh>
#include <protocols/loop_modeling/loggers/PdbLogger.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Protocol headers
#include <protocols/moves/MonteCarlo.fwd.hh>

namespace protocols {
namespace loop_modeling {
namespace loggers {

/// @brief Record snapshots of the simulation in the PDB format.
class PdbLogger : public Logger {

protected:

	/// @brief Record the input pose.
	void log_beginning_(Pose const & pose);

	/// @brief Record the pose after each Monte Carlo step.
	void log_monte_carlo_(protocols::moves::MonteCarlo const & monte_carlo);

private:

	/// @brief Write the given pose to disk. The file name will incorporate the 
	/// current iteration.
	void log_pose(Pose const & pose) const;

	/// @brief Write the given pose to disk. The file name will incorporate the 
	/// given iteration.
	void log_pose(Pose const & pose, Size iteration) const;

};

}
}
}

#endif
