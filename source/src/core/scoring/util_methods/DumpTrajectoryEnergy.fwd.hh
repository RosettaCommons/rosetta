// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/util_methods/DumpTrajectoryEnergy.fwd.hh
/// @brief Forward declarations for an EnergyMethod that dumps trajectories to file.
/// @details Dumps trajectories of the minimizer and packer to file when the dump_trajectory
/// ScoreType is enable. Output may be controlled through the dump_trajectory:* flags.
/// @author Brian Coventry (bcov@uw.edu).


#ifndef INCLUDED_core_scoring_util_methods_DumpTrajectoryEnergy_fwd_hh
#define INCLUDED_core_scoring_util_methods_DumpTrajectoryEnergy_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace util_methods {

class DumpTrajectoryEnergy;

typedef utility::pointer::shared_ptr< DumpTrajectoryEnergy > DumpTrajectoryEnergyOP;


} // util_methods
} // scoring
} // core


#endif
