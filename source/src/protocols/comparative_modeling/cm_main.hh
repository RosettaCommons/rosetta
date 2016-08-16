// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/comparative_modeling/cm_main.hh
/// @author James Thompson

#ifndef INCLUDED_protocols_comparative_modeling_cm_main_hh
#define INCLUDED_protocols_comparative_modeling_cm_main_hh

namespace protocols {
namespace comparative_modeling {

/// @brief Initiates LoopRelaxThreadingMover using the job distributor (jd2)
void cm_main();

}
}

#endif  // INCLUDED_protocols_comparative_modeling_cm_main_hh
