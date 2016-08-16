// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/medal/MedalMain.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_PROTOCOLS_MEDAL_MEDALMAIN_HH
#define INCLUDED_PROTOCOLS_MEDAL_MEDALMAIN_HH

namespace protocols {
namespace medal {

/// @brief Entry point for Medal protocol
void* Medal_main(void*);

/// @brief Entry point for MedalExchange protocol
void* MedalExchange_main(void*);

}  // namespace medal
}  // namespace protocols

#endif  // PROTOCOLS_MEDAL_MEDAL_MAIN_HH_
