// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/random/RandomGeneratorSettings.fwd.hh
/// @brief Forward declarations for a class to store settings from the options system for the random generator.  Moved from core to basic.
/// @author Original author unknown.
/// @modified Moved from core to basic by Vikram K. Mulligan (vmulligan@flatironinstitute.org).

#ifndef INCLUDED_basic_random_RandomGeneratorSettings_fwd_hh
#define INCLUDED_basic_random_RandomGeneratorSettings_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace basic {
namespace random {

class RandomGeneratorSettings;

typedef utility::pointer::shared_ptr< RandomGeneratorSettings > RandomGeneratorSettingsOP;
typedef utility::pointer::shared_ptr< RandomGeneratorSettings const > RandomGeneratorSettingsCOP;

} //basic
} //random

#endif //INCLUDED_basic_random_RandomGeneratorSettings_fwd_hh
