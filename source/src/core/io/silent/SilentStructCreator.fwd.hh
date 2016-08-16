// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/silent/SilentStructCreator.hh
/// @brief  Base class for SilentStructCreators for the SilentStruct load-time factory registration scheme
/// @author James Thompson

#ifndef INCLUDED_core_io_silent_SilentStructCreator_fwd_hh
#define INCLUDED_core_io_silent_SilentStructCreator_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace io {
namespace silent {

/// @brief Abstract base class for a SilentStruct factory; the Creator class is responsible for
/// creating a particular SilentStruct class.
class SilentStructCreator;

typedef utility::pointer::shared_ptr< SilentStructCreator > SilentStructCreatorOP;
typedef utility::pointer::shared_ptr< SilentStructCreator const > SilentStructCreatorCOP;

} //namespace silent
} //namespace io
} //namespace core

#endif
