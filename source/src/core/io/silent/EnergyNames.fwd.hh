// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/silent/EnergyNames.fwd.hh
///
/// @brief silent input file reader for mini.
/// @author James Thompson

#ifndef INCLUDED_core_io_silent_EnergyNames_fwd_hh
#define INCLUDED_core_io_silent_EnergyNames_fwd_hh

// mini headers
#include <utility/pointer/owning_ptr.hh>

//Auto Headers


namespace core {
namespace io {
namespace silent {

	class EnergyNames;

	typedef utility::pointer::shared_ptr< EnergyNames > EnergyNamesOP;
} // namespace silent
} // namespace io
} // namespace core

#endif
