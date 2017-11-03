// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/id/bogus.hh
/// @author Andrew Leaver-Fay
/// @details Fix the static-initialization-order fiasco w/ global variable initialization
/// Previously, the global variables in core:id were not guaranteed to be initialized in the
/// right order; by calling the function defined in this header, you can guarantee that
/// the global variables defined in the core::id namespace are initialized after the function
/// call completes.

#ifndef INCLUDED_core_id_bogus_HH
#define INCLUDED_core_id_bogus_HH

namespace core {
namespace id {

/// @brief Calling this function ensures that all of the global variables in core::id are initialize
/// before they are needed
void
initialize_core_id_globals();

} // namespace id
} // namespace core

#endif
