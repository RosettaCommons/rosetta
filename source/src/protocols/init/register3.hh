// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/init/register3.hh
/// @brief  Initialization function to ensure all load-time factory registration occurs
///         for classes that live in the protocols library.  protocols::init::init() calls core::init::init(),
///         and devel::init() call protocols::init::init().
///         We split this off such that the compiler isn't attempting to compile *all* of the registration classes at once.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Rocco Moretti (rmorettiase@gmail.com)
//

#ifndef INCLUDED_protocols_init_register3_HH
#define INCLUDED_protocols_init_register3_HH

namespace protocols {
namespace init {

//Call the registration
void register3();

} //namespace init
} //namespace protocols

#endif //INCLUDED_protocols_init_init_HH
