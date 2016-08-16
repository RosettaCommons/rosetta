// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/basic/init.cc
/// @brief  Some statics/funtions helpers for Rosetta init functions.
/// @author Sergey Lyskov


#ifndef INCLUDED_basic_init_hh
#define INCLUDED_basic_init_hh

namespace basic {

/// @brief return true if core::init was already called and false otherwise
bool was_init_called();


/// @brief set global 'init_was_called' to true
void init();

} // namespace basic


#endif // INCLUDED_basic_init_hh
