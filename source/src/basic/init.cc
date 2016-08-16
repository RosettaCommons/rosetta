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


namespace basic {

/// Global variable to check if core::init was ever called. Used to detect situations when user
/// forgot to call init and display helpful message instead of silently crashing.
/// We put it here instead of core/init/init.cc so code from all libraries level have access to it
static bool init_was_called = false;


bool was_init_called() { return init_was_called; }


void init()
{
	init_was_called = true;
}

} // namespace basic
