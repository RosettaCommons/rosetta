// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/init/register4.cc
/// @brief  Declare WidgetRegistrators as static (global) variables in this .cc file
///         so that at load time, they will be initialized, and the Creator classes
///         they register will be handed to the appropriate WidgetFactory.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <protocols/init/register4.hh>

// Just the movers, as they're likely the biggest

#include <protocols/init/init.MoverCreators.ihh>
#include <protocols/init/init.MoverRegistrators.ihh>

namespace protocols {
namespace init {

void register4()
{
	// A no-op, but needed to make sure that the Creators/Registrators above get statically allocated.
	// (As they might not before a function in this compilation unit is called.
}

} //namespace init
} //namespace protocols


