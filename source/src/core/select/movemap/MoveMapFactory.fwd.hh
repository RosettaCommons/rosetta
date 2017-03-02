// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/movemap/MoveMapFactory.fwd.hh
/// @brief  Forward declaration of MoveMapFactory class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_select_movemap_MoveMapFactory_FWD_HH
#define INCLUDED_core_select_movemap_MoveMapFactory_FWD_HH

// Utility Headers
#include <utility/pointer/owning_ptr.hh>


namespace core {
namespace select {
namespace movemap {

class MoveMapFactory;

typedef utility::pointer::shared_ptr< MoveMapFactory > MoveMapFactoryOP;
typedef utility::pointer::shared_ptr< MoveMapFactory const > MoveMapFactoryCOP;

} //namespace movemap
} //namespace select
} //namespace core


#endif
