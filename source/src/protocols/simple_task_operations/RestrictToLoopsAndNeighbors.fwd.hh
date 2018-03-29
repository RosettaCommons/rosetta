// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_task_operations/RestrictToLoopsAndNeighbors.fwd.hh
/// @brief  Forward declaration of the RestrictToLoopsAndNeighbors class
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)

#ifndef INCLUDED_protocols_task_operations_RestrictToLoopsAndNeighbors_FWD_HH
#define INCLUDED_protocols_task_operations_RestrictToLoopsAndNeighbors_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace simple_task_operations {

class RestrictToLoopsAndNeighbors;

typedef utility::pointer::shared_ptr< RestrictToLoopsAndNeighbors > RestrictToLoopsAndNeighborsOP;

} //namespace simple_task_operations
} //namespace protocols

#endif // INCLUDED_protocols_task_operations_RestrictToLoopsAndNeighbors_FWD_HH
