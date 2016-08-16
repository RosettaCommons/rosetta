// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/ResfileModes.fwd.hh
/// @brief  fwd header of classes for resfile options
/// @author

#ifndef INCLUDED_core_pack_task_ResfileReader_fwd_hh
#define INCLUDED_core_pack_task_ResfileReader_fwd_hh

//Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace task {

class ResfileCommand;

typedef utility::pointer::shared_ptr< ResfileCommand > ResfileCommandOP;
typedef utility::pointer::shared_ptr< ResfileCommand const > ResfileCommandCOP;

class ResfileContents;

typedef utility::pointer::shared_ptr< ResfileContents > ResfileContentsOP;
typedef utility::pointer::shared_ptr< ResfileContents const > ResfileContentsCOP;

} //namespace task
} //namespace pack
} //namespace core


#endif
