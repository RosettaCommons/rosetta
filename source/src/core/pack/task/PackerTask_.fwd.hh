// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/PackerTask_.fwd.hh
/// @brief  Implementation class for task class to describe packer's behavior forward declaration
/// Almost all of rosetta needs to use packer tasks, but very little of rosetta needs
/// to see how it behaves internally.
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_pack_task_PackerTask__fwd_hh
#define INCLUDED_core_pack_task_PackerTask__fwd_hh

//Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace task {

class PackerTask_;

typedef utility::pointer::shared_ptr< PackerTask_ > PackerTask_OP;
typedef utility::pointer::shared_ptr< PackerTask_ const > PackerTask_COP;

} //namespace task
} //namespace pack
} //namespace core

#endif
