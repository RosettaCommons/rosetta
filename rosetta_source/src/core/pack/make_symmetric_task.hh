// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file  core/pack/make_symmetric_task.hh
/// @brief utility functions for handling with symmetric conformations
/// @author Ingemar Andre

#ifndef INCLUDED_core_pack_make_symmetric_task_hh
#define INCLUDED_core_pack_make_symmetric_task_hh


// Project headers headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>

namespace core {
namespace pack {

void
make_symmetric_PackerTask(
  pose::Pose const & pose,
  pack::task::PackerTaskOP task
);

} // pack
} // core



#endif
