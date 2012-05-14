// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   STMStoredTask.fwd.hh
/// @brief
/// @author Neil King (neilking@uw.edu)


#ifndef INCLUDED_devel_matdes_STMStoredTask_fwd_hh
#define INCLUDED_devel_matdes_STMStoredTask_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace devel {
namespace matdes {

class STMStoredTask;
typedef utility::pointer::owning_ptr< STMStoredTask > STMStoredTaskOP;
typedef utility::pointer::owning_ptr< STMStoredTask const > STMStoredTaskCOP;

} // matdes
} // devel

#endif
