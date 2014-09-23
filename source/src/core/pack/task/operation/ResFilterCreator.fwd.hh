// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief  Base class for ResFilter factory-registration and creation classes.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author ashworth


#ifndef INCLUDED_core_pack_task_operation_ResFilterCreator_fwd_hh
#define INCLUDED_core_pack_task_operation_ResFilterCreator_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace task {
namespace operation {

class ResFilterCreator;

typedef utility::pointer::shared_ptr< ResFilterCreator > ResFilterCreatorOP;
typedef utility::pointer::shared_ptr< ResFilterCreator const > ResFilterCreatorCOP;

}
}
}
}

#endif
