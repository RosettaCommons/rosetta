// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/operation/NoRepackDisulfides.fwd.hh
/// @brief  forward declaration for NoRepackDisulfides
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_pack_task_operation_NoRepackDisulfides_fwd_hh
#define INCLUDED_core_pack_task_operation_NoRepackDisulfides_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace task {
namespace operation {


/// @brief forward declaration for NoRepackDisulfides
class NoRepackDisulfides;


/// @brief NoRepackDisulfides owning pointer
typedef utility::pointer::shared_ptr< NoRepackDisulfides > NoRepackDisulfidesOP;


/// @brief NoRepackDisulfides const owning pointer
typedef utility::pointer::shared_ptr< NoRepackDisulfides const > NoRepackDisulfidesCOP;


/// @brief NoRepackDisulfides owning pointer
typedef utility::pointer::weak_ptr< NoRepackDisulfides > NoRepackDisulfidesAP;


/// @brief NoRepackDisulfides const owning pointer
typedef utility::pointer::weak_ptr< NoRepackDisulfides const > NoRepackDisulfidesCAP;


} // namespace operation
} // namespace task
} // namespace pack
} // namespace core


#endif /* INCLUDED_core_pack_task_operation_NoRepackDisulfides_FWD_HH */
