// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/operation/OptCysHG.fwd.hh
/// @brief  forward declaration for OptCysHG
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_pack_task_operation_OptCysHG_fwd_hh
#define INCLUDED_core_pack_task_operation_OptCysHG_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace task {
namespace operation {


/// @brief forward declaration for OptCysHG
class OptCysHG;


/// @brief OptCysHG owning pointer
typedef utility::pointer::shared_ptr< OptCysHG > OptCysHGOP;


/// @brief OptCysHG const owning pointer
typedef utility::pointer::shared_ptr< OptCysHG const > OptCysHGCOP;


/// @brief OptCysHG owning pointer
typedef utility::pointer::weak_ptr< OptCysHG > OptCysHGAP;


/// @brief OptCysHG const owning pointer
typedef utility::pointer::weak_ptr< OptCysHG const > OptCysHGCAP;


} // namespace operation
} // namespace task
} // namespace pack
} // namespace core


#endif /* INCLUDED_core_pack_task_operation_OptCysHG_FWD_HH */
