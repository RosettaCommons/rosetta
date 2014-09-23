// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/cluster/APCluster.fwd.hh
///
/// @brief
/// @author Ian W. Davis


#ifndef INCLUDED_protocols_cluster_APCluster_fwd_hh
#define INCLUDED_protocols_cluster_APCluster_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace cluster {


class APCluster; // fwd declaration
typedef utility::pointer::shared_ptr< APCluster > APClusterOP;
typedef utility::pointer::shared_ptr< APCluster const > APClusterCOP;


} // namespace cluster
} // namespace protocols

#endif // INCLUDED_protocols_cluster_APCluster_FWD_HH
