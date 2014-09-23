// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_toolbox_KCluster_fwd_hh
#define INCLUDED_protocols_toolbox_KCluster_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace toolbox {

class KClusterElement;
typedef utility::pointer::shared_ptr< KClusterElement > KClusterElementOP;
typedef utility::pointer::shared_ptr< KClusterElement const > KClusterElementCOP;

class KClusterData;
typedef utility::pointer::shared_ptr< KClusterData > KClusterDataOP;
typedef utility::pointer::shared_ptr< KClusterData const > KClusterDataCOP;

class KCluster;
typedef utility::pointer::shared_ptr< KCluster > KClusterOP;
typedef utility::pointer::shared_ptr< KCluster const > KClusterCOP;

class KMedoid;
typedef utility::pointer::shared_ptr< KMedoid > KMedoidOP;
typedef utility::pointer::shared_ptr< KMedoid const > KMedoidCOP;

class GreedyKCenter;
typedef utility::pointer::shared_ptr< GreedyKCenter > GreedyKCenterOP;
typedef utility::pointer::shared_ptr< GreedyKCenter const > GreedyKCenterCOP;

}
}

#endif
