// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody_design/CDRClusterFeatures.fwd.hh
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_clusters_CDRClusterFeatures_fwd_hh
#define INCLUDED_protocols_antibody_clusters_CDRClusterFeatures_fwd_hh

#include  <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace antibody {
namespace clusters {

class CDRClusterFeatures;
typedef utility::pointer::shared_ptr< CDRClusterFeatures > CDRClusterFeaturesOP;
typedef utility::pointer::shared_ptr< CDRClusterFeatures const > CDRClusterFeaturesCOP;

}
}
}

#endif //INCLUDED_protocols_antibody_clusters_CDRClusterFeatures_fwd_hh
