// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/ShortestPathInFoldTree.fwd.hh
/// @brief  kinematics::ShortestPathInFoldTree forward declarations header
/// @author Oliver Lange


#ifndef INCLUDED_protocols_topology_broker_claims_SequenceClaim_fwd_hh
#define INCLUDED_protocols_topology_broker_claims_SequenceClaim_fwd_hh


// Utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

#include <utility/vector1.fwd.hh>


namespace protocols {
namespace topology_broker {
namespace claims {

class SequenceClaim;

// Types
typedef  utility::pointer::shared_ptr< SequenceClaim >  SequenceClaimOP;
typedef  utility::pointer::shared_ptr< SequenceClaim const >  SequenceClaimCOP;

typedef  utility::pointer::weak_ptr< SequenceClaim >  SequenceClaimAP;
typedef  utility::pointer::weak_ptr< SequenceClaim const >  SequenceClaimCAP;

typedef utility::vector1< SequenceClaimOP > SequenceClaims;


}
}
}

#endif
