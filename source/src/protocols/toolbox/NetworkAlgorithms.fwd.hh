// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/NetworkAlgorithms.fwd.hh
/// @brief  FOrward header for class that analyzes properties of residue networks
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_toolbox_NetworkAlgorithms_fwd_hh
#define INCLUDED_protocols_toolbox_NetworkAlgorithms_fwd_hh

#include <utility/pointer/owning_ptr.hh>


namespace protocols{
namespace toolbox{
	class Node;
	typedef utility::pointer::shared_ptr< Node > NodeOP;
	typedef utility::pointer::shared_ptr< Node const > NodeCOP;

	class ResidueNetwork;
	typedef utility::pointer::shared_ptr< ResidueNetwork > ResidueNetworkOP;
	typedef utility::pointer::shared_ptr< ResidueNetwork const > ResidueNetworkCOP;

	class DistanceResidueNetwork;
	typedef utility::pointer::shared_ptr< DistanceResidueNetwork > DistanceResidueNetworkOP;
	typedef utility::pointer::shared_ptr< DistanceResidueNetwork const > DistanceResidueNetworkCOP;

	class CovalentResidueNetwork;
	typedef utility::pointer::shared_ptr< CovalentResidueNetwork > CovalentResidueNetworkOP;
	typedef utility::pointer::shared_ptr< CovalentResidueNetwork const > CovalentResidueNetworkCOP;

}
}

#endif
