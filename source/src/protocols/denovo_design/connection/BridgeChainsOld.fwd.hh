// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/denovo_design/connection/BridgeChains.fwd.hh
/// @brief  BridgeChains forward header
/// @author Tom Linsky

#if 0 // depricated

#ifndef INCLUDED_protocols_denovo_design_connection_BridgeChains_old_fwd_hh
#define INCLUDED_protocols_denovo_design_connection_BridgeChains_old_fwd_hh

// Utility headers
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace denovo_design {
namespace connection {

// Forward
class BridgeChains;

// Types
typedef  utility::pointer::shared_ptr< BridgeChains >  BridgeChainsOP;
typedef  utility::pointer::shared_ptr< BridgeChains const >  BridgeChainsCOP;

typedef  utility::pointer::weak_ptr< BridgeChains >  BridgeChainsAP;
typedef  utility::pointer::weak_ptr< BridgeChains const >  BridgeChainsCAP;

} // namespace connection
} // namespace denovo_design
} // namespace protocols

#endif

#endif // depricated
