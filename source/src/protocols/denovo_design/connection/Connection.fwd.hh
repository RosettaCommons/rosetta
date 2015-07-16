// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/denovo_design/Connection.fwd.hh
/// @brief  Connection forward header
/// @author Tom Linsky


#ifndef INCLUDED_protocols_denovo_design_Connection_fwd_hh
#define INCLUDED_protocols_denovo_design_Connection_fwd_hh


// Utility headers
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>


namespace protocols {
namespace denovo_design {
namespace connection {

// Forward
class Connection;

// Types
typedef utility::pointer::shared_ptr< Connection > ConnectionOP;
typedef utility::pointer::shared_ptr< Connection const > ConnectionCOP;

typedef utility::pointer::weak_ptr< Connection > ConnectionAP;
typedef utility::pointer::weak_ptr< Connection const > ConnectionCAP;

// Forward
class StapleTomponents;

// Types
typedef utility::pointer::shared_ptr< StapleTomponents > StapleTomponentsOP;
typedef utility::pointer::shared_ptr< StapleTomponents const > StapleTomponentsCOP;

typedef utility::pointer::weak_ptr< StapleTomponents > StapleTomponentsAP;
typedef utility::pointer::weak_ptr< StapleTomponents const > StapleTomponentsCAP;

// Forward
class BridgeTomponents;

// Types
typedef utility::pointer::shared_ptr< BridgeTomponents > BridgeTomponentsOP;
typedef utility::pointer::shared_ptr< BridgeTomponents const > BridgeTomponentsCOP;

typedef utility::pointer::weak_ptr< BridgeTomponents > BridgeTomponentsAP;
typedef utility::pointer::weak_ptr< BridgeTomponents const > BridgeTomponentsCAP;

} // namespace connection
} // namespace denovo_design
} // namespace protocols

#endif
