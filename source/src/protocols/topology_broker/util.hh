// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_topology_broker_util_hh
#define INCLUDED_protocols_topology_broker_util_hh

// Package Headers
#include <protocols/topology_broker/TopologyBroker.fwd.hh>

namespace protocols {
namespace topology_broker {

/// @brief read broker::setup flag and add all claims to Broker
void add_cmdline_claims(TopologyBroker&, bool do_I_need_frags=true);

}
}

#endif
