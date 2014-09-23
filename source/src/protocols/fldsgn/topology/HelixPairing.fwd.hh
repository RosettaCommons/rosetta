 // -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// Copyright in the Rosetta software belongs to the developers and their institutions.
// For more information, see www.rosettacommons.org.

/// @file   ./src/protocols/fldsgn/topology/HelixPairing.fwd.hh
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@u.washington.edu )

#ifndef INCLUDED_protocols_fldsgn_topology_HelixPairing_fwd_hh
#define INCLUDED_protocols_fldsgn_topology_HelixPairing_fwd_hh

// utitlity headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace fldsgn {
namespace topology {

	class HelixPairing;
	class HelixPairingSet;

	typedef utility::pointer::shared_ptr< HelixPairing > HelixPairingOP;
	typedef utility::pointer::shared_ptr< HelixPairingSet > HelixPairingSetOP;
	typedef utility::pointer::shared_ptr< HelixPairing const > HelixPairingCOP;
	typedef utility::pointer::shared_ptr< HelixPairingSet const > HelixPairingSetCOP;
	typedef utility::vector1< HelixPairingOP > HelixPairings;

} // namespace topology
} // namespace fldsgn
} // namespace protocols

#endif
