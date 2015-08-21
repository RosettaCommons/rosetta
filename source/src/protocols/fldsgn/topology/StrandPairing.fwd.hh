// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   ./src/protocols/fldsgn/topology/StrandPairing.fwd.hh
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@u.washington.edu )

#ifndef INCLUDED_protocols_fldsgn_topology_StrandPairing_fwd_hh
#define INCLUDED_protocols_fldsgn_topology_StrandPairing_fwd_hh

// utitlity headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace fldsgn {
namespace topology {

class StrandPairing;
class StrandPairingSet;
class Residues;

typedef utility::pointer::shared_ptr< StrandPairing > StrandPairingOP;
typedef utility::pointer::shared_ptr< StrandPairingSet > StrandPairingSetOP;
typedef utility::pointer::shared_ptr< StrandPairing const > StrandPairingCOP;
typedef utility::pointer::shared_ptr< StrandPairingSet const > StrandPairingSetCOP;
typedef utility::vector1< StrandPairingOP > StrandPairings;

// typedef utility::pointer::owning_ptr< Residues > ResiduesOP;
//typedef utility::pointer::owning_ptr< StrandPairing const > StrandPairingCOP;
//typedef utility::pointer::owning_ptr< StrandPairingSet const > StrandPairingSetCOP;

} // namespace topology
} // namespace fldsgn
} // namespace protocols

#endif
