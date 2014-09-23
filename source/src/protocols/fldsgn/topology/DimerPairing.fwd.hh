// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


#ifndef INCLUDED_protocols_fldsgn_topology_DimerPairing_fwd_hh
#define INCLUDED_protocols_fldsgn_topology_DimerPairing_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace fldsgn {
namespace topology {

class DimerPairing;
typedef utility::pointer::shared_ptr< DimerPairing > DimerPairingOP;
typedef utility::pointer::shared_ptr< DimerPairing const > DimerPairingCOP;

class DimerPairings;

} // ns topology
} // ns fldsgn
} // ns protocols

#endif
