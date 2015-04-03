// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief  fwd headers for ns dssp
/// @author olange: ported from original bblum-rosetta++ version $


#ifndef INCLUDED_core_scoring_dssp_PairingsList_fwd_hh
#define INCLUDED_core_scoring_dssp_PairingsList_fwd_hh


// Utility headers
//#include <utility/pointer/access_ptr.fwd.hh>

#include <utility/vector1.fwd.hh>


namespace core {
namespace scoring {
namespace dssp {

// Forward
//class BaseJumpSetup;
/// @brief list of pairings
class Pairing;
typedef utility::vector1<Pairing> PairingsList;

// new better version of the name --- try to phase out PairingsList
typedef utility::vector1<Pairing> PairingList;


// Types
// NO PairingsList is not a Ref-Counted class
//typedef  utility::pointer::owning_ptr< PairingsList >  PairingsListOP;
//typedef  utility::pointer::owning_ptr< PairingsList const >  PairingsListCOP;


} //dssp
} //scoring
} //core

#endif
