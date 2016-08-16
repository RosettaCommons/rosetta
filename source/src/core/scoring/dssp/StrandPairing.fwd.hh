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


#ifndef INCLUDED_core_scoring_dssp_StrandPairing_fwd_hh
#define INCLUDED_core_scoring_dssp_StrandPairing_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>


namespace core {
namespace scoring {
namespace dssp {

// Forward
class StrandPairing;
class StrandPairingSet;

// Types
typedef  utility::pointer::shared_ptr< StrandPairing >  StrandPairingOP;
typedef  utility::pointer::shared_ptr< StrandPairing const >  StrandPairingCOP;

typedef  utility::pointer::shared_ptr< StrandPairingSet > StrandPairingSetOP;
typedef  utility::pointer::shared_ptr< StrandPairingSet const >  StrandPairingSetCOP;

} // namespace dssp
} // namespace scoring
} // namespace core

#endif
