// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file	    protocols/membrane/MembranePositionFromTopologyMover.fwd.hh
///
/// @brief      Computes and sets the initial position of the membrane
/// @details	Computes and sets the initial position of the membrane from
///				sequence or structure (can be specified by the user at construction
///				or as a setup cmd flag).
///				CAUTION: ONLY FOR FLEXIBLE MEMBRANE AND FIXED PROTEIN!!!
///
///				NOTE: Requires a membrane pose!
///				NOTE: sequence not yet implemented
///				Last Modified: 6/21/14
///
/// @author		Rebecca Alford (rflaford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_MembranePositionFromTopologyMover_fwd_hh
#define INCLUDED_protocols_membrane_MembranePositionFromTopologyMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace membrane {
		
class MembranePositionFromTopologyMover;
typedef utility::pointer::shared_ptr< MembranePositionFromTopologyMover > MembranePositionFromTopologyMoverOP;
typedef utility::pointer::shared_ptr< MembranePositionFromTopologyMover const > MembranePositionFromTopologyMoverCOP;
		
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_MembranePositionFromTopologyMover_fwd_hh
