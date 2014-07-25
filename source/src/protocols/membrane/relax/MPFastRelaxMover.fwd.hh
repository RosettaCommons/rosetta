// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/relax/MPFastRelaxMover.fwd.hh
///
/// @brief      Basic FastRelax Protocol for Membrane Proteins
/// @details    Refinement and minimization of membrane protein structures using an adapted
///				version of the FastRelax protocol. Uses the membrane framework and adapted
///				minimization settings. 			
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (7/20/14)

#ifndef INCLUDED_protocols_membrane_relax_MPFastRelaxMover_fwd_hh
#define INCLUDED_protocols_membrane_relax_MPFastRelaxMover_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh> 

namespace protocols {
namespace membrane {
namespace relax {

class MPFastRleaxMover;
typedef utility::pointer::owning_ptr< MPFastRleaxMover > MPFastRleaxMoverOP; 
typedef utility::pointer::owning_ptr< MPFastRleaxMover const > MPFastRleaxMoverCOP; 

} // relax
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_relax_MPFastRelaxMover_fwd_hh

