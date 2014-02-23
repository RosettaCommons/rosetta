// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/CreateMembranePoseMover.fwd.hh
///
/// @brief      Create Membrane Pose - Mover Class
/// @details    This mover creates a multi-chain membrane pose from the membrane framework
///             using a series of JD2 resource manager initialized resources. The objective
///             of the mover is to call the membrane protein factory, as for a new membrane
///             pose regardless of what has been passed on the commandline, and override that
///             pose. This should also be the top level interface to the membrane protein
///             framework.
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (2/22/14)  

#ifndef INCLUDED_protocols_membrane_CreateMembranePoseMover_fwd_hh
#define INCLUDED_protocols_membrane_CreateMembranePoseMover_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace membrane {
    
    /// @brief Class: Membrane Mover
    /// @details Helping Test Mover for initializing a membrane protein
    class CreateMembranePoseMover;
    typedef utility::pointer::owning_ptr< CreateMembranePoseMover > CreateMembranePoseMoverOP;
    typedef utility::pointer::owning_ptr< CreateMembranePoseMover const > CreateMembranePoseMoverCOP;
    
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_CreateMembranePoseMover_fwd_hh

