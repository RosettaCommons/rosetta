// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/mover/MembraneUnitTestMover.fwd.hh
///
/// @brief      Top-Level Unit Test for the Membrane Protein Factory (Mover Class)
/// @details    The purpose of this application is to test the membrane protein factory
///             initialization code including all external dependencies which cannot be tested
///             in JD2, Resource Manager, and the Pose cache. This can also serve as the integration
///             test for memrane protein initialization.
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (2/22/14)

#ifndef INCLUDED_protocols_membrane_MembraneUnitTestMover_fwd_hh
#define INCLUDED_protocols_membrane_MembraneUnitTestMover_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace membrane {
    
    /// @brief Class: Membrane Unit Testing Mover
    /// @details Helping Test Mover for initializing a membrane protein
    class MembraneUnitTestMover;
    typedef utility::pointer::owning_ptr< MembraneUnitTestMover > MembraneUnitTestMoverOP;
    typedef utility::pointer::owning_ptr< MembraneUnitTestMover const > MembraneUnitTestMoverCOP;
    
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_MembraneUnitTestMover_fwd_hh

