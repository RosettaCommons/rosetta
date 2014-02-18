// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/geometry/MembraneResidueFactory.fwd.hh
///
/// @brief 		Methods to check that poses are arranged within reasonable bounds of a membrane
/// @detailed 	Checks refs between centers, normal, depth, thickness, and spanning definitions
///				Provides exception handling
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_geometry_MembraneResidueFactory_fwd_hh
#define INCLUDED_core_membrane_geometry_MembraneResidueFactory_fwd_hh

// Utility Headers
#include <utility/pointer/owining_ptr.hh>

namespace core {
namespace membrane {
namespace geometry {
    
    /// @brief   Class: MembraneBoundsChecking
    /// @detail  Reasonable bounds checking and exception handling for membrane
    ///          framework definition
    class MembraneResidueFactory;
    typedef utility::pointer::owining_ptr< MembraneResidueFactory > MembraneResidueFactoryOP;
    typedef utility::pointer::owning_ptr< MembraneResidueFactory const > MembraneResidueFactoryCOP;
    
    
} // geometry
} // membrane
} // core

#endif // INCLUDED_core_membrane_geometry_MembraneResidueFactory_hh


