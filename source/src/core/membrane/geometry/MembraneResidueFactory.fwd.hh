// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/geometry/MembraneResidueFactory.hh
///
/// @brief 		Membrane Residue Factory - Creates Residue of Type MEM
/// @details 	Creates a membrane residue of type MEM with AA type virtual residue.
///             Uses the residue params framework for specifying parameters of a membrane residue.
///             Contains data for the memrbane normal, center and thicnkess.
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_geometry_MembraneResidueFactory_fwd_hh
#define INCLUDED_core_membrane_geometry_MembraneResidueFactory_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace membrane {
namespace geometry {
    
    /// @brief 	Class: MembraneResidueFactory
    /// @defail Build an implicit definition of a membrane as a set of virtual
    /// 		atoms representing the normal and center of the membrane within the
    /// 		pose framework.
    class MembraneResidueFactory;
    typedef utility::pointer::owning_ptr< MembraneResidueFactory > MembraneResidueFactoryOP;
    typedef utility::pointer::owning_ptr< MembraneResidueFactory const > MembraneResidueFactoryCOP;

    
} // geometry
} // membrane
} // core

#endif // INCLUDED_core_membrane_geometry_MembraneResidueFactory_fwd_hh

