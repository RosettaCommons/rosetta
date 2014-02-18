// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 	 MembraneProteinFactory.cc
///
/// @brief 	 MembraneProteinFactory
/// @details The membrane protein factory creates a single pose from various membrane proteins
///			 loaded on the front end and initialized as membrane proteins. This single framework
///			 will then be passed off to the MembraneHub (which coordinates I/O) and sent back to the protocol it was
///			 called from (usually in pose loading)
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_MembraneProteinFactory_cc
#define INCLUDED_core_membrane_MembraneProteinFactory_cc

// Utility Headers
#include <core/membrane/MembraneProteinFactory.fwd.hh>

/// @brief      Membrane Protein Factory
/// @details    Initializes a pose as a membrane protein
namespace core {
namespace membrane {

    class MembraneProteinFactory;
    typedef utility::pointer::owning_ptr< MembraneProteinFactory > MembraneProteinFactoryOP;
    typedef utility::pointer::owning_ptr< MembraneProteinFactory const > MembraneProteinFactoryCOP;

} // membrane
} // core

#endif // INCLUDED_core_membrane_MembraneProteinFactory_fwd_hh

