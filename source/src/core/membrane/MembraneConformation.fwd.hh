// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 	 core/membrane/MembraneConformation.fwd.hh
///
/// @brief 	 Membrane Conformation
/// @details The Membrane Conformation is responsible for:
///             - maintaining a correct membrane foldtree
///             - maintaining references to the membrane and embedding residues
///             - providing access to membrane related data
///
/// @note    Last Modified 1/9/14
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_MembraneConformation_fwd_hh
#define INCLUDED_core_membrane_MembraneConformation_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh> 

namespace core {
namespace membrane {
    
    /// @brief Class: Membrane Conformation
    /// @details Handles memrbane conformation, foldtree, and maintains memrbane info
    class MembraneConformation;
    typedef utility::pointer::owning_ptr< MembraneConformation > MembraneConformationOP;
    typedef utility::pointer::owning_ptr< MembraneConformation const > MembraneConformationCOP;
    
} // membrane
} // core

#endif // INCLUDED_core_membrane_MembraneConformation_fwd_hh

