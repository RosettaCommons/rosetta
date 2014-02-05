// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/chemical/ResidueSupport.hh
/// @brief support functions for class residue; functions that
/// should not be included as part of the class.
/// @author Phil Bradley


#ifndef INCLUDED_core_chemical_ResidueSupport_hh
#define INCLUDED_core_chemical_ResidueSupport_hh

// Package Headers
// AUTO-REMOVED #include <core/chemical/ResidueType.hh>

// ObjexxFCL Headers
// AUTO-REMOVED #include <ObjexxFCL/FArray2D.hh>

#include <core/chemical/ResidueType.fwd.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>


namespace core {
namespace chemical {



// Find a better place to declare this function
/// @brief relies on class Graph to find all pairs shortest path information
ObjexxFCL::FArray2D_int
get_residue_path_distances( ResidueType const & res );
    
    
    
}
}

#endif
