// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/geometry/MembraneResidueFactory.cc
///
/// @brief 		Membrane Residue Factory - Creates residues of type MEM and EMB
/// @details 	Creates a membrane residue of type MEM with AA type virtual residue. Adds the residue by
///             either a specified jump or to the end of the fold tree and then makes the new residue the root.
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_geometry_MembraneResidueFactory_hh
#define INCLUDED_core_membrane_geometry_MembraneResidueFactory_hh

// Unit headers
#include <core/membrane/geometry/MembraneResidueFactory.fwd.hh>

// Project Headers
#include <core/conformation/membrane/definitions.hh>
#include <core/conformation/membrane/Exceptions.hh>

// Package Headers
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.hh>

// C++ Headers
#include <algorithm>
#include <string>
#include <cstdlib>
#include <cmath>

namespace core {
namespace membrane {
namespace geometry {
            
    /// @brief 	Class: MembraneResidueFactory
    /// @defail Build an implicit definition of a membrane as a set of virtual
    /// 		residues representing the normal and center of the membrane and membrane embedding
    ///         within the pose framework
    class MembraneResidueFactory : public utility::pointer::ReferenceCount {
        
    public: // typedefs
        
        typedef core::chemical::ResidueTypeSetCAP ResidueTypeSet;
        typedef core::chemical::ResidueTypeCOPs ResidueTypeList;
        
    public: // methods
        
        /// @brief Constructor for membrane residue factory
        MembraneResidueFactory();
        
        /// @brief Destructor
        ~MembraneResidueFactory();
        
        /// @brief    Add a membrane definiiton residue to the pose
        /// @details  Add a residue of type MEM to the pose at a Jump nres+1 and make the
        ///           mmebrane residue the root of the fold tree
        void add_membrane_residue(
                                                          core::Vector & center,
                                                          core::Vector & normal,
                                                          core::Real & thickness,
                                                          core::pose::Pose & pose,
                                                          bool fullatom
                                                          );
 
        /// @brief      Add an embedding definition residue to the pose
        /// @details    Generate an embedding residue of tpe EMB given params and make it the root of the fold tree.
        ///             Should add by specified jump which is the chain beginning for the specific chain
        void add_embedding_residue(
                                                          core::Vector & center,
                                                          core::Vector & normal,
                                                          core::Real & depth,
                                                          core::pose::Pose & pose,
                                                          Size jump,
                                                          bool fullatom
                                                          );
        
    }; // MembraneResidueFactory
            
} // geometry
} // membrane
} // core

#endif // INCLUDED_core_membrane_geometry_MembraneResidueFactory_hh

