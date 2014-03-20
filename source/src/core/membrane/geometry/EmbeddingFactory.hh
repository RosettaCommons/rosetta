// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/geometry/EmbeddingFactory.hh
///
/// @brief 		Factory Method for creating a membrane embedding from user options and tags
/// @details 	Creates reference center, normal, depth/spanning definition and topology
///				definition from the user's specifications
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_geometry_EmbeddingFactory_hh
#define INCLUDED_core_membrane_geometry_EmbeddingFactory_hh

// Unit Headers
#include <core/membrane/geometry/EmbeddingFactory.fwd.hh>

// Project Headers
#include <core/membrane/properties/SpanningTopology.hh>
#include <core/membrane/util/definitions.hh>
#include <core/membrane/util/Exceptions.hh>

#include <core/membrane/geometry/util.hh>

// Package Headers
#include <core/pose/Pose.hh>

#include <core/conformation/Residue.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray4D.hh>

#include <basic/Tracer.hh>
#include <basic/prof.hh>

// Numeric Headers
#include <numeric/xyzVector.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector1.hh>

// C++ Headers
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <string>

using namespace core::membrane::util;

namespace core {
namespace membrane {
namespace geometry {
    
/// @brief Create an embedding definition from opts
/// @details Generate the correct protein emebdding using a user-specified
/// method
class EmbeddingFactory : public utility::pointer::ReferenceCount {

public: // methods
    
    /// @brief      Default Constructor
    EmbeddingFactory(
                     core::pose::PoseOP pose,
                     EmbedConfigInfoOP config,
                     core::membrane::properties::SpanningTopologyOP topology
                     );
    
    /// @brief      Default Destructor
    ~EmbeddingFactory();

    /// @brief      Create Embedding Method
    /// @details    Creates a membrane protein embedding definition from embedidng config
    ///             object loaded at resource definiiton (cmd). The embedding definition is a
    ///             final definition of the starting pose's definition in the membrane
    ///
    /// @param      description
    ///                 resource description
    /// @param      fullatom
    ///                 residue typeset
    void create_and_add_embedding( bool fullatom );

    /// @brief      Create Embedding Method
    /// @details    Creates a membrane protein embedding definition from embedidng config
    ///             object loaded at resource definiiton (cmd). The embedding definition is a
    ///             final definition of the starting pose's definition in the membrane
    ///
    /// @param      description
    ///                 resource description
    /// @param      fullatom
    ///                 residue typeset
    void create_and_add_embedding( bool fullatom, core::Size jump );
    
    /// @brief      Embed from Defaults
	/// @details    Define the membrane embedding from a pre-determined set of parameters
	///             Normal (0, 0, 1), Center (0, 0, 0)
	void embed_from_pdb( core::membrane::util::EmbedConfigInfoOP embedding );
    
	/// @brief     Embed from membrane spanning topology
	/// @details   Determines embedding by finding the average point between Ca inner and Ca outer (distance
	///            between the intracellualr and extracellular region).
	void embed_from_topology( core::membrane::util::EmbedConfigInfoOP embedding );
    
	/// @brief      Embedding Search and Score Method
	/// @details    Starting from embedding defined from topology (see method above), score each conformation, perturb,
	///             score, and aver n trials and MC accept lowest. Takes the embed params resource loaded via resource
    ///             loader (RM).
	void embed_from_search_and_score( core::membrane::util::EmbedConfigInfoOP embedding );
    
	/// @brief		Compute Fullatom Projections for Fa Embedding
	///	@details	Computes Projections, depth, and stores coordinates. Computes fullatom
	///				center and corresponding normal from these computations
	void embed_for_FA( core::membrane::util::EmbedConfigInfoOP embedding );
	
    /// @brief   Non Membrane Embedding Method
    /// @details The embedding config info object will contain the data required to initialize this
    ///          information:
    ///             center = center of mass of the non membrane chain
    ///             normal = vector sum of residue-com
    ///             detph = 30A + some user specified factory
    void
    embed_for_non_membrane( core::membrane::util::EmbedConfigInfoOP embedding );
	

private: // data
    
    // embedding config info
    core::membrane::util::EmbedConfigInfoOP embed_info_;
    
    // topology info
    core::membrane::properties::SpanningTopologyOP topology_;
    
    // pose
    core::pose::PoseOP pose_;

}; // class EmbeddingFactory
    
} // core
} // membrane
} // geometry

#endif // INCLUDED_core_membrane_geometry_EmbeddingFactory_hh

