// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 	 MembraneProteinFactory.hh
///
/// @brief 	 MembraneProteinFactory
/// @details The membrane protein factory creates a single pose from various membrane proteins
///			 loaded on the front end and initialized as membrane proteins. This single framework
///			 will then be passed off to the MembraneHub (which coordinates I/O) and sent back to the protocol it was
///			 called from (usually in pose loading)
///
/// @author  Rebecca Alford (rfalford12@gmail.com)


/// @brief Class: Membrane Protein Factory
/// @details Construct a membrane protein in Rosetta
/// @author Rebecca Alford (rfalford12@gmail.com)
///
/// The membrane protein factory constructs a membrane protein in Rosetta through the following steps:
///     (1) Construct fold tree with pre-defined membrane fold tree topology
///     (2) Add a MembraneConformation which maintains info about the membrane protein
///     (3) Add a virtual residue to the root of the foldtree to define the membrane
///     (4) Add a new virtual residue representing the embedding of each protein chain in the membrane
///
/// The membrane protein factory provides three constructors - single chain construction, multi chain construction via the resource
/// manager and multi chain construction via the options system (existing system). The new membrane protein framework
/// is not backwards compatible with the previous RosettaMembrane code (Yarov-Yaravoy et al. 2006) and uses a different set of options
///
/// Maintains the following invariants:
///   (1) For an n chain pose, will return a pose with n+1 chains (n chains plus a membrane chain)
///   (2) Initializes a correct membrane foldtree topology (initial)
///   (3) Membrane chain is always the n+1 chain
///   (4) Membrane fold tree is a valid fold tree
///
/// Notes:
///     - This class is not a mover. It will not move the backbone whatsoever (see MembraneRigidInitialMover)
///

#ifndef INCLUDED_core_membrane_MembraneProteinFactory_hh
#define INCLUDED_core_membrane_MembraneProteinFactory_hh

// Unit Headers
#include <core/membrane/MembraneProteinFactory.fwd.hh>

// Project Headers
#include <core/membrane/properties/SpanningTopology.hh>
#include <core/membrane/properties/LipidAccInfo.hh>
#include <core/membrane/util/definitions.hh>
#include <core/membrane/util/Exceptions.hh>

#include <core/membrane/geometry/MembraneResidueFactory.hh>
#include <core/membrane/geometry/EmbeddingFactory.hh>
#include <core/membrane/geometry/util.hh>

#include <core/membrane/MembraneInfo.hh>
#include <core/conformation/Conformation.hh>

// Package Headers
#include <core/pose/Pose.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>
#include <algorithm>

using namespace core::membrane::properties;
using namespace core::membrane;
using namespace core::pose;

namespace core {
namespace membrane {

    /// @brief      Membrane Protein Factory
    /// @details    Initializes a pose as a membrane protein in Rosetta
    class MembraneProteinFactory : public utility::pointer::ReferenceCount {
        
    public: // constructors
        
        /// @brief   Resource Manager Constructor (No Options Specified)
        /// @details Construct a Membrane Protein Factory loading all required resources
        ///          from the resource manager using the default option settings
        ///
        /// @throws EXCN_Illegal_Arguments
        ///             if file is "" - cannot construct a membrane protein without chain refs for the resource manager
        ///
        /// @return MembraneProteinFactory
        MembraneProteinFactory( std::string membrane_chains );
        
        /// @brief   Resource Manager Constructor (Non-Default)
        /// @details Construct a Membrane Protein Factory loading all required resources
        ///          from the resource manager
        ///
        /// @param  include_lips
        ///             load and include lipid accessibility data in scoring
        /// @param  membrane_chains
        ///             text file storing references to membrane chains in the membrane protein
        /// @param  fullatom
        ///             specified fullatom residue typeset
        ///
        /// @return MembraneProteinFactory
        MembraneProteinFactory(
                               bool include_lips,
                               std::string membrane_chains,
                               bool fullatom
                               );
        
        /// @brief   Interactive Constructor (With Default MP Options)
        /// @details Construct a Membrane Protein Factory loading all required resources
        ///          from user specified inputs
        ///
        /// @param chians
        ///         list of pose chains to compose into a single pose
        /// @param topologies
        ///         list of spannign topology objects for each chain
        /// @param embeddings
        ///         list of embedding data objects for each chain
        /// @param lipid_acc
        ///         list of lipid accessibility data objects for each chain
        ///
        /// @return MembraneProteinFactory
        MembraneProteinFactory(
                               utility::vector1< PoseOP > chains,
                               utility::vector1< SpanningTopologyOP > topologies,
                               utility::vector1< core::membrane::util::EmbedConfigInfoOP > embeddings,
                               utility::vector1< LipidAccInfoOP > lipid_acc
                               );
        
        /// @brief   Interactive Constructor (with non default options)
        /// @details Construct a Membrane Protein Factory loading all required resources
        ///          from user specified inputs
        ///
        /// @param  include_lips
        ///             load and include lipid accessibility data in scoring
        /// @param  fullatom
        ///             specified fullatom residue typeset
        /// @param  chians
        ///             list of pose chains to compose into a single pose
        /// @param  topologies
        ///             list of spannign topology objects for each chain
        /// @param  embeddings
        ///             list of embedding data objects for each chain
        /// @param  lipid_acc
        ///             list of lipid accessibility data objects for each chain
        ///
        /// @return MembraneProteinFactory
        MembraneProteinFactory(
                               bool include_lips,
                               bool fullatom,
                               utility::vector1< PoseOP > chains,
                               utility::vector1< SpanningTopologyOP > topologies,
                               utility::vector1< core::membrane::util::EmbedConfigInfoOP > embeddings,
                               utility::vector1< LipidAccInfoOP > lipid_acc
                               );
        
        /// @brief    Default Destructor
        /// @details
        ///
        /// @note
        ~MembraneProteinFactory();
        
        /// @brief 	 Create Membrane Protein
        /// @details Create a membrne proteins from a series of loaded membrane proteins
        ///
        /// @return  Pose (as starting structure)
        core::pose::PoseOP create_membrane_pose();
        
    private: // default constructors

        /// @brief  Default Constructor (Private
        /// @details Private Default constructor (DO NOT USE)
        MembraneProteinFactory();
        
    private: // methods
        
        /// @brief      Initialize Chains
        /// @details    Initialize membrane chains from initialized prefix file provided in the constructor
        ///
        /// @throws     Argument exception if chain list not specified (also this is well docuemnted, no excuse)
        void initialize_chains();
        
        /// @brief Build Pose
        /// @details Create pose containing membrane/embedding residues from multi-
        /// chain input.
        /// @throws EXCN_Resource_Manager, EXCN_Membrane_Bounds
        void build_pose( core::pose::PoseOP pose );
        
        /// @brief Initialize Spanning Topology
        /// @details Initialize spanning topology in the final pose
       void initialize_topology( core::pose::PoseOP pose );
        
        /// @brief Initialize Lipds Exposure Data
        /// @details Initialize lipid exposure data in the final pose
       void initialize_lips_exp( core::pose::PoseOP pose );
        
        /// @brief Show Membrane Pose Info
        /// @details Show all ifnormation relevant to the new membrane pose
        void show( std::ostream & out );
        
        /// @brief      Initialize Resources from the Resource Manager
        /// @details    Load required resources for initializing a membrane protein
        ///
        /// @note Precondition: Initialized prefix list
        /// @throws EXCN_Resource_Manager (missing reuqired resource)
        void initialize_resources();
                
    private: // data
        
        // Map a list of chains to their corresponding base path
        std::map< core::Size, std::string > chains_map_;
        
        // Option Settings
        std::string prefix_file_;
        bool fullatom_;
        
        // Include Lips data
        bool include_lips_;
        
        // List of 'by chain' resources
        utility::vector1< core::pose::PoseOP > chains_;
        utility::vector1< core::membrane::properties::SpanningTopologyOP > topologies_;
        utility::vector1< core::membrane::util::EmbedConfigInfoOP > embeddings_;
        utility::vector1< core::membrane::properties::LipidAccInfoOP > lipid_acc_;
        
        // Factory instances
        core::membrane::geometry::MembraneResidueFactory mrf_;
    
    }; // class MembraneProteinFactory

} // membrane
} // core

#endif // INCLUDED_core_membrane_MembraneProteinFactory_hh



