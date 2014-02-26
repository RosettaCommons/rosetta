// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/geometry/EmbeddingFactory.cc
///
/// @brief 		Factory Method for creating a membrane embedding from user options and tags
/// @details 	Creates reference center, normal, depth/spanning definition and topology
///				definition from the user's specifications
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_geometry_EmbeddingFactory_cc
#define INCLUDED_core_membrane_geometry_EmbeddingFactory_cc

// Unit headers
#include <core/membrane/geometry/EmbeddingFactory.hh>

// Project Headers
#include <core/membrane/properties/SpanningTopology.hh>
#include <core/membrane/geometry/util.hh>
#include <core/membrane/util/definitions.hh>
#include <core/membrane/util/Exceptions.hh>

#include <core/membrane/geometry/MembraneResidueFactory.hh>
#include <core/membrane/geometry/util.hh>

#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/util.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>

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

using namespace core;
using namespace core::membrane::properties;

using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.membrane.geometry.EmbeddingFactory");

namespace core {
namespace membrane {
namespace geometry {
    
    /**
     
     Ways to fix the logic in the file format
     
     POSE <PDB_ID>
     <type tag> /// <== MEMBRANE or NON_MEMBRANE
     center x y z tag
     normal x y z tag
     depth x y z tag
     
     
     **/
            

    /// @brief  Default Constructor
    EmbeddingFactory::EmbeddingFactory(
                                       core::pose::PoseOP pose,
                                       EmbedConfigInfoOP config,
                                       core::membrane::properties::SpanningTopologyOP topology
                                       ) :
        utility::pointer::ReferenceCount(),
        embed_info_(NULL),
        topology_(NULL),
        pose_(NULL)
    {
        embed_info_ = config;
        topology_ = topology;
        pose_ = pose;
    }
    
    /// @brief      Default Destructor
    EmbeddingFactory::~EmbeddingFactory() {}
    
    /// @brief      Create Embedding Method
    /// @details    Creates a membrane protein embedding definition from embedidng config
    ///             object loaded at resource definiiton (cmd). The embedding definition is a
    ///             final definition of the starting pose's definition in the membrane
    ///
    /// @param      description
    ///                 resource description
    /// @param      fullatom
    ///                 residue typeset
    void EmbeddingFactory::create_and_add_embedding( bool fullatom ) {
        
        // Grab the jump from chain begin
        core::conformation::Conformation & conf = pose_->conformation();
        core::Size jump = conf.chain_begin(1); // assuming single chain for this (initialize by chain is the rule)
        
        // Create and add embedding
        create_and_add_embedding( fullatom, jump );
    }
    
    /// @brief      Create Embedding Method
    /// @details    Creates a membrane protein embedding definition from embedidng config
    ///             object loaded at resource definiiton (cmd). The embedding definition is a
    ///             final definition of the starting pose's definition in the membrane
    ///
    /// @param      description
    ///                 resource description
    /// @param      fullatom
    ///                 residue typeset
    void
    EmbeddingFactory::create_and_add_embedding( bool fullatom, core::Size jump ) {
        
        using namespace core::membrane::geometry;
        using namespace core::membrane::util;
        
        // throw exception for unreasonable tag???
        
        // Create an instance of the Membrane Residue Factory
        MembraneResidueFactory mrf;
        
        // Declare accept parameters
        EmbedConfigInfoOP accept = new core::membrane::util::EmbedConfigInfo;
        accept->center = embed_info_->center;
        accept->normal = embed_info_->normal;
        accept->depth = embed_info_->depth;
        
        // Grab Method Indicator Tags froM Embed Config file
        std::string normal_tag = embed_info_->normal_tag;
        std::string center_tag = embed_info_->center_tag;

        // Compute new normal and center based on indicator tags
        if ( normal_tag.compare("override_pdb") == 0 || center_tag.compare("override_pdb") == 0 ) {
        
            EmbedConfigInfoOP temp = new EmbedConfigInfo;
            embed_from_pdb( temp );
        
            if ( normal_tag.compare("override_pdb") == 0 ) {
                accept->normal = temp->normal;
            } else if ( center_tag.compare("override_pdb") == 0 ) {
                accept->center = temp->center;
            }
        }
        
        if ( normal_tag.compare("topology") == 0 || center_tag.compare("topology") == 0 ) {
            
            EmbedConfigInfoOP temp = new EmbedConfigInfo;
            embed_from_topology( temp );
            
            if ( normal_tag.compare("topology") == 0 ) {
                accept->normal = temp->normal;
            } else if ( center_tag.compare("topology") ) {
                accept->center = temp->center;
            }
        }
        
        if ( normal_tag.compare("search") == 0 || center_tag.compare("search") == 0 ) {
            
            EmbedConfigInfoOP temp = new EmbedConfigInfo;
            embed_from_search_and_score( temp );
            
            if ( normal_tag.compare("search") == 0 ) {
                accept->normal = temp->normal;
            } else if ( center_tag.compare("search") ) {
                accept->center = temp->center;
            }
        }
        
        TR << "Embedding factory accepting the following coordinates: " << std::endl;
        TR << "Accept center: " << accept->center.x() << " " << accept->center.y() << " " << accept->center.z() << std::endl;
        TR << "Accept normal: " << accept->normal.x() << " " << accept->normal.y() << " " << accept->normal.z() << std::endl;
        TR << "Accept depth: " << accept->depth << std::endl;
        
        // Cnstruct an embedding definition from config
        mrf.add_embedding_residue( accept->center, accept->normal, accept->depth, *pose_, jump, fullatom);

    }
    
    /////////////////// Embedding Helper Methods for Factory ///////////
    
    /// @brief      Embed from Defaults
	/// @details    Define the membrane embedding from a pre-determined set of parameters
	///             Normal (0, 0, 1), Center (0, 0, 0)
	///
	/// @param      EmbedConfigInfoOP
    ///                 ref to an embed config info obj passed by the embedding factory
	void
    EmbeddingFactory::embed_from_pdb( core::membrane::util::EmbedConfigInfoOP embedding ) {
        
		TR << "Embed from PDB: Calculating embedding from pdb" << std::endl;
        
		// Initialize depth, normal and center for a 2x30 membrane
		embedding->center.assign(0, 0, 0);
		embedding->normal.assign(0, 0, 1);
        
        TR << "Embed from PDB: Calculation complete!" << std::endl;
        
		// Done!
		return;
	}
    
	/// @brief     Embed from membrane spanning topology
	/// @details   Determines embedding by finding the average point between Ca inner and Ca outer (distance
	///            between the intracellualr and extracellular region).
    ///
	/// @param     EmbedConfigInfoOP
    ///             embedding config info to modify
	void
    EmbeddingFactory::embed_from_topology( core::membrane::util::EmbedConfigInfoOP embedding  ) {
        
		TR << "Embed from Topology: Calculating embedding from topology" << std::endl;
        
		// Initialize inside and outside
		Vector inside(0);
		Vector outside(0);
        
		// For each helix in pose, get spanning region
		for ( Size i = 1; i <= topology_->total_tmhelix(); ++i ) {
			Vector const & start = pose_->residue( topology_->span()(i, 1) ).atom( 2 ).xyz();
			Vector const & end = pose_->residue( topology_->span()(i, 2) ).atom( 2 ).xyz();
            
			// Add odd helices from outside in
			if ( topology_->helix_id()(i) % 2 == 0 ) {
				inside += start;
				outside += end;
			} else {
				outside += start;
				inside += end;
			}
		}
        
		// Calculate new normal and center from spanning
		embedding->normal = outside - inside;
		embedding->normal.normalize();
		embedding->center = 0.5*(outside+inside)/topology_->tmh_inserted();
        
        TR << "Embed from Topology: Calculation Complete!" << std::endl;
        
		// Done!
		return;
	}
    
	/// @brief      Embedding Search and Score Method
	/// @details    Starting from embedding defined from topology (see method above), score each conformation, perturb,
	///             score, and aver n trials and MC accept lowest. Takes the embed params resource loaded via resource
    ///             loader (RM).
	void
    EmbeddingFactory::embed_from_search_and_score( core::membrane::util::EmbedConfigInfoOP )
	{
		TR << "Embed from Search and Score: Calculating..." << std::endl;
        
        using namespace core::membrane::util;
        
        // DO NOT USE THIS CODE UNTIL THE OTHER SCORING CLASSES ARE TESTED!!!!!
        /**
        // Assign Starting Values from Embedding Configuration
        core::Vector normal = embed_info_->normal;
        core::Vector center = embed_info_->center;
        
        // Call Search normal center method
        MembraneSearch msearch = new MembraneSearch(desc, chain);
        EmbedConfigInfoOP new_config = msearch.search_normal_center( pose_, normal, center );
        
        // Set New Normal and center values
        embedding->normal = new_config->normal;
        embedding->center = new_config->center;
        
        **/
        
        TR << "Embed from Search and Score: Calculation Complete!" << std::endl;
        
        // Done!
        return;
	}
    
    /// @brief   Non Membrane Embedding Method
    /// @details The embedding config info object will contain the data required to initialize this
    ///          information:
    ///             center = center of mass of the non membrane chain
    ///             normal = vector sum of residue-com
    ///             detph = 30A + some user specified factory
    void
    EmbeddingFactory::embed_for_non_membrane( core::membrane::util::EmbedConfigInfoOP )
    {
        
        TR << "Embed for membrane associated chian: Calculating..." << std::endl;
        
        using namespace core::membrane::util;
        using namespace core::membrane::geometry;
        
        /**
        // Calculate center vector (and construct core::Vector)
        numeric::xyzVector< core::Real > com = center_of_mass( *pose_, 1, pose_->total_residue() );
        core::Vector center( com.x(), com.y(), com.z() ); // safe casting
        
        
        // Method of Normal Calculation
        if ( embedding->method )
        // Taking uninitialized data here for now to think about what i need
        utility::vector1< core::Size > residues; // residues to include for com calculation (get from embed config info)
        numeric::xyzVector< core::Real > net_norm = calc_
        core::Vector normal( net_norm.x(), net_norm.y(), net_norm.z() ); // safe casting
        
        // Calculate Depth
        core::Real depth = 30;
        
        // Preset membrane associated factor
        if ( embedding->depth_tag.compare("memrbane_associated") == 0 ) {
            depth = depth;
        }
        
        // Preset tag here
        if ( embedding->depth_tag.compare("tag here") == 0 ) {
            // do something to depth
        }
        
        // Associate rsd with enum type
        
        
        // Calculate Residue Center of mass at CA for the selected CA
        embedding->norma = nnormal;
        embedding->center = center;
        embedding->depth = depth;
    
         **/
        TR << "Embed for membrane associated chain: Complete!" << std::endl;
        
        // Done!
        return;
    }
    
} // core
} // membrane
} // geometry

#endif // INCLUDED_core_membrane_geometry_EmbeddingFactory_cc

