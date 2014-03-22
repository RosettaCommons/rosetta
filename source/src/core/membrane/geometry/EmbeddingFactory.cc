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
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/membrane/geometry/util.hh>
#include <core/conformation/membrane/definitions.hh>
#include <core/conformation/membrane/Exceptions.hh>

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
using namespace core::conformation::membrane;

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
                                       SpanningTopologyOP topology
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
        
        // throw exception for unreasonable tag???
        
        // Create an instance of the Membrane Residue Factory
        MembraneResidueFactory mrf;
        
        // Declare accept parameters
        EmbedConfigInfoOP accept = new EmbedConfigInfo;
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
    EmbeddingFactory::embed_from_pdb( EmbedConfigInfoOP embedding ) {
        
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
    EmbeddingFactory::embed_from_topology( EmbedConfigInfoOP embedding  ) {
        
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
    EmbeddingFactory::embed_from_search_and_score( EmbedConfigInfoOP )
	{
		TR << "Embed from Search and Score: Computing..." << std::endl;
        
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
    
	/// @brief		Compute Fullatom Projections for Fa Embedding
	///	@details	Computes Projections, depth, and stores coordinates. Computes fullatom
	///				center and corresponding normal from these computations
	void EmbeddingFactory::embed_for_FA( EmbedConfigInfoOP ) {
		
	}
	
	///////// Helper Methods for Computing Embeddings /////////////////////////////
	/**
	/// @brief Score Pose from Normal and Center
	void score_normal_center( core::Vector & normal, core::Vector & center, core::Real & score ) {

		// Compute Pose Total residue and get topology from pose

		
		
	}
	

		Size const nres=pose.total_residue();
		MembraneTopology const & topology( MembraneTopology_from_pose(pose) );
		score=0;
		Real residue_score(0);
		Real tm_projection(0);
		Real non_helix_pen(0);
		Real termini_pen(0);
		for ( Size i = 1; i <= nres; ++i ) {
			Size rsdSeq(i);
			if ( core::pose::symmetry::is_symmetric( pose ) ) {
				using namespace core::conformation::symmetry;
				SymmetricConformation const & symm_conf (
														 dynamic_cast< SymmetricConformation const & > ( pose.conformation() ) );
				SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );
				if (!symm_info->bb_is_independent(pose.residue(i).seqpos())) {
					rsdSeq = symm_info->bb_follows(pose.residue(i).seqpos());
				}
				if (symm_info->is_virtual(i)) {
					rsdSeq = 0;
				}
			}
			if (rsdSeq ==0 ) continue;
			
			//CA coords
			if ( pose.residue(rsdSeq).aa() == core::chemical::aa_vrt ) continue;
			if(!topology.allow_scoring(rsdSeq)) continue;
			
			Vector const & xyz( pose.residue( i ).atom( 2 ).xyz());
			core::Real depth=dot(xyz-center,normal)+30;
			evaluate_env(pose,pose.residue(i),depth,residue_score);
			score+=residue_score;
		}
		if(Menv_penalties_) {
			tm_projection_penalty(pose,normal,center,tm_projection);
			non_helix_in_membrane_penalty(pose,normal,center,non_helix_pen);
			termini_penalty(pose,normal,center,termini_pen);
			score+=tm_projection+non_helix_pen+termini_pen; // bw skipping term_penalty+50.0*term_penalty; //bw 0.5*c++ version.
		}
	}
	
	void
	MembranePotential::search_memb_normal(Vector & n,
										  Real const & alpha,
										  Real const & theta) const
	{
		Real r_alpha = numeric::conversions::radians(alpha);
		Real r_theta = numeric::conversions::radians(theta);
		Vector u(std::sin(r_alpha) * std::cos(r_theta), std::sin(r_alpha) * std::sin(r_theta), std::cos(r_alpha));
		n=rotation_matrix_degrees(u,alpha)*n;
	}
	
	void
	MembranePotential::search_memb_center(Vector & c,
										  Vector & n,
										  Real const & delta) const
	{
		c=c+delta*n;
	}
	
	void
	MembranePotential::rot_perturb_vector(Vector & v,
										  Real const & std_dev) const
	{
		Vector u(numeric::random::gaussian(),numeric::random::gaussian(),numeric::random::gaussian()); //bw rotation_matrix will normalize.
		Real alpha(numeric::random::gaussian()*std_dev);
		v=rotation_matrix(u,alpha)*v;
	}
	
	void
	MembranePotential::rigid_perturb_vector(Vector & v,
											Real const & std_dev) const
	{
		Vector u(numeric::random::gaussian(),numeric::random::gaussian(),numeric::random::gaussian());
		u.normalize();
	**/
    
} // core
} // membrane
} // geometry

#endif // INCLUDED_core_membrane_geometry_EmbeddingFactory_cc

