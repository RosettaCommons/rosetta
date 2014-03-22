// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 	 core/conformation/membrane/MembraneInfo.hh
///
/// @brief 	 Membrane Conformation Info Obejct
/// @details The Membrane Conformation Info Object is responsible for:
///             - maintaining a correct membrane foldtree
///             - maintaining references to the membrane and embedding residues
///             - providing access to membrane related data
///
/// @note    Last Modified 3/12/14
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_conformation_membrane_MembraneInfo_hh
#define INCLUDED_core_conformation_membrane_MembraneInfo_hh

// Unit headers
#include <core/conformation/membrane/MembraneInfo.fwd.hh>
#include <core/conformation/Conformation.hh>

// Project Headers
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/LipidAccInfo.hh>

// Package Headers
#include <core/kinematics/FoldTree.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>

using namespace core::conformation;
using namespace core::conformation::membrane;

/// The goal of the membrane info object is to coordinate information related to
/// the conformaiton of a memrbane protein - as it relates to both kinematics and scoring
/// of a pose. The object lives inside of Conformaiton, but maintains an access pointer
/// to the conformation in order to access its data members and check for updates.
/// This class does not inherit from an observer, but implements the observer pattern

namespace core {
namespace conformation {
namespace membrane {
		
/// @brief A Membrane conformation: Additional data for maintaining a mebrane
/// @details Handles membrane proteins
class MembraneInfo : public utility::pointer::ReferenceCount {
		
	public: // construction, copying, access, invariants
		
		//// Constructors //////////////////////////
		
		/// @brief Standard Constructor
		MembraneInfo(
							 ConformationOP conf_in,
							 utility::vector1< std::pair< int, int > > embres_map,
							 int membrane
							 );
		
		/// @brief Copy Constructor
		MembraneInfo( MembraneInfo const & src );
		
		/// @brief Return Membrane Embedding Map
		utility::vector1< std::pair< int, int > > embres_map() const;
		
		/// @brief Return Membrane root
		int membrane() const;
		
		/// @brief Assert Membrane Conformation Invariants
		bool is_membrane() const;
   
   private: // class helper methods
	
		/// @brief Initialize Membrane Info Class
		void init();
	
		/// @brief Compute WHole Pose Centroid Embedding Parameters
		void init_for_lowres();
	
		/// @brief Compute Whole Pose Fullatom Embedding Parameters
		void init_for_highres();
	
		/// @brief Watch for Membrane Relevant changes in the conformation every time
		///		   membrane info is updated. 
		bool conformation_changed();
	
	public: // membrane embedding data
		
		/// Membrane Postion Info ///////////
		
		/// @brief Get Membrane Center Coords
		core::Vector membrane_center() const;
		
		/// @brief Get Membrane Normal Coords
		core::Vector membrane_normal() const;
		
		/// @brief Get Membrane Thickness Parameter
		core::Real membrane_thickness() const;
	
		/// Chain-Specific Kinematic Embedding Info ///////////
		
		/// @brief Get Chain Embedding Center Coords
		core::Vector embedding_center( core::Size chain ) const;
			
		/// @brief Get Chain Embedding Normal coords
		core::Vector embedding_normal( core::Size chain ) const;
	
		/// @brief Get Chain Embedding Depth Parameter
		core::Real embedding_depth( core::Size chain ) const;
	
		/// Whole-Chain Non-Kinematic Embedding Info ////////////////
	
		/// @brief WHole Pose Center
		core::Vector pose_embedding_center();
		
		/// @brief WHole Pose Normal
		core::Vector pose_embedding_normal();
	
		/// @brief Depth
		core::Real residue_depth( core::Size seqpos );
		
		/// Fullatom Whole-Chain Non-Kinematic Membrane Info ///////////////
	
		/// @brief WHole Pose Steepness
		core::Real pose_embedding_steepness();
	
		/// @brief Get Fullatom TM Projection at Seqpos and Atom Pos
		core::Real fa_proj( core::Size seqpos, core::Size atom );
	
		/// @brief Get Depth of FA Coordinate
		core::Real fa_depth( core::Size seqpos, core::Size atom );
		
		/// @brief Get Fullatom Projection Derivative
		core::Real fa_proj_deriv( core::Size seqpos, core::Size atom );
	
		/// @brief Return Coordinate of the FA Proj
		core::Vector fa_proj_coord( core::Size seqpos, core::Size atom );
			
		/// Update Fullatom or Centroid Info Set ///////////////
	
		/// @brief Update Whole Pose for High Resolution Typesets
		void update_fullatom_info();
	
		/// @brief Update HWole Pose Embedding Info for Low Resolution Typesets
		void update_lowres_info();
				
	public: // membrane conformation data
	
		/// @brief Return the total number of polymer reisdues in the pose (no mp residues)
		core::Size total_polymer_residue() const;
	
		/// @brief Return the total number of polymer chains in the pose (no mp chain)
		core::Size num_polymer_chains() const;
		
		/// @brief Add Spanning Topology
		void add_topology_by_chain( SpanningTopology sp, core::Size chain );
		utility::vector1< SpanningTopology > spanning_topology() const;
		
		/// @brief Add Lipid Accessibility Info
		void add_lips_by_chain( LipidAccInfo const & sp, core::Size chain );
		utility::vector1< LipidAccInfo > lipid_acc_data() const;
		
	private: // methods
	
		/// @brief Default Constructor
		MembraneInfo();
		
		//// FoldTree ////////////////////////////
		
		/// @brief Build Membrane FoldTree
	void build_membrane_foldtree();
		
	private: // data
	
		// Store an Access Pointer to Conformation
		ConformationOP conf_;
	
		/// Non-Kinematic Embedding Parameters /////////////
		
		// Normal/Center
		core::Vector center_;
		core::Vector normal_;
	
		// Steepness
		core::Real steepness_;
	
		// Lowres Depth Info
		utility::vector1< core::Real > depth_;
	
		// Fullatom TM Projection Info
		utility::vector1 < utility::vector1 < core::Real > > fa_proj_;
		utility::vector1 < utility::vector1 < core::Real > > fa_depth_;
		utility::vector1 < utility::vector1 < core::Vector > > fa_proj_coord_;
		utility::vector1 < utility::vector1 < core::Real > > fa_proj_deriv_;
		
		/// Data for Object Syncrhonization ////////////////
	
		// Determine if we are working with a fullatom pose
		bool fullatom_;
		core::Size nres_;
	
		/// Kinematic Embedding Info //////////////////
	
		// Store Embedding Residue Map
		utility::vector1< std::pair< int, int > > embres_map_;
		
		// Store index of the membrane origin
		int membrane_;
	
		/// Membrane Conformation Info ///////////////
		
		// Lipid accessibility data
		utility::vector1< LipidAccInfo > lipid_acc_data_;
		utility::vector1< SpanningTopology > spanning_topology_;
		
	}; // MembraneInfo
		
} // membrane
} // conformation
} // core

#endif // INCLUDED_core_membrane_MembraneInfo_hh

