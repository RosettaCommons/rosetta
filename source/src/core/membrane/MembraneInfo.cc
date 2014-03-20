// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 	 core/membrane/MembraneInfo.cc
///
/// @brief 	 Membrane Conformation Info Object
/// @details The Membrane Conformation Info Object is responsible for:
///             - maintaining a correct membrane foldtree
///             - maintaining references to the membrane and embedding residues
///             - providing access to membrane related data
///
/// @note    Last Modified 3/12/14
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_MembraneInfo_cc
#define INCLUDED_core_membrane_MembraneInfo_cc

// Unit Headers
#include <core/membrane/MembraneInfo.hh>

// Project Headers
#include <core/membrane/properties/SpanningTopology.hh>
#include <core/membrane/properties/LipidAccInfo.hh>

// Package Headers
#include <core/conformation/Conformation.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/membrane/util/Exceptions.hh>

#include <core/conformation/util.hh>
#include <core/types.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>

using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.membrane.MembraneInfo");

using namespace core::conformation;

namespace core {
namespace membrane {

	//// Constructors //////////////////////////
	
	/// @brief Default Constructor
	MembraneInfo::MembraneInfo() :
		utility::pointer::ReferenceCount(),
		conf_(NULL),
		center_(0, 0, 0),
		normal_(0, 0, 0),
		steepness_(0),
		membrane_(1)
	{
		init();
	}

	/// @brief Standard Constructor
	MembraneInfo::MembraneInfo(
		ConformationOP conf_in,
		utility::vector1< std::pair< int, int > > embres_map,
		int membrane
		) :
		utility::pointer::ReferenceCount(),
		conf_( conf_in ),
		center_(0, 0, 0),
		normal_(0, 0, 0),
		steepness_(0),
		embres_map_( embres_map ),
		membrane_( membrane )
	{
		init();
	}

	/// @brief Copy Constructor
	MembraneInfo::MembraneInfo( MembraneInfo const & src ) :
		utility::pointer::ReferenceCount(),
		conf_( src.conf_ ),
		center_( src.center_ ),
		normal_( src.normal_ ),
		steepness_( src.steepness_ ),
		embres_map_( src.embres_map_ ),
		membrane_( src.membrane_ )
	{
		init();
	}


	/// @brief Return Membrane Embedding Map
	utility::vector1< std::pair< int, int > >
	MembraneInfo::embres_map() const { return embres_map_; }

	/// @brief Return Membrane root
	int
	MembraneInfo::membrane() const { return membrane_; }

	/// @brief Assert Membrane Conformation Invariants
	bool MembraneInfo::is_membrane() const {
		
		// Check Class Invariants
		
		// Check membrane root is within the bounds of the pose
		if ( membrane_ <= 0 || membrane_ > (int) conf_->size() ) {
			TR << "Membrane root specified" << membrane_ << "must be between 1 and " << conf_->size() << std::endl;
			return false;
		}
		
		// Check that there is an embedding residue anchored by jump for each protein chain
		if ( embres_map_.size() != conf_->num_chains()-1 ) {
			TR << "There must exist an embedding residue for every protein chain in the pose" << std::endl;
			TR << "Specified " << embres_map_.size() << " embedding residues but needs " << conf_->num_chains()-1 << std::endl;
			return false;
		}
		
		// Check that all of the jump anchors are within the bounds of the pose
		for ( core::Size i = 1; i <= conf_->num_chains()-1; ++i ) {
			if ( embres_map_[ i ].first <= 0 || embres_map_[ i ].first > (int) conf_->size() ) {
				TR << "Embedding residue specified must be within the pose (between 1 and " << conf_->size() << ")" << std::endl;
				TR << "Specified " << embres_map_[ i ].first << std::endl;
				return false;
			}
		}
		
		// Check Fold tree
		if ( !conf_->fold_tree().check_fold_tree() ) {
			TR << "Membrane foldtree is an invalid foldtree!" << std::endl;
			return false;
		}
		
		// Check the membrane embedding residue is the root of the pose
		if ( !conf_->fold_tree().is_root( membrane_ ) ) {
			TR << "Membrane root " << membrane_ << " must be the root of the foldtree!" << std::endl;
			return false;
		}
		
		// Else return true
		return true;
	}
	
	/////////////////////////////// Membrane Info Class Helper Methods ///////////////////////////
	
	/// @brief Initialize Membrane Info Class
	void MembraneInfo::init() {
		
		// Initialize Base info from Conformation
		fullatom_ = conf_->is_fullatom();
		nres_ = conf_->size();
		
		// Based on Residue Typeset, update whole pose embedding info
		if ( fullatom_ ) {
			
			// Initialize Base Members
			fa_proj_.resize(nres_);
			fa_depth_.resize(nres_);
			fa_proj_coord_.resize(nres_);
			fa_proj_deriv_.resize(nres_);
			
			// not yet implemented
			// should initialize: center, normal, steepness, all fa params
			
		} else {
			
			depth_.resize(nres_);
			
			// not yet implemented
			// should initialize: center, normal, depth array
			
		}
		
		// Build fold tree, check is a membrane conf
		build_membrane_foldtree();
		assert( is_membrane() );
		assert( !conformation_changed() );
	}
	
	/// @brief Compute WHole Pose Centroid Embedding Parameters
	void MembraneInfo::init_for_lowres() {
		// not yet implemented
	}
	
	/// @brief Compute Whole Pose Fullatom Embedding Parameters
	void MembraneInfo::init_for_highres() {
		// not yet implemented
	}
	
	/// @brief Watch for Membrane Relevant changes in the conformation every time
	///		   membrane info is updated.
	bool MembraneInfo::conformation_changed() {
		
		using namespace core::membrane::util;
		
		// Keep track of changes
		bool changed = false;
		
		// Grab new conformation info
		core::Size new_nres_ = conf_->size();
		bool new_fa_ = conf_->is_fullatom();
		
		// Check for differences and respond appropriately
		if ( !new_fa_ && fullatom_ ) {
			init_for_lowres();
			changed = true;
		} else if ( new_fa_ && !fullatom_ ) {
			init_for_highres();
			changed = true;
		}
		
		// Check for pose length changed
		if ( nres_ != new_nres_ ) {
			throw new EXCN_NonMembrane("Cannot change length of a membrane pose!");
		}
		
		return changed;
	}

	
	////////////////////////////////// Membrane Embedding Information ///////////////////////////////
	
	/// Membrane Postion Info ///////////
	
	/// @brief Get Membrane Center Coords
	core::Vector
	MembraneInfo::membrane_center() const {
		return conf_->residue( membrane_ ).atom( 2 ).xyz();
	}
	
	/// @brief Get Membrane Normal Coords
	core::Vector
	MembraneInfo::membrane_normal() const {
		return conf_->residue( membrane_ ).atom( 1 ).xyz();
	}
	
	/// @brief Get Membrane Thickness Parameter
	core::Real
	MembraneInfo::membrane_thickness() const {
		return conf_->residue( membrane_ ).atom( 3 ).xyz().y();
	}
	
	//// Chain Specific Embedding Residues
	
	/// @brief Get Chain Embedding Center Coords
	core::Vector
	MembraneInfo::embedding_center( core::Size chain ) const {
		core::Size resnum = embres_map_[ chain ].second;
		return conf_->residue( resnum ).atom( 2 ).xyz();
	}
	
	/// @brief Get Chain Embedding Normal coords
	core::Vector
	MembraneInfo::embedding_normal( core::Size chain ) const {
		core::Size resnum = embres_map_[ chain ].second;
		return conf_->residue( resnum ).atom( 1 ).xyz();
	}
	
	/// @brief Get Chain Embedding Depth Parameter
	core::Real
	MembraneInfo::embedding_depth( core::Size chain ) const {
		core::Size resnum = embres_map_[ chain ].second;
		return conf_->residue( resnum ).atom( 3 ).xyz().y();
	}
	
	/// Whole-Chain Non-Kinematic Embedding Info ////////////////
	
	/// @brief WHole Pose Center
	core::Vector
	MembraneInfo::pose_embedding_center() { return center_; }
	
	/// @brief WHole Pose Normal
	core::Vector
	MembraneInfo::pose_embedding_normal() { return normal_;}
	
	/// @brief Depth
	core::Real 
	MembraneInfo::residue_depth( core::Size seqpos ) { return depth_[seqpos]; }
	
	/// Fullatom Whole-Chain Non-Kinematic Membrane Info ///////////////
	
	/// @brief WHole Pose Steepness
	core::Real
	MembraneInfo::pose_embedding_steepness() {
		if ( !fullatom_ ) {
			TR << "Warning: Fullatom Membrane Conformation Data Uninitialized or Updated for a Centroid Pose!" << std::endl;
		}
		return steepness_;
	}
	
	/// @brief Get Fullatom TM Projection at Seqpos and Atom Pos
	core::Real
	MembraneInfo::fa_proj( core::Size seqpos, core::Size atom ) {
		if ( !fullatom_ ) {
			TR << "Warning: Fullatom Membrane Conformation Data Uninitialized or Updated for a Centroid Pose!" << std::endl;
		}
		return fa_proj_[seqpos][atom];
	}
	
	/// @brief Get Depth of FA Coordinate
	core::Real MembraneInfo::fa_depth( core::Size seqpos, core::Size atom ) {
		if ( !fullatom_ ) {
			TR << "Warning: Fullatom Membrane Conformation Data Uninitialized or Updated for a Centroid Pose!" << std::endl;
		}
		return fa_depth_[seqpos][atom];
	}
	
	/// @brief Get Fullatom Projection Derivative
	core::Real
	MembraneInfo::fa_proj_deriv( core::Size seqpos, core::Size atom ) {
		if ( !fullatom_ ) {
			TR << "Warning: Fullatom Membrane Conformation Data Uninitialized or Updated for a Centroid Pose!" << std::endl;
		}
		return fa_proj_deriv_[seqpos][atom];
	}
	
	/// @brief Return Coordinate of the FA Proj
	core::Vector
	MembraneInfo::fa_proj_coord( core::Size seqpos, core::Size atom ) {
		if ( !fullatom_ ) {
			TR << "Warning: Fullatom Membrane Conformation Data Uninitialized or Updated for a Centroid Pose!" << std::endl;
		}
		return fa_proj_coord_[seqpos][atom];
	}
	
	/// Update Fullatom or Centroid Info Set ///////////////
	
	/// @brief Update Whole Pose for High Resolution Typesets
	void MembraneInfo::update_fullatom_info() {
		// not yet implemented
	}
	
	/// @brief Update HWole Pose Embedding Info for Low Resolution Typesets
	void MembraneInfo::update_lowres_info() {
		// not yet implemented
	}
	
	/////////////////////////// Access Membrane Conformation Info ///////////////////////////
	
	/// @brief Return the total number of polymer residues in the pose (excluding mp residues)
	core::Size
	MembraneInfo::total_polymer_residue() const {
		core::Size nres = conf_->num_chains();
		return conf_->size()-nres;
	}
	
	/// @brief Return the total number of polymer chains in the pose (excluding mp residues)
	core::Size
	MembraneInfo::num_polymer_chains() const {
		return conf_->num_chains()-1;
	}
	
	/// @brief Add Spanning Topology
	void
	MembraneInfo::add_topology_by_chain( core::membrane::properties::SpanningTopology sp, core::Size ) {
		spanning_topology_.push_back( sp );
	}
	
	utility::vector1< properties::SpanningTopology >
	MembraneInfo::spanning_topology() const {
		return spanning_topology_;
	}
	
	/// @brief Add Lipid Accessibility Info
	void
	MembraneInfo::add_lips_by_chain( properties::LipidAccInfo const & sp, core::Size ) {
		lipid_acc_data_.push_back( sp );
	}
	
	utility::vector1< properties::LipidAccInfo >
	MembraneInfo::lipid_acc_data() const {
		return lipid_acc_data_;
	}
	
	//// FoldTree ////////////////////////////
	
	/// @brief Build Membrane FoldTree
	void
	MembraneInfo::build_membrane_foldtree() {
		
		using namespace core::kinematics;
		
		// Construct a new foldtree
		FoldTreeOP ft = new FoldTree();
		
		// Add Initial membrane edge and maintain a jump counter
		int jump_counter = 1;
		ft->add_edge( 1, membrane_, jump_counter );
		jump_counter++;
		
		// Add peptide edges
		for ( core::Size i = 1; i <= conf_->num_chains()-1; ++i ) {
			
			// Add edge for each chain
			ft->add_edge( conf_->chain_begin(i), conf_->chain_end(i), -1 );
			
			// Add a chain-connecting edge if between chains
			if ( i != conf_->num_chains()-1 ) {
				ft->add_edge( conf_->chain_begin(i), conf_->chain_begin(i+1), jump_counter );
				jump_counter++;
			}
		}
		
		// Set the membrane to the root of the foldtree. This must occur before adding embedding edges
		// because the embedding edges are added to a non-vertex posiiton and will be lost in the reorder arbitrarily
		ft->reorder(membrane_);
		
		// Add membrane edges
		for ( core::Size i = 1; i <= embres_map_.size(); ++i ) {
			
			// Add edge between each chain beginning and its corresponding embedding
			ft->add_edge( embres_map_[ i ].first, embres_map_[ i ].second, jump_counter );
			jump_counter++;
		}
		
		// Setup the new foldtree
		conf_->fold_tree( *ft );
	}
	
} // membrane
} // core

#endif // INCLUDED_core_membrane_MembraneInfo_cc
