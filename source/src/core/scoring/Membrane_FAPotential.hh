// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		core/scoring/Membrane_FAPotential.cc
///
/// @brief		Membrane FA Potential - Class for Fullatom Membrane Scoring Methods
/// @details	Compute High resolution energy terms and high resolution embedding corrections
///				for penalties. Also contains pass-through methods for accessing and updating
///				mp framework supported data in a membrane conformation.
///				Last Modified: 3/11/14
///
///	@author		Rebecca Faye Alford (rfalford12@gmail.com)
/// @author		Patrick Barth (original)

// Unit headers
#include <core/scoring/Membrane_FAPotential.fwd.hh>
#include <core/scoring/MembranePotential.hh>

// Project headers
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/EnvPairPotential.hh>
#include <core/conformation/Residue.hh>

// Package headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

#include <basic/datacache/CacheableData.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>

namespace core {
	namespace scoring {
		
		////////////////////////////////// Membrane Fullatom Embedding ///////////////////////////////////
		
		/// @brief	 Membrane Fullatom embedding info
		/// @details Cacheable Data - Stores Full atom embedding information including
		///			 projection from z axis, fa depth, center, penalty, membrane thicnkess
		///			 steepness and normal
		class Membrane_FAEmbed : public basic::datacache::CacheableData {
			
		public:
			
			/// @brief Default Constructor
			Membrane_FAEmbed(): calculated_(false) {};
			
			/// @brief Copy Constructor
			Membrane_FAEmbed( Membrane_FAEmbed const & src );
			
			/// @brief Cacheable Data base Mehtod - Clone Object
			basic::datacache::CacheableDataOP clone() const {
				return new Membrane_FAEmbed( *this );
			}
			
			/// @brief Compute FA Proj to Z Axis
			Real & fa_proj( Size const seqpos, Size const atom ) { return fa_proj_[seqpos][atom]; }
			Real fa_proj( Size const seqpos, Size const atom ) const { return fa_proj_[seqpos][atom]; }
			
			/// @brief Compute Depth of Position in Membrane
			Real & fa_depth(Size const seqpos, Size const atom) { return fa_depth_[seqpos][atom]; }
			Real fa_depth(Size const seqpos, Size const atom) const { return fa_depth_[seqpos][atom]; }
			
			/// @brief Compute Derivative of Fa Proj.
			Real & fa_proj_deriv( Size const seqpos, Size const atom ) { return fa_proj_deriv_[seqpos][atom]; }
			Real fa_proj_deriv(Size const seqpos, Size const atom) const { return fa_proj_deriv_[seqpos][atom]; }
			
			/// @brief Get Coordinates (I think it is storing these)
			Vector & fa_proj_coord( Size const seqpos, Size const atom ) { return fa_proj_coord_[seqpos][atom]; }
			Vector fa_proj_coord( Size const seqpos, Size const atom ) const { return fa_proj_coord_[seqpos][atom]; }
			
			/// @brief Compute Fullatom Center
			Real fa_center() const { return fa_center_; }
			Real & fa_center() { return fa_center_; }
			
			Real fa_penalty() const { return fa_penalty_;}
			
			Real & fa_penalty(){ return fa_penalty_;}
			
			Real thickness() const{ return thickness_;}
			
			Real & thickness(){ return thickness_;}
			
			// @brief Probably the same as tm_projection
			Real steepness() const{ return steepness_;}
			
			Real & steepness(){ return steepness_;}
			
			bool calculated() const{ return calculated_;}
			
			bool & calculated(){ return calculated_;}
			
			bool Fa_Membed_update() const{ return Fa_Membed_update_;}
			
			bool & Fa_Membed_update(){ return Fa_Membed_update_;}
			
			void initialize( pose::Pose const & pose );
			
		private: // methods
			
			/// @brief Allocate Memory needed in pose cache?? (rebecca thinks we don't need this)
			void allocate_appropriate_memory( pose::Pose const & pose ) const;
			
		private: // data
			
			// fa projection + depth info
			mutable utility::vector1 < utility::vector1 < Real > > fa_proj_; //pba
			mutable utility::vector1 < utility::vector1 < Real > > fa_depth_; //pba
			mutable utility::vector1 < utility::vector1 < Vector > > fa_proj_coord_; //pba
			mutable utility::vector1 < utility::vector1 < Real > > fa_proj_deriv_; //pba
			
			// Probably for debugging, currently never changed anywhere
			bool calculated_;
			
			// Center, Penalty, Thickness, Steepnessgm
			mutable Real fa_center_;
			mutable Real fa_penalty_;
			Real thickness_; //pba
			Real steepness_; //pba
			
			// Turn on updates
			bool Fa_Membed_update_;
		};
		
		/// @brief		Mmebrane Fullatom Potential - Scoring Class
		/// @details	Helper methods for computing fullatom energy terms
		///				in the membrane scoring function
		class Membrane_FAPotential : public EnvPairPotential {
			
		public:
			
			/// @brief Default Constructor (initialize base class from membrane potential sfxn)
			Membrane_FAPotential():
			membrane_potential_( ScoringManager::get_instance()->get_MembranePotential() ) {};
			
			void compute_fa_projection(pose::Pose & pose) const;
			
			/// @brief Base Class Method - Finalize Scoring Setup
			virtual	void finalize( pose::Pose & pose ) const;
			
		private: // methods
			
			/// @brief Compute Fullatom Projection
			void fa_projection(
							   pose::Pose & pose,
							   Vector const & normal,
							   Vector const & center,
							   Real const & thickness,
							   Real const & steepness,
							   Real const & penalty
							   ) const;
			
		private: // data
			
			bool calculated_; // currently not sure where/why this is called - possibly in loop
			MembranePotential const & membrane_potential_; // const base class
			
		};
		
		/// @brief Grab Const MP Fa Embedding data from the pose cache
		Membrane_FAEmbed const & Membrane_FAEmbed_from_pose( pose::Pose const & );
		
		/// @brief Grab Const MP Fa embedding data from the pose cache
		Membrane_FAEmbed & nonconst_Membrane_FAEmbed_from_pose( pose::Pose & );
		
	} // scoring
} // core

