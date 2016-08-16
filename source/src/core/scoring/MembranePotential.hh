// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/scoring/MembranePotential.hh
///
/// @brief  Membrane Potential - Base Scoring Methods for LowRes Energy Function
/// @details Compute Low Res membrane energy terms: Menv, MPair, MCBeta and Membrane
///    penalties. Also contains pass-through methods for accessing and updating
///    mp framework supported data in a membrane conformation.
///    Last Modified: 3/11/14
///
/// @author  Rebecca Faye Alford (rfalford12@gmail.com)
/// @author  Julia Koehler Leman (julia.koehler1982@gmail.com)
/// @author  Bjorn Wallner (Original)

#ifndef INCLUDED_core_scoring_MembranePotential_hh
#define INCLUDED_core_scoring_MembranePotential_hh

// Unit Headers
#include <core/scoring/MembranePotential.fwd.hh>
#include <core/scoring/EnvPairPotential.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <core/conformation/Residue.hh>
#include <basic/datacache/CacheableData.hh>

// Utility headers
#include <numeric/xyzVector.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray4D.hh>

#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {

////////////////////////////////// Membrane Embedding ///////////////////////////////////

/// @brief   Whole Pose Membrane Embedding
/// @details Define the embedding of the membrane pose based on computed
///    normal and center parameters. These are initialzed in the membrane protein
///    framework and then recomputed based upon the structured and stored in MP residues
///    (see MP Framework code)
class MembraneEmbed : public basic::datacache::CacheableData {

public:

	/// @brief Default Constructor
	MembraneEmbed(): calculated_(false), spanning_(false) {};

	/// @brief Copy Constructor
	MembraneEmbed( MembraneEmbed const & src );

	/// @brief Clone Cacheable Data
	basic::datacache::CacheableDataOP clone() const {
		return basic::datacache::CacheableDataOP( new MembraneEmbed( *this ) );
	}

	/// @brief Compute Size of MP (??)
	Size size() const{ return depth_.size(); }

	/// @brief Compute depth of residue in the membrane
	Real & depth( Size const seqpos ) { return depth_[seqpos]; }

	/// @brief Set Pose Embedding Normal (should not use this method, deprecated 3/11/14)
	void set_normal( Vector const & v ) { normal_ = v; }

	/// @brief Set Pose Embedding Center (should not use this method, deprecated 3/11/14)
	void set_center( Vector const & v ) { center_ = v; }

	/// @brief Set Penalty
	void set_penalty( Real const & p ) { penalty_ = p; }

	/// @brief Get NonConst Depth
	Real depth( Size const seqpos ) const {
		return depth_[seqpos];
	}

	/// @brief Return Ref Spanning Parameter
	bool & spanning() { return spanning_; }

	/// @brief Return Non_Ref Spanning Parameter
	bool spanning() const { return spanning_; }

	/// @brief Return Calculated (no idea what this does)
	bool calculated() const { return calculated_; }

	/// @brief Return Non_Ref Calculated (no idea what this does - maybe observer)
	bool & calculated() { return calculated_; }

	/// @brief Get Normal Parameter
	Vector const & normal() const { return normal_; }

	/// @brief Get Center Parameter
	Vector const & center() const { return center_; }

	/// @brief Get MP Penalty
	Real const & penalty() const { return penalty_; }

	/// @brief Initialize Membrane Embedding From Pose??
	void initialize( pose::Pose const & pose );

private: // data

	utility::vector1 < Real > depth_;
	Vector normal_;
	Vector center_;

	bool calculated_;
	bool spanning_;

	// unused Size tm_projection_;
	Real penalty_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // MembraneEmbed

////////////////////////////////// Membrane Potential ////////////////////////////////////

/// @brief   Rosetta Membrane Low Resolution Scoring Methods
/// @details Compute scoring terms part of the Membrane Low resolution energy function. Developed
///    by Vladmir Yarov-Yaravoy et al. 2006. Includes Menv, MPair, MCBeta, and membrane
///    alpha helical specific penalties. Framework tied.
class MembranePotential : public EnvPairPotential {

public:

	/// @brief Default Constructor
	MembranePotential();

	/// @brief Evalaute Membrane Environment
	void evaluate_env(
		pose::Pose const & pose,
		conformation::Residue const & rsd,
		Real const MembraneDepth,
		Real & membrane_env_score
	) const;

	void
	evaluate_env(
		pose::Pose const & pose,
		conformation::Residue const & rsd,
		Real & membrane_env_score
	) const;


	/// @brief Evaluate CBeta Score (no idea...?)
	void evaluate_cbeta(
		pose::Pose const & pose,
		conformation::Residue const & rsd,
		Real & membrane_cb_score
	) const;

	/// @brief Evaluate Energy For Two Residues
	void evaluate_pair(
		pose::Pose const & pose,
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		Real const cendist,
		Real & membrane_pair_score
	) const;

	/// I think these methods should be private ///

	/// @brief Finalize Setup of MP Potential Class
	virtual void finalize( pose::Pose & pose ) const;

	/// @brief Compute Membrane Embedding from pose and add Membed to Pose Cache
	void compute_membrane_embedding( pose::Pose & pose ) const;

	/// @brief Initialize Membrane Center/Normal
	void init_membrane_center_normal(
		pose::Pose const & pose,
		Vector & normal,
		Vector & center
	) const;


	/// @brief Compute Transmembrane Spanning Projection Penalty (documentation for what this is - se MP potential from refactor)
	void tm_projection_penalty(
		pose::Pose const & pose,
		Real & tm_proj
	) const;

	void tm_projection_penalty(
		pose::Pose const & pose,
		Vector const & normal,
		Vector const & center,
		Real & tm_proj
	) const;

	/// @brief Compute penaly for alpha helices that are not tm spanning
	void non_helix_in_membrane_penalty(
		pose::Pose const & pose,
		Real & non_helix_in_membrane_penalty
	) const;

	void non_helix_in_membrane_penalty(
		pose::Pose const & pose,
		Vector const & normal,
		Vector const & center,
		Real & non_helix_in_membrane_penalty
	) const;

	/// @brief Compute penalty for ???
	void termini_penalty(
		pose::Pose const & pose,
		Real & termini_penalty
	) const;

	void termini_penalty(
		pose::Pose const & pose,
		Vector const & normal,
		Vector const & center,
		Real & termin_penalty
	) const;

	/// @brief User Specified use penalties
	bool Menv_penalties() const
	{
		return Menv_penalties_;
	}

	/// @brief Initialize Membrane Embedding
	bool Membed_init() const { return Membed_init_; }


private:

	/// @brief Score normal and center with environment score
	void score_normal_center(
		pose::Pose const & pose,
		Vector const & normal,
		Vector const & center,
		Real & score
	) const;

	/// @brief Helper function to determine normal vector
	void search_memb_normal(
		Vector & n,
		Real const & alpha,
		Real const & theta
	) const;

	/// @brief Helper function to determine membrane center
	void search_memb_center(
		Vector & c,
		Vector & n,
		Real const & delta) const; //vmyy

	/// @brief Randomnly Rotate Vector
	void rot_perturb_vector(
		Vector & v,
		Real const & std_dev
	) const;

	/// @brief Randomnly Translate Vector
	void rigid_perturb_vector(
		Vector & v,
		Real const & std_dev
	) const;

	/// @brief Check Spanning ( should use new adapted check spanning in metrics or geom)
	bool check_spanning(
		pose::Pose const & pose,
		Vector const & normal,
		Vector const & center
	) const;

private: // data

	/// Membrane Etable Data
	ObjexxFCL::FArray3D< Real > mem_env_log6_;
	ObjexxFCL::FArray3D< Real > mem_env_log10_;
	ObjexxFCL::FArray1D< Real > mem_cbeta_den6_;
	ObjexxFCL::FArray1D< Real > mem_cbeta_den12_;
	ObjexxFCL::FArray1D< Real > mem_cbeta_2TM_den6_;
	ObjexxFCL::FArray1D< Real > mem_cbeta_2TM_den12_;
	ObjexxFCL::FArray1D< Real > mem_cbeta_4TM_den6_;
	ObjexxFCL::FArray1D< Real > mem_cbeta_4TM_den12_;
	ObjexxFCL::FArray4D< Real > mem_pair_log_;

	/// Centroid Database
	Real const cen_dist5_pad;
	//unused Real const cen_dist6_pad;
	Real const cen_dist7_pad;
	Real const cen_dist10_pad;
	Real const cen_dist12_pad;

	Real const cen_dist5_pad_plus ;
	// unused Real const cen_dist6_pad_plus ;
	Real const cen_dist7_pad_plus ;
	Real const cen_dist10_pad_plus;
	// unused Real const cen_dist12_pad_plus;

	Real const cen_dist5_pad_minus ;
	Real const cen_dist7_pad_minus ;
	Real const cen_dist10_pad_minus;
	Real const cen_dist12_pad_minus;

	Real const cen_dist5_pad_hinv ;
	// unused Real const cen_dist6_pad_hinv ;
	Real const cen_dist7_pad_hinv ;
	Real const cen_dist10_pad_hinv;
	Real const cen_dist12_pad_hinv;

	/// User Specified options for Scoring
	// unused bool calculated_;
	bool no_interpolate_Mpair_;
	bool Menv_penalties_;
	bool Membed_init_;
	bool memb_center_search_;
	bool memb_normal_search_;

	/// MCM Search Parameters for normal/embed caluclations
	Size membrane_center_max_delta_;
	Size membrane_normal_start_angle_;
	Size membrane_normal_delta_angle_;
	Size membrane_normal_max_angle_;
	Size membrane_normal_cycles_;
	Real membrane_normal_magnitude_;
	Real membrane_center_magnitude_;
	Real smooth_move_frac_;
};

/// @brief Add Const Membrane Embedding to the pose cache
MembraneEmbed const & MembraneEmbed_from_pose( pose::Pose const & pose );

/// @brief Add Non Const Membrane Embedding to the pose cache
MembraneEmbed & nonconst_MembraneEmbed_from_pose( pose::Pose & pose );

} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_MembranePotential )
#endif // SERIALIZATION


#endif // INCLUDED_core_scoring_MembranePotential_hh
