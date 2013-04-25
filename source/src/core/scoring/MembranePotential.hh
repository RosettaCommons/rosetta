// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/MembranePotential.hh
/// @brief  Membrane Potential
/// @author Bjorn Wallner


#ifndef INCLUDED_core_scoring_MembranePotential_hh
#define INCLUDED_core_scoring_MembranePotential_hh

#include <core/types.hh>

// Unit headers
#include <core/scoring/MembranePotential.fwd.hh>
// AUTO-REMOVED #include <core/scoring/MembraneTopology.fwd.hh>
// AUTO-REMOVED #include <core/scoring/MembraneTopology.hh>
#include <core/scoring/EnvPairPotential.hh>

// Package headers
// AUTO-REMOVED #include <core/scoring/EnergyGraph.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <basic/datacache/CacheableData.hh>

// Utility headers
// AUTO-REMOVED #include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray4D.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {

class MembraneEmbed : public basic::datacache::CacheableData {

public:
	MembraneEmbed(): calculated_(false), spanning_(false) {};
	MembraneEmbed( MembraneEmbed const & src );

	basic::datacache::CacheableDataOP
	clone() const
	{
		return new MembraneEmbed( *this );
	}

	Size
	size() const{
		return depth_.size();
	}

	Real &
	depth(Size const seqpos)
	{
		return depth_[seqpos];
	}

	void
	set_normal(Vector const & v)
	{
		normal_=v;
	}
	void
	set_center(Vector const & v)
	{
		center_=v;
	}
  void
  set_penalty(Real const & p)
  {
    penalty_=p;
  }

	Real
	depth(Size const seqpos) const
	{
		return depth_[seqpos];
	}

	bool &
	spanning()
	{
		return spanning_;
	}

	bool
	spanning() const
	{
		return spanning_;
	}

	bool
	calculated() const
	{
		return calculated_;
	}

	bool &
	calculated()
	{
		return calculated_;
	}

	Vector const &
	normal() const
	{
		return normal_;
	}

	Vector const &
	center() const
	{
		return center_;
	}

  Real const &
  penalty() const
  {
    return penalty_;
  }

	void
	initialize( pose::Pose const & pose );

private:
	utility::vector1 < Real > depth_;
	Vector normal_;
	Vector center_;
	bool calculated_;
	bool spanning_;
	Size tm_projection_;
  Real penalty_;
};

class MembranePotential : public EnvPairPotential {

public:
	MembranePotential();

	///
	void
	evaluate_env(
		pose::Pose const & pose,
		conformation::Residue const & rsd,
		Real const MembraneDepth,
		Real & membrane_env_score
	) const;
	//

	void
	evaluate_env(
		pose::Pose const & pose,
		conformation::Residue const & rsd,
		Real & membrane_env_score
	) const;

		///
	void
	evaluate_cbeta(
		pose::Pose const & pose,
		conformation::Residue const & rsd,
		Real & membrane_cb_score
	) const;

	///Where the action happens. Calcluate the energy for the two residues
	void
	evaluate_pair(
		pose::Pose const & pose,
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		Real const cendist,
		Real & membrane_pair_score
	) const;

	virtual void
	finalize( pose::Pose & pose ) const;

	void
	compute_membrane_embedding( pose::Pose & pose) const;

	void
	init_membrane_center_normal(pose::Pose const & pose,
															Vector & normal,
															Vector & center) const;


	void tm_projection_penalty(pose::Pose const & pose,Real & tm_proj) const;
	void tm_projection_penalty(pose::Pose const & pose, Vector const & normal,Vector const & center,Real & tm_proj) const;

	void non_helix_in_membrane_penalty(pose::Pose const & pose, Real & non_helix_in_membrane_penalty) const;
	void non_helix_in_membrane_penalty(pose::Pose const & pose, Vector const & normal,Vector const & center,Real & non_helix_in_membrane_penalty) const;

	void termini_penalty(pose::Pose const & pose, Real & termini_penalty) const;
	void termini_penalty(pose::Pose const & pose, Vector const & normal,Vector const & center,Real & termin_penalty) const;

	bool
	Menv_penalties() const
	{
		return Menv_penalties_;
	}

  bool
  Membed_init() const
  {
    return Membed_init_;
  }

protected:

private:

	void score_normal_center(pose::Pose const & pose,
													 Vector const & normal,
													 Vector const & center,
													 Real & score) const;
	void search_memb_normal(Vector & n,
													Real const & alpha,
													Real const & theta) const; //vmyy
	void search_memb_center(Vector & c,
													Vector & n,
													Real const & delta) const; //vmyy
	void rot_perturb_vector(Vector & v, Real const & std_dev) const;  //bw does not belong here
	void rigid_perturb_vector(Vector & v, Real const & std_dev) const;  //bw does not belong here
	bool check_spanning(pose::Pose const & pose, Vector const & normal,Vector const & center) const;

private: // data

	ObjexxFCL::FArray3D< Real > mem_env_log6_;
	ObjexxFCL::FArray3D< Real > mem_env_log10_;
	ObjexxFCL::FArray1D< Real > mem_cbeta_den6_;
	ObjexxFCL::FArray1D< Real > mem_cbeta_den12_;
	ObjexxFCL::FArray1D< Real > mem_cbeta_2TM_den6_;
	ObjexxFCL::FArray1D< Real > mem_cbeta_2TM_den12_;
	ObjexxFCL::FArray1D< Real > mem_cbeta_4TM_den6_;
	ObjexxFCL::FArray1D< Real > mem_cbeta_4TM_den12_;
	ObjexxFCL::FArray4D< Real > mem_pair_log_;

	Real const cen_dist5_pad;
	Real const cen_dist6_pad;
	Real const cen_dist7_pad;
	Real const cen_dist10_pad;
	Real const cen_dist12_pad;

	Real const cen_dist5_pad_plus ;
	Real const cen_dist6_pad_plus ;
	Real const cen_dist7_pad_plus ;
	Real const cen_dist10_pad_plus;
	Real const cen_dist12_pad_plus;

	Real const cen_dist5_pad_minus ;
	Real const cen_dist7_pad_minus ;
	Real const cen_dist10_pad_minus;
	Real const cen_dist12_pad_minus;

	Real const cen_dist5_pad_hinv ;
	Real const cen_dist6_pad_hinv ;
	Real const cen_dist7_pad_hinv ;
	Real const cen_dist10_pad_hinv;
	Real const cen_dist12_pad_hinv;
	bool calculated_;
	bool no_interpolate_Mpair_;
	bool Menv_penalties_;
  bool Membed_init_;
	bool memb_center_search_;
	bool memb_normal_search_;
	Size membrane_center_max_delta_;
	Size membrane_normal_start_angle_;
	Size membrane_normal_delta_angle_;
	Size membrane_normal_max_angle_;
	Size membrane_normal_cycles_;
	Real membrane_normal_magnitude_;
	Real membrane_center_magnitude_;
	Real smooth_move_frac_;
};

  MembraneEmbed const & MembraneEmbed_from_pose( pose::Pose const & pose );
  MembraneEmbed & nonconst_MembraneEmbed_from_pose( pose::Pose & pose );

} // ns scoring
} // ns core

#endif
