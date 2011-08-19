// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file   core/scoring/methods/Membrane_FAPotential.hh
/// @brief  Membrane FA Potential
/// @author Patrick Barth


#ifndef INCLUDED_core_scoring_Membrane_FAPotential_hh
#define INCLUDED_core_scoring_Membrane_FAPotential_hh

#include <core/types.hh>

// Unit headers
#include <core/scoring/Membrane_FAPotential.fwd.hh> //pba
#include <core/scoring/MembranePotential.hh> //pba
#include <core/scoring/MembranePotential.fwd.hh>
#include <core/scoring/MembraneTopology.fwd.hh>
#include <core/scoring/MembraneTopology.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/EnvPairPotential.hh>

// Package headers
#include <core/scoring/EnergyGraph.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <basic/datacache/CacheableData.hh>

// Utility headers
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray4D.hh>

// C++

namespace core {
namespace scoring {


class Membrane_FAEmbed : public basic::datacache::CacheableData {


public:
	Membrane_FAEmbed(): calculated_(false) {};
	Membrane_FAEmbed( Membrane_FAEmbed const & src );


	basic::datacache::CacheableDataOP
	clone() const
	{
		return new Membrane_FAEmbed( *this );
	}

  Real &
  fa_proj(Size const seqpos, Size const atom)
  {
    return fa_proj_[seqpos][atom];
  }

  Real
  fa_proj(Size const seqpos, Size const atom) const
  {
    return fa_proj_[seqpos][atom];
  }

  Real &
  fa_depth(Size const seqpos, Size const atom)
  {
    return fa_depth_[seqpos][atom];
  }

  Real
  fa_depth(Size const seqpos, Size const atom) const
  {
    return fa_depth_[seqpos][atom];
  }

  Real &
  fa_proj_deriv(Size const seqpos, Size const atom)
  {
    return fa_proj_deriv_[seqpos][atom];
  }

  Real
  fa_proj_deriv(Size const seqpos, Size const atom) const
  {
    return fa_proj_deriv_[seqpos][atom];
  }

  Vector &
  fa_proj_coord(Size const seqpos, Size const atom)
  {
    return fa_proj_coord_[seqpos][atom];
  }

  Vector
  fa_proj_coord(Size const seqpos, Size const atom) const
  {
    return fa_proj_coord_[seqpos][atom];
  }

  Real
  fa_center() const //pba
  {
    return fa_center_;
  }

  Real &
  fa_center() //pba
  {
    return fa_center_;
  }

  Real
  fa_penalty() const //pba
  {
    return fa_penalty_;
  }

  Real &
  fa_penalty() //pba
  {
    return fa_penalty_;
  }

  Real
  thickness() const //pba
  {
    return thickness_;
  }

  Real &
  thickness() //pba
  {
    return thickness_;
  }

  Real
  steepness() const //pba
  {
    return steepness_;
  }

  Real &
  steepness() //pba
  {
    return steepness_;
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

  bool
  Fa_Membed_update() const
  {
    return Fa_Membed_update_;
  }

  bool &
  Fa_Membed_update()
  {
    return Fa_Membed_update_;
  }

	void
	initialize( pose::Pose const & pose );

private:

  void
  allocate_appropriate_memory( pose::Pose const & pose ) const;

/////////////////////////////////////////////////////////////////////////////
// data
/////////////////////////////////////////////////////////////////////////////

private:
	mutable utility::vector1 < utility::vector1 < Real > > fa_proj_; //pba
  mutable utility::vector1 < utility::vector1 < Real > > fa_depth_; //pba
  mutable utility::vector1 < utility::vector1 < Vector > > fa_proj_coord_; //pba
  mutable utility::vector1 < utility::vector1 < Real > > fa_proj_deriv_; //pba
	bool calculated_;
  mutable Real fa_center_; //pba
  mutable Real fa_penalty_;
  Real thickness_; //pba
  Real steepness_; //pba
  bool Fa_Membed_update_;
};

//class Membrane_FAPotential : public MembranePotential {
class Membrane_FAPotential : public EnvPairPotential {

public:
	Membrane_FAPotential():
    membrane_potential_( ScoringManager::get_instance()->get_MembranePotential() ) {};

  void
  compute_fa_projection(pose::Pose & pose) const;
//  compute_fa_projection(pose::Pose const & pose) const;

	void
	finalize( pose::Pose & pose ) const;

protected:
/*
  Membrane_FAEmbed const & Membrane_FAEmbed_from_pose( pose::Pose const & ) const;
  Membrane_FAEmbed & nonconst_Membrane_FAEmbed_from_pose( pose::Pose & ) const;
	MembraneEmbed const & MembraneEmbed_from_pose( pose::Pose const & ) const;
	MembraneEmbed & nonconst_MembraneEmbed_from_pose( pose::Pose & ) const;
	//MembraneTopology const & MembraneTopology_from_pose( pose::Pose const & ) const;
	//MembraneTopology & nonconst_MembraneTopology_from_pose( pose::Pose & ) const;
  MembranePotential & nonconst_MembranePotential_from_pose( pose::Pose & ) const;
*/
private: // function

  void
  fa_projection(
    pose::Pose & pose,
    Vector const & normal,
    Vector const & center,
    Real const & thickness,
    Real const & steepness,
    Real const & penalty
  ) const;

private: // data

//	ObjexxFCL::FArray3D< Real > mem_env_log6_;
//	ObjexxFCL::FArray1D< Real > cenpack_log_;
	bool calculated_; //pba needed ??
  MembranePotential const & membrane_potential_;

};

  Membrane_FAEmbed const & Membrane_FAEmbed_from_pose( pose::Pose const & );
  Membrane_FAEmbed & nonconst_Membrane_FAEmbed_from_pose( pose::Pose & );


} // ns scoring
} // ns core

#endif
