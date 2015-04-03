// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/GoapEnergy.cc
/// @brief  C++ implementaion of GOAP(Generalized Orientation-dependent, All-atom statistical Potential)
///         by Zhou H & Skolnick J, Biophys J 2011, 101(8):2043-52.
/// @author Hahnbeom Park

#ifndef INCLUDED_core_scoring_methods_GoapEnergy_HH
#define INCLUDED_core_scoring_methods_GoapEnergy_HH

// Unit Headers
#include <core/scoring/methods/GoapEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <cmath>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray5D.hh>
#include <numeric/xyz.functions.hh>
#include <utility/vector1.hh>

#include <stdio.h>
#include <iostream>

namespace core{
namespace scoring{
namespace methods{


////////////////////////////////////////////////////////
class GoapRsdType : public utility::pointer::ReferenceCount
{
public:

  GoapRsdType();

  ~GoapRsdType();

  void setup_rsdtype( chemical::ResidueTypeCOP rsd );
  void setup_connectivity( chemical::ResidueType const &rsd );

	// Set
	void set_atmid( Size const i, Size const j ){ atmid_[i] = j; }
	void add_atmname_using( std::string const str ){ atmname_using_.push_back( str ); }
	void set_root_atom( Size const i, Size const j ){ root_atom_[i] = j; }
	void set_branch_atom( Size const i, Size const j ){ branch_atom_[i] = j; }
	void set_angle_atom( Size const i, Size const j ){ angle_atom_[i] = j; }
	void set_using( Size const i, bool const val ){ is_using_[i] = val; }

  // Accessors

	Size natom() const { return natom_;}
	Size nusing() const { return atmname_using_.size(); }
  std::string atmname_using( Size const i ) const { return atmname_using_[i]; }
  bool is_using( Size const i ) const { return is_using_[i]; }
  bool connected_by_twobonds( Size const i ) const { return connected_by_twobonds_[i]; }
  Size i2( Size const i ) const { return i2_[i]; }
  Size i3( Size const i ) const { return i3_[i]; }
	Size atmid( Size const i ) const { return atmid_[i]; }
	std::string name() const { return name_; }

private:
	Size natom_;
  utility::vector1< std::string > atmname_using_;
  utility::vector1< bool > is_using_;
  utility::vector1< bool > connected_by_twobonds_;
  utility::vector1< Size > i2_; // atom index for angle definition
  utility::vector1< Size > i3_; // atom index for angle definition
	utility::vector1< Size > atmid_;
	std::map< Size, Size > root_atom_;
	std::map< Size, Size > branch_atom_;
	std::map< Size, Size > angle_atom_;
	std::string name_;

}; // GoapRsdType

////////////////////////////////////////////////////////
class GoapEnergy : public ContextIndependentTwoBodyEnergy
{
public:
  typedef ContextIndependentTwoBodyEnergy  parent;

  GoapEnergy( EnergyMethodOptions const & options );
  GoapEnergy( GoapEnergy const & src );

  ~GoapEnergy();

  /// clone
  virtual
  EnergyMethodOP
  clone() const { return EnergyMethodOP( new GoapEnergy(*this) ); }

  virtual
  void
  setup_for_packing(
     pose::Pose & ,
     utility::vector1< bool > const & ,
     utility::vector1< bool > const &
     ) const;

  virtual
  void
  setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

  virtual
  void
  setup_for_derivatives( pose::Pose & ,
												 ScoreFunction const & ) const;
  
  virtual
  void
  residue_pair_energy(
     conformation::Residue const & rsd1,
     conformation::Residue const & rsd2,
     pose::Pose const &, //pose,
     ScoreFunction const &,
     EnergyMap & emap
     ) const;

  virtual
  bool
  defines_intrares_energy( EnergyMap const & /*weights*/ ) const { return false; }

  virtual
  bool
  minimize_in_whole_structure_context( pose::Pose const & ) const { return false; }

  virtual
  Distance
  atomic_interaction_cutoff() const { return max_dis(); }

  virtual
  void indicate_required_context_graphs( utility::vector1< bool > & ) const {};

  virtual
  void
  eval_intrares_energy(
	 conformation::Residue const & ,//rsd,
    pose::Pose const & ,//pose,
    ScoreFunction const & ,//sfxn,
    EnergyMap & //emap
   ) const {} // Just do nothing - no intrares 

	virtual
	void
	eval_residue_pair_derivatives(
		conformation::Residue const & ,//rsd1,
		conformation::Residue const & ,//rsd2,
		ResSingleMinimizationData const &,
		ResSingleMinimizationData const &,
		ResPairMinimizationData const & ,//min_data,
		pose::Pose const &,
		EnergyMap const & ,//weights,
		utility::vector1< DerivVectorPair > & ,//r1_atom_derivs,
		utility::vector1< DerivVectorPair > & //r2_atom_derivs
  ) const;

  virtual
  Size version() const { return 1; }

  // Accessors
  Vector xn( Size const resno, Size const atmno ) const { return xn_[resno][atmno]; }
  Vector xd( Size const resno, Size const atmno ) const { return xd_[resno][atmno]; }
  Real max_dis() const { return max_dis_; }
	bool continuous() const { return continuous_; }
	bool eval_res( Size const resno ) const { return eval_res_[resno]; }

	Size distbin_map( Size const i ) const { 
		std::map< Size const, Size >::const_iterator it;
		it = distbin_map_.find(i);
		if( it == distbin_map_.end() ){
			return 0;
		} else {
			return it->second;
		}
	}

private:


  void set_default();

	void
	read_Goap_parameters();

	void
	read_angle_definitions( std::string const connection_file );

	void
	read_potential_values( std::string const distance_file,
												 std::string const angle_file );


	bool
  calculate_dipoles( pose::Pose const &pose, 
		     conformation::Residue const &rsd1,
		     GoapRsdTypeCOP rsdtype,
		     Size const atm1,
		     Vector &xn, 
		     Vector &xd
		     ) const;

  Real 
  get_angle_score( Real const dist,
		   Size const atype1,
		   Size const atype2,
		   Vector const xn1,
		   Vector const xd1,
		   Vector const xn2,
		   Vector const xd2,
		   Vector const xyz1,
		   Vector const xyz2
		   ) const;

  Real
  get_distance_score( Real const dist,
		      Size const atype1,
		      Size const atype2 ) const;
  

  // inline functions
  inline Size distbin( Real const dis ){ return distbin_map_[ (Size)(dis*2) ]; }

  inline Real
  calc_cosineang( Vector const v1, Vector const v2 ) const
  { 
    Real cosang = v1.dot(v2);
    cosang /= std::sqrt( v1.dot(v1) * v2.dot(v2) );

    return cosang;
  }

  inline Real
  calc_phi(Vector const v1,
	   Vector const v2,
	   Vector const v3) const 
  {
    Real const rad2deg( 57.296 );

    Real cos1 = v1.dot(v3);
    cos1 /= std::sqrt( v1.dot(v1) * v3.dot(v3) );

    Vector v4 = v3 - cos1*v1;
    Real cos2 = v4.dot( v2 ) / std::sqrt( v4.dot(v4) );

		Vector v5 = v1.cross( v2 );
    Real cos3 = v5.dot( v4 ); 

    if( cos2 <= -1.0) cos2 = -0.9999;
    if( cos2 >= 1.0)  cos2 = 0.9999;

    Real phi = std::acos( cos2 )*rad2deg;
    if( cos3 < 0 ) phi = -phi;

    return phi;
  }

private:

	mutable GoapRsdTypeMap rsdtypemap_;
  mutable utility::vector1< utility::vector1< Vector > > xn_;
  mutable utility::vector1< utility::vector1< Vector > > xd_;
	mutable utility::vector1< bool > eval_res_;

	ObjexxFCL::FArray3D< Real > distance_table_; // 167, 167, 20
	ObjexxFCL::FArray5D< int > angle_table_; // 167, 167, 20, 12, 5: store as integer by mulitplying 1e4
  
  std::map< Size const, Size > distbin_map_;

  Size N_ATMTYPES;
  Size N_DISTANCE_BINS;
  Size N_ANGLE1_BINS;
  Size N_ANGLE2_BINS;
	Size MIN_SEQ_SEPARATIONS;

  Real max_dis_;

	bool continuous_;

}; // class GoapScore

} // namespace methods
} // namespace scoring
} // namespace core

#endif
