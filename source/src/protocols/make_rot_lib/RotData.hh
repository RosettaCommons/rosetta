// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2008 University of Washington
// (C) 199x-2008 University of California Santa Cruz
// (C) 199x-2008 University of California San Francisco
// (C) 199x-2008 Johns Hopkins University
// (C) 199x-2008 University of North Carolina, Chapel Hill
// (C) 199x-2008 Vanderbilt University

#ifndef INCLUDED_protocols_make_rot_lib_rotdata_HH
#define INCLUDED_protocols_make_rot_lib_rotdata_HH

// core headers
#include <core/types.hh>

// utility headers
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh> //for arg_min(vector1<>)

// numeric headers
#include <numeric/angle.functions.hh>

// C++ headers
#include <iostream>

namespace protocols {
namespace MakeRotLib {

using namespace core;
using namespace utility;

// Define Rotamer Data class and variable
class RotData {
 private:
  Real phi_, psi_, omega_, min_omega_, epsilon_, min_epsilon_;
  Real energy_;
  Real probability_;
  Size num_chi_;            // number of chi angles in AA
	Size num_clusters_;       // number of clusters
  Size cluster_num_;        // cluster id
  vector1<Real> inp_chi_;  // starting chi angles
  vector1<Real> min_chi_;  // minimized chi angles
	vector1<int> lib_chi_val_; // rotamer number for dunbrack format
  vector1<Real> std_dev_;  // standard deviation of chi angles
  vector1<Real> cen_dst_;  // distance from each centroid

	// for debuging
	Real twist_;
	Real inter_rep_;
	Real inter_atr_;
	Real intra_rep_;
	Real intra_atr_;
	Real solvation_;

public:
	void set_twist( Real twist ) {
		twist_ = twist;
	}
	void set_inter_rep( Real inter_rep ) {
		inter_rep_ = inter_rep;
	}
	void set_inter_atr( Real inter_atr ) {
		inter_atr_ = inter_atr;
	}
	void set_intra_rep( Real intra_rep ) {
		intra_rep_ = intra_rep;
	}
	void set_intra_atr( Real intra_atr ) {
		intra_atr_ = intra_atr;
	}
	void set_solvation( Real solvation ) {
		solvation_ = solvation;
	}

	Real get_twist() {
		return twist_;
	}
	Real get_inter_rep() {
		return inter_rep_;
	}
	Real get_inter_atr() {
		return inter_atr_;
	}
	Real get_intra_rep() {
		return intra_rep_;
	}
	Real get_intra_atr() {
		return intra_atr_;
	}
	Real get_solvation() {
		return solvation_;
	}

 public:
  // ctor
  RotData( Size NumChi, Size NumCluster );

  // setters and getters
  void set_phi (Real Phi) {
    phi_ = Phi;
  }
  void set_psi (Real Psi) {
    psi_ = Psi;
  }
  void set_omega (Real Omega) {
    omega_ = Omega;
  }
	void set_min_omega (Real MinOmega) {
    min_omega_ = MinOmega;
  }
  void set_epsilon (Real Epsilon) {
    epsilon_ = Epsilon;
  }
  void set_min_epsilon (Real MinEpsilon) {
    min_epsilon_ = MinEpsilon;
  }
  Real get_phi() {
    return phi_;
  }
  Real get_psi() {
    return psi_;
  }
	Real get_omega() {
		return omega_;
	}
	Real get_min_omega() {
		return min_omega_;
	}
	Real get_epsilon() {
		return epsilon_;
	}
	Real get_min_epsilon() {
		return min_epsilon_;
	}

  void set_energy (Real Energy) {
    energy_ = Energy;
  }

  Real get_energy () {
    return energy_;
  }

  void set_probability (Real Probability) {
    probability_ = Probability;
  }

  Real get_probability () {
    return probability_;
  }

  void set_num_chi(Size Num_Chi) {
    num_chi_ = Num_Chi;
  }

  Size get_num_chi() {
    return num_chi_;
  }

	void set_num_clusters( Size num ) {
		num_clusters_ = num;
	}

	Size get_num_clusters() {
		return num_clusters_;
	}

  void set_cluster_num (Size Cluster_Num) {
    cluster_num_ = Cluster_Num;
  }

  Size get_cluster_num () {
    return cluster_num_;
  }

  void set_inp_chi( Real angle, Size num ) {
    inp_chi_[ num ] =  numeric::nonnegative_principal_angle_degrees( angle );
  }

  Real get_inp_chi( Size num ) {
    return inp_chi_[ num ];
  }

  void set_min_chi ( Real angle, Size num) {
    min_chi_[ num ] =  numeric::nonnegative_principal_angle_degrees( angle );
  }

  Real get_min_chi( Size num ) {
    return min_chi_[ num ];
  }

  void set_lib_chi_val ( int val, Size num) {
    lib_chi_val_[ num ] =  val;
  }

  Real get_lib_chi_val( Size num ) {
    return lib_chi_val_[ num ];
  }

  void set_std_dev (Real STD, Size num) {
    std_dev_[ num ] = STD;
  }

  Real get_std_dev( Size num ) {
    return std_dev_[ num ];
  }

  void set_cen_dist(Real dist , Size num) {
    cen_dst_[ num ]= dist;
  }

  Real get_cen_dist( Size num ) {
    return cen_dst_[ num ];
  }

	Size get_min_cent_dist(){
		return arg_min( cen_dst_ );  //gets index of closest centroid
	}

	void show( std::ostream & out ) const;
	bool load( std::istream & in );
};

} // namespace MakeRotLib
} // namespace protocols

#endif // INCLUDED_protocols_makerotlib_rotdata_HH
