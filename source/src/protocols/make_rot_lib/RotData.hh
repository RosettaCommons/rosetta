// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/make_rot_lib/MakeRotLibMover.hh
/// @brief Header file for RotData
/// @author P. Douglas Renfrew ( renfrew@nyu.edu )

#ifndef INCLUDED_protocols_make_rot_lib_rotdata_hh
#define INCLUDED_protocols_make_rot_lib_rotdata_hh

// unit headers
#include <protocols/make_rot_lib/RotData.fwd.hh>

// core headers
#include <core/types.hh>

// utility headers
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>

// numeric headers
#include <numeric/angle.functions.hh>

// C++ headers
#include <iostream>

namespace protocols {
namespace make_rot_lib {

// Define Rotamer Data class and variable
class RotData
{
private:
	utility::vector1<core::Real> bbs_;
	utility::vector1<core::Size> bb_ids_;
	core::Size num_bbs_;
	core::Real omega_, min_omega_, epsilon_, min_epsilon_;
	core::Real energy_;
	core::Real probability_;
	core::Size num_chi_;						// number of chi angles in AA
	core::Size num_clusters_;			 // number of clusters
	core::Size cluster_num_;				// cluster id
	
	utility::vector1< core::Real > semi_energy_dist_;
	utility::vector1< core::Real > semi_prob_dist_;
	
	utility::vector1< core::Real > inp_chi_;	// starting chi angles
	utility::vector1< core::Real > min_chi_;	// minimized chi angles
	utility::vector1< core::Size > lib_chi_val_; // rotamer number for dunbrack format
	utility::vector1< core::Real > std_dev_;	// standard deviation of chi angles
	utility::vector1< core::Real > cen_dst_;	// distance from each centroid

	// for debuging
	core::Real twist_;
	core::Real inter_rep_;
	core::Real inter_atr_;
	core::Real intra_rep_;
	core::Real intra_atr_;
	core::Real solvation_;
	
	bool semirotameric_;

public:
	
	void set_semi_energy_dist( core::Size i, core::Real setting ) {
		semi_energy_dist_[ i ] = setting;
	}
	
	void set_semi_prob_dist( core::Size i, core::Real setting ) {
		semi_prob_dist_[ i ] = setting;
	}
	
	void resize_semi_vectors( core::Size i ) {
		semi_energy_dist_.resize( i, 0 );
		semi_prob_dist_.resize( i, 0 );
	}
	
	void set_twist( core::Real twist ) {
		twist_ = twist;
	}
	void set_inter_rep( core::Real inter_rep ) {
		inter_rep_ = inter_rep;
	}
	void set_inter_atr( core::Real inter_atr ) {
		inter_atr_ = inter_atr;
	}
	void set_intra_rep( core::Real intra_rep ) {
		intra_rep_ = intra_rep;
	}
	void set_intra_atr( core::Real intra_atr ) {
		intra_atr_ = intra_atr;
	}
	void set_solvation( core::Real solvation ) {
		solvation_ = solvation;
	}

	bool get_semirotameric() {
		return semirotameric_;
	}
	
	core::Real get_semi_prob_dist( core::Size i ) {
		return semi_prob_dist_[ i ];
	}
	
	core::Real get_semi_energy_dist( core::Size i ) {
		return semi_energy_dist_[ i ];
	}
	
	core::Real get_twist() {
		return twist_;
	}
	core::Real get_inter_rep() {
		return inter_rep_;
	}
	core::Real get_inter_atr() {
		return inter_atr_;
	}
	core::Real get_intra_rep() {
		return intra_rep_;
	}
	core::Real get_intra_atr() {
		return intra_atr_;
	}
	core::Real get_solvation() {
		return solvation_;
	}

 public:
	// ctor
	RotData( core::Size NumChi, core::Size NumCluster );
	RotData( core::Size NumChi, core::Size NumCluster, bool semirotameric );
	RotData( core::Size NumChi, core::Size NumBBs, core::Size NumCluster );
	RotData( core::Size NumChi, core::Size NumBBs, core::Size NumCluster, bool semirotameric );

	friend bool operator== ( RotData & r1, RotData & r2 );
	
	// setters and getters
	void set_bb( core::Size i, core::Real BB ) {
		bbs_[ i ] = BB;
	}
    
    // setters and getters
	void set_bb_id( core::Size i, core::Size bbid ) {
		bb_ids_[ i ] = bbid;
	}
    
    /* deprecated - implemented for the MakeRotLib unit test */
	void set_phi (core::Real Phi) {
		bbs_[ 1 ] = Phi;
    }
    
	void set_psi (core::Real Psi) {
		bbs_[ 2 ] = Psi;
    }
	///////////////////////////////////////////////////////////
    
    void set_omg (core::Real Omega) {
		omega_ = Omega;
	}
	void set_min_omg (core::Real MinOmega) {
		min_omega_ = MinOmega;
	}
	void set_eps (core::Real Epsilon) {
		epsilon_ = Epsilon;
	}
	void set_min_eps (core::Real MinEpsilon) {
		min_epsilon_ = MinEpsilon;
	}
    
	void set_num_bbs(core::Size i) {
		num_bbs_ = i;
	}
	core::Size get_num_bbs() {
		return num_bbs_;
	}
	utility::vector1< core::Real > get_bbs() {
		return bbs_;
	}
	void resize_bbs( core::Size i ) {
		bbs_.resize( i );
	}
	void resize_bb_ids( core::Size i ) {
		bb_ids_.resize( i );
	}
	core::Real get_bb( core::Size i ) {
		return bbs_[ i ];
	}
	core::Size get_bb_id( core::Size i ) {
		return bb_ids_[ i ];
	}
    
    /* deprecated - implemented for the MakeRotLib unit test */
	core::Real get_phi() {
		return bbs_[1];
	}
    
	core::Real get_psi() {
		return bbs_[2];
	}
    ///////////////////////////////////////////////////////////
    
    
	core::Real get_omg() {
		return omega_;
	}
	core::Real get_min_omg() {
		return min_omega_;
	}
	core::Real get_eps() {
		return epsilon_;
	}
	core::Real get_min_eps() {
		return min_epsilon_;
	}

	void set_energy (core::Real Energy) {
		energy_ = Energy;
	}

	core::Real get_energy () {
		return energy_;
	}

	void set_probability (core::Real Probability) {
		probability_ = Probability;
	}

	core::Real get_probability () {
		return probability_;
	}

	void set_num_chi(core::Size Num_Chi) {
		num_chi_ = Num_Chi;
	}

	core::Size get_num_chi() {
		return num_chi_;
	}

	void set_num_clusters( core::Size num ) {
		num_clusters_ = num;
	}

	core::Size get_num_clusters() {
		return num_clusters_;
	}

	void set_cluster_num (core::Size Cluster_Num) {
		cluster_num_ = Cluster_Num;
	}

	core::Size get_cluster_num () {
		return cluster_num_;
	}

	void set_inp_chi( core::Real angle, core::Size num ) {
		inp_chi_[ num ] =	numeric::nonnegative_principal_angle_degrees( angle );
	}

	core::Real get_inp_chi( core::Size num ) {
		return inp_chi_[ num ];
	}

	void set_min_chi ( core::Real angle, core::Size num) {
		min_chi_[ num ] =	numeric::nonnegative_principal_angle_degrees( angle );
	}

	core::Real get_min_chi( core::Size num ) {
		return min_chi_[ num ];
	}

	void set_lib_chi_val ( int val, core::Size num) {
		lib_chi_val_[ num ] =	val;
	}

	core::Real get_lib_chi_val( core::Size num ) {
		return lib_chi_val_[ num ];
	}

	void set_std_dev (core::Real STD, core::Size num) {
		std_dev_[ num ] = STD;
	}

	core::Real get_std_dev( core::Size num ) {
		return std_dev_[ num ];
	}

	void set_cen_dist(core::Real dist , core::Size num) {
		cen_dst_[ num ]= dist;
	}

	core::Real get_cen_dist( core::Size num ) {
		return cen_dst_[ num ];
	}

	core::Size get_min_cent_dist(){
		return arg_min( cen_dst_ );	//gets index of closest centroid
	}

	void show( std::ostream & out ) const;
	bool load( std::istream & in );
};

} // namespace make_rot_lib
} // namespace protocols

#endif // INCLUDED_protocols_makerotlib_rotdata_HH
