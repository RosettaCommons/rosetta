// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

//////////////////////////////////////////////
///
/// @file protocols/scoring/methods/pcs2/GridSearchIteratorCA.cc
///
/// @brief
///
/// @details
///
/// @param
///
/// @return
///
/// @remarks
///
/// @references
///
/// @authorv Christophe Schmitz
///
////////////////////////////////////////////////


// Unit headers
#include <protocols/scoring/methods/pcs2/GridSearchIteratorCA.hh>

// Package headers

// Project headers

// Utility headers
#include <utility/exit.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

// Objexx headers

// C++ headers
#include <iostream>

#include <utility/vector1.hh>


namespace protocols {
namespace scoring {
namespace methods {
namespace pcs2 {


GridSearchIteratorCA::GridSearchIteratorCA()
{
	utility_exit_with_message( "You shouldn't call the empty constructor for GridSearchIteratorCA class" );
}

GridSearchIteratorCA::~GridSearchIteratorCA(){
}

GridSearchIteratorCA::GridSearchIteratorCA(GridSearchIteratorCA const & other)
{
	x_vec_ = other.x_vec_;
	y_vec_ = other.y_vec_;
	z_vec_ = other.z_vec_;
	res_num_cur_ = other.res_num_cur_;
	res_num_total_ = other.res_num_total_;
}

GridSearchIteratorCA &
GridSearchIteratorCA::operator=(GridSearchIteratorCA const & other){
	if ( this != &other ) {
		x_vec_ = other.x_vec_;
		y_vec_ = other.y_vec_;
		z_vec_ = other.z_vec_;
		res_num_cur_ = other.res_num_cur_;
		res_num_total_ = other.res_num_total_;
	}
	return *this;
}

void
GridSearchIteratorCA::reset(){
	res_num_cur_ = 1;
}


void
GridSearchIteratorCA::set_vec(utility::vector1<core::Real> & x_vec,
	utility::vector1<core::Real> & y_vec,
	utility::vector1<core::Real> & z_vec,
	core::Size index,
	core::pose::Pose const & pose){

	numeric::xyzVector< core::Real > coo1 = pose.residue(index).atom("CA").xyz();
	numeric::xyzVector< core::Real > coo2 = pose.residue(index).atom("C").xyz();

	//I can get a divistion by zero if I use only CA

	x_vec.push_back((coo1.x() + coo2.x())/2.0);
	y_vec.push_back((coo1.y() + coo2.y())/2.0);
	z_vec.push_back((coo1.z() + coo2.z())/2.0);

}

GridSearchIteratorCA::GridSearchIteratorCA(core::pose::Pose const & pose)
{
	core::Size i;
	core::Size n_res(pose.size());

	res_num_cur_ = 1;


	if ( n_res <= 8 ) {
		for ( i = 1; i <= n_res; i++ ) {
			set_vec(x_vec_, y_vec_, z_vec_, i, pose);
		}
		res_num_total_ = x_vec_.size();
		return;
	}


	set_vec(x_vec_, y_vec_, z_vec_, 1, pose);
	set_vec(x_vec_, y_vec_, z_vec_, n_res, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res + 1) / 2, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 1 + 3) / 4, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 3 + 1) / 4, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 1 + 7) / 8, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 7 + 1) / 8, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 3 + 5) / 8, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 5 + 3) / 8, pose);

	/*
	if(n_res <= 16){
	res_num_total_ = x_vec_.size();
	return;
	}


	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 1 + 15) / 16, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 15 + 1) / 16, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 5 + 11) / 16, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 11+ 5) / 16, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 3 + 13) / 16, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 13 + 3) / 16, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 7 + 9) / 16, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 9 + 7) / 16, pose);


	if(n_res <= 32){
	res_num_total_ = x_vec_.size();
	return;
	}


	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 1 + 31) / 32, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 31 + 1) / 32, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 9 + 23) / 32, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 23 + 9) / 32, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 13 + 19) / 32, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 19 + 13) / 32, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 5 + 27) / 32, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 27 + 5) / 32, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 15 + 17) / 32, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 17 + 15) / 32, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 7 + 25) / 32, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 25 + 7) / 32, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 11 + 21) / 32, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 21 + 11) / 32, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 3 + 29) / 32, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 3 + 29) / 32, pose);


	if(n_res <= 64){
	res_num_total_ = x_vec_.size();
	return;
	}


	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 1 + 63) / 64, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 63 + 1) / 64, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 17 + 47) / 64, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 47 + 17) / 64, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 33 + 31) / 64, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 31 + 33) / 64, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 9 + 55) / 64, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 55 + 9) / 64, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 25 + 39) / 64, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 39 + 25) / 64, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 13 + 51) / 64, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 51 + 13) / 64, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 21 + 43) / 64, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 43 + 21) / 64, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 35 + 29) / 64, pose);
	set_vec(x_vec_, y_vec_, z_vec_, (n_res * 29 + 35) / 64, pose);
	*/


	if ( n_res <= 128 ) {
		res_num_total_ = x_vec_.size();
		return;
	}

	res_num_total_ = x_vec_.size();
	return;

}


bool
GridSearchIteratorCA::next_center(core::Real &x,
	core::Real &y,
	core::Real &z){
	if ( res_num_cur_ > res_num_total_ ) {
		return false;
	}

	x = x_vec_[res_num_cur_];
	y = y_vec_[res_num_cur_];
	z = z_vec_[res_num_cur_];

	res_num_cur_ ++;
	return true;
}

} //namespace pcs2
} //namespace methods
} //namespace scoring
} //namespace methods
