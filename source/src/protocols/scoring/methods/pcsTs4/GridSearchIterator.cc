// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

 //////////////////////////////////////////////
 ///
 /// @file GridSearchIterator.cc
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
 /// @references C Schmitz et.al. J Mol Biol. Mar 9, 2012; 416(5): 668–677 ; Yagi H et.al Structure, 2013, 21(6):883-890
 ///
 /// @authorv Christophe Schmitz , Kala Bharath Pilla
 /// ,
 ////////////////////////////////////////////////


// Unit headers
#include <protocols/scoring/methods/pcsTs4/GridSearchIterator.hh>

// Package headers

// Project headers

// Utility headers
#include <utility/exit.hh>

// Numeric headers
#include <numeric/constants.hh>
#include <numeric/xyzVector.hh>

// Objexx headers

// C++ headers
#include <iostream>

//#include <limits>

namespace protocols{
namespace scoring{
namespace methods{
namespace pcsTs4{


GridSearchIterator_Ts4::GridSearchIterator_Ts4():
	x_center_(0),
	y_center_(0),
	z_center_(0),
	step_(0),
	edge_(0),
	delta_(0),
	small_cutoff_square_(0),
	large_cutoff_square_(0),
	x_vector_(0),
	y_vector_(0),
	z_vector_(0),
	norme_vector_(0),
	cone_angle_cos_(0)
{
	utility_exit_with_message( "You shouldn't call the empty constructor for GridSearchIterator_Ts4 class" );
}

GridSearchIterator_Ts4::~GridSearchIterator_Ts4(){
}

GridSearchIterator_Ts4::GridSearchIterator_Ts4(GridSearchIterator_Ts4 const & other):
	ReferenceCount(),
	x_center_(other.x_center_),
	y_center_(other.y_center_),
	z_center_(other.z_center_),
	step_(other.step_),
	edge_(other.edge_),
	delta_(other.delta_),
	small_cutoff_square_(other.small_cutoff_square_),
	large_cutoff_square_(other.large_cutoff_square_),
	x_vector_(other.x_vector_),
	y_vector_(other.y_vector_),
	z_vector_(other.z_vector_),
	norme_vector_(other.norme_vector_),
	cone_angle_cos_(other.cone_angle_cos_)
{
	x_current_ = other.x_current_;
	y_current_ = other.y_current_;
	z_current_ = other.z_current_;
	step_x_ = other.step_x_;
	step_y_ = other.step_y_;
	step_z_ = other.step_z_;
	next_to_return_ = other.next_to_return_;
}

GridSearchIterator_Ts4 &
GridSearchIterator_Ts4::operator=(GridSearchIterator_Ts4 const & other){
	if ( this != &other ) {
		if((x_center_ != other.x_center_)||
			 (y_center_ != other.y_center_)||
			 (z_center_ != other.z_center_)||
			 (step_ != other.step_)||
			 (edge_ != other.edge_)||
			 (delta_ != other.delta_)){
			utility_exit_with_message( "You can't call the = operator on GridSearchIterator_Ts4 object of different size" );
		}
		x_current_ = other.x_current_;
		y_current_ = other.y_current_;
		z_current_ = other.z_current_;
		step_x_ = other.step_x_;
		step_y_ = other.step_y_;
		step_z_ = other.step_z_;
		next_to_return_ = other.next_to_return_;
	}
	return *this;
}


GridSearchIterator_Ts4::GridSearchIterator_Ts4(numeric::xyzVector< core::Real > const coo1,
	numeric::xyzVector< core::Real > const coo2,
	core::Real const k,
	core::Real const edge_size,
	core::Real const step_size,
	core::Real const small_cutoff,
	core::Real const large_cutoff,
	core::Real const cone_angle
):
	x_center_(coo2.x() + k*( coo2.x() - coo1.x() )),
	y_center_(coo2.y() + k*( coo2.y() - coo1.y() )),
	z_center_(coo2.z() + k*( coo2.z() - coo1.z() )),
	step_(step_size),
	edge_(edge_size),
	delta_(step_/2.0),
	small_cutoff_square_(small_cutoff * small_cutoff),
	large_cutoff_square_(large_cutoff * large_cutoff),
	x_vector_(coo1.x() - coo2.x()),
	y_vector_(coo1.y() - coo2.y()),
	z_vector_(coo1.z() - coo2.z()),
	norme_vector_(sqrt( (coo2.x()-coo1.x())*(coo2.x()-coo1.x()) + (coo2.y()-coo1.y())*(coo2.y()-coo1.y()) + (coo2.z()-coo1.z())*(coo2.z()-coo1.z()) )),
	cone_angle_cos_(cos(cone_angle/180.0*core::Real( numeric::constants::d::pi )))
{
	if( edge_size < 0){
		utility_exit_with_message("Edge size of the cube search is negative and has to be positive");
	}
	if( step_size <= 0){
		utility_exit_with_message("step_size of the cube search has to be strictly positive");
	}

	reset();
}


bool
GridSearchIterator_Ts4::next_center(core::Real &x,
																core::Real &y,
																core::Real &z){

	//double r2;
	while(next(x, y, z) == true){
		double r2 = (x-x_center_)*(x-x_center_) + (y-y_center_)* (y-y_center_) + (z-z_center_)*(z-z_center_);
		if((r2 <= large_cutoff_square_) && (r2 >= small_cutoff_square_ )){ //we test small and large cutoff
			core::Real cos_angle = ((x-x_center_)*x_vector_ + (y-y_center_)*y_vector_ + (z-z_center_)*z_vector_ );
			if((r2 == 0 ) || (norme_vector_ == 0)){
				return true;
			}
			cos_angle = cos_angle / sqrt(r2) / norme_vector_;
			if(cos_angle >= cone_angle_cos_){ //we test cone angle
				return true;
			}
		}
	}
	return false;
}

// This iterator sample a cube in 3D space.
// The following code looks strange and could be written in an easier way.
// However, it is strange because the iterator ensure that
// the points are visited in a specific order in such a way
// that the next point visited is a direct neighboor of the previous one
// It will be important if I decide to work on unassigned PCS Data.
bool
GridSearchIterator_Ts4::next(core::Real &x,
												 core::Real &y,
												 core::Real &z){
	x = x_current_;
	y = y_current_;
	z = z_current_;

	if(!next_to_return_){
		reset();
		return(false);
	}

	if( z_current_ - z_center_ > edge_/2.0){
		next_to_return_ = false;
		return(true);
	}

	x_current_ += step_x_;

	if( x_current_ - x_center_> edge_/2.0 + delta_){
			step_x_ = -step_;
			x_current_ += step_x_;
			y_current_ += step_y_;
	}

	if( x_current_ - x_center_ < -edge_/2.0 - delta_){
			step_x_ = step_;
			x_current_ += step_x_;
			y_current_ += step_y_;
	}

	if( y_current_ - y_center_> edge_/2.0 + delta_){
			step_y_ = -step_;
			y_current_ += step_y_;
			z_current_ += step_z_;
	}

	if( y_current_ - y_center_ < -edge_/2.0 - delta_){
			step_y_ = step_;
			y_current_ += step_y_;
			z_current_ += step_z_;
	}

	if( z_current_ - z_center_ > edge_/2.0 + delta_){
		next_to_return_ = false;
		return(true);
	}

	return(true);
}

void
GridSearchIterator_Ts4::reset(){
		x_current_ = -edge_/2.0 + x_center_;
		y_current_ = -edge_/2.0 + y_center_;
		z_current_ = -edge_/2.0 + z_center_;
		step_x_ = step_;
		step_y_ = step_;
		step_z_ = step_;
		next_to_return_ = true;
}

} //namespace pcsTs4
} //namespace methods
} //namespace scoring
} //namespace methods
