// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief  Symmetry data container
/// @file   core/conformation/symmetry/SymmData.hh
/// @author Ingemar Andre


#ifndef INCLUDED_core_conformation_symmetry_VirtualCoordinates_hh
#define INCLUDED_core_conformation_symmetry_VirtualCoordinates_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <core/types.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

// C++ headers
#include <string>
#include <vector>

namespace core {
namespace conformation {
namespace symmetry {

class VirtualCoordinates
{

	public:

	VirtualCoordinates(){}

	/// @brief copy constructor
 	VirtualCoordinates( VirtualCoordinates const & src );

//	{
//		axis_x_ = src.axis_x_;
//		axis_y_ = src.axis_y_;
//		axis_origin_ = src.axis_origin_;
//	}

	VirtualCoordinates &
  operator=( VirtualCoordinates const & src ) {
		axis_x_ = src.axis_x_;
		axis_y_ = src.axis_y_;
		axis_origin_ = src.axis_origin_;
		return *this;
	}

	~VirtualCoordinates(){}

	void
	push_back_x( numeric::xyzVector< core::Real> & vector )
	{
		axis_x_.push_back( vector );
	}

	void
	push_back_y( numeric::xyzVector< core::Real> & vector )
  {
    axis_y_.push_back( vector );
  }

	void
	push_back_origin( numeric::xyzVector< core::Real> & vector )
  {
    axis_origin_.push_back( vector );
  }

	void
	push_back( numeric::xyzVector< core::Real> & x,
						 numeric::xyzVector< core::Real> & y,
						 numeric::xyzVector< core::Real> & origin )
	{
		push_back_x( x );
		push_back_y( y );
		push_back_origin( origin );
	}

	numeric::xyzVector< core::Real> &
	get_x( core::Size num )
	{
	debug_assert( axis_y_.size() >= num );
  	return axis_x_[num-1];

	}

	numeric::xyzVector< core::Real> &
	get_y( core::Size num )
	{
	debug_assert( axis_y_.size() >= num );
  	return axis_y_[num-1];
	}

	numeric::xyzVector< core::Real> &
	get_origin( core::Size num )
	{
	 debug_assert( axis_y_.size() >= num );
  	return axis_origin_[num-1];
	}

	bool
	kosher()
	{
		if ( axis_x_.size() == axis_y_.size() && axis_x_.size() == axis_origin_.size() ) return true;
		else return false;
	}

	Size
	size()
	{
	debug_assert( kosher() );
		return axis_x_.size();
	}


	void
	add_coordinate_from_string( std::vector< std::string > coords,
															core::Size coord_start=2 );


	private:

	std::vector < numeric::xyzVector< core::Real> > axis_x_;
	std::vector < numeric::xyzVector< core::Real> > axis_y_;
	std::vector < numeric::xyzVector< core::Real> > axis_origin_;
};

} // symmetry
} // conformation
} // core
#endif
