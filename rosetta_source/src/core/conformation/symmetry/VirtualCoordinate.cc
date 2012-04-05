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

/// @brief  VirtualCoordinate container
/// @file   core/conformation/symmetry/VirtualCoordinate.cc
/// @author Ingemar Andre

// Unit headers
#include <core/conformation/symmetry/VirtualCoordinate.hh>

// C++ headers
#include <iostream>

// Utility header
// AUTO-REMOVED #include <utility/io/izstream.hh>
#include <utility/string_util.hh>

#include <utility/vector1.hh>


namespace core {
namespace conformation {
namespace symmetry {

VirtualCoordinate::VirtualCoordinate( VirtualCoordinate const & src )
{
    axis_x_ = src.axis_x_;
    axis_y_ = src.axis_y_;
    axis_origin_ = src.axis_origin_;
}

VirtualCoordinate::VirtualCoordinate(){}

	/// @brief copy constructor
//VirtualCoordinate::VirtualCoordinate( VirtualCoordinate const & src );

VirtualCoordinate::VirtualCoordinate(
		numeric::xyzVector< core::Real> axis_x,
		numeric::xyzVector< core::Real> axis_y,
		numeric::xyzVector< core::Real> axis_origin
	)
	{
		axis_x_ = axis_x;
		axis_y_ = axis_y;
		axis_origin_ = axis_origin;
	}

	VirtualCoordinate &
  VirtualCoordinate::operator=( VirtualCoordinate const & src ) {
		axis_x_ = src.axis_x_;
		axis_y_ = src.axis_y_;
		axis_origin_ = src.axis_origin_;
		return *this;
	}

	VirtualCoordinate::~VirtualCoordinate(){}

	// @details accessor functions
	numeric::xyzVector< core::Real> &
	VirtualCoordinate::get_x()
	{
		return axis_x_;
	}

	numeric::xyzVector< core::Real> &
	VirtualCoordinate::get_y()
	{
		return axis_y_;
	}

	numeric::xyzVector< core::Real> &
	VirtualCoordinate::get_origin()
	{
		return axis_origin_;
	}
// @details read the coordinates of a virtual residues from string. Start reading
// coordinates from coord_start. The coordinates correspond to the unit vectors for
// X, Y axis and a origin. Vectors are not automatically normalized here. Should we
// do that?
void
VirtualCoordinate::add_coordinate_from_string(
										utility::vector1< std::string > coords,
                    core::Size coord_start )
{
		assert( coords.size() >= 3 );
		utility::vector1< std::string> split ( utility::string_split( coords[ coord_start  ], ',' ) );
		assert( split.size() == 3 );
    axis_x_ = Vector( ( static_cast<core::Real>( std::atof( split[1].c_str() ) ) ),
                                             ( static_cast<core::Real>( std::atof( split[2].c_str() ) ) ),
                                             ( static_cast<core::Real>( std::atof( split[3].c_str() ) ) ) );
		split = utility::string_split( coords[ coord_start +1 ], ',' );
		axis_y_ = Vector( ( static_cast<core::Real>( std::atof( split[1].c_str() ) ) ),
                                             ( static_cast<core::Real>( std::atof( split[2].c_str() ) ) ),
                                             ( static_cast<core::Real>( std::atof( split[3].c_str() ) ) ) );
		axis_origin_ = Vector(0,0,0);
		if ( coords.size() == 5 ) {
			split = utility::string_split( coords[ coord_start +2 ], ',' );
			axis_origin_ = Vector( ( static_cast<core::Real>( std::atof( split[1].c_str() ) ) ),
                                             ( static_cast<core::Real>( std::atof( split[2].c_str() ) ) ),
                                             ( static_cast<core::Real>( std::atof( split[3].c_str() ) ) ) );
		}
}

} // symmetry
} // conformation
} // core
