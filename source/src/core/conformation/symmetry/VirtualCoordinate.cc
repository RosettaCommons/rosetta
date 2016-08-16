// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief  VirtualCoordinate container
/// @file   core/conformation/symmetry/VirtualCoordinate.cc
/// @author Ingemar Andre

// Unit headers
#include <core/conformation/symmetry/VirtualCoordinate.hh>

// C++ headers
#include <iostream>

// Utility header
#include <utility/string_util.hh>

#include <utility/vector1.hh>


namespace core {
namespace conformation {
namespace symmetry {

VirtualCoordinate::VirtualCoordinate( VirtualCoordinate const & src ):
	axis_x_(src.axis_x_),
	axis_y_(src.axis_y_),
	axis_origin_(src.axis_origin_),
	mirror_Z_(src.mirror_Z_)
{}

VirtualCoordinate::VirtualCoordinate():
	axis_x_(),
	axis_y_(),
	axis_origin_(),
	mirror_Z_(false)
{}


VirtualCoordinate::VirtualCoordinate(
	numeric::xyzVector< core::Real> axis_x,
	numeric::xyzVector< core::Real> axis_y,
	numeric::xyzVector< core::Real> axis_origin
):
	axis_x_(axis_x),
	axis_y_(axis_y),
	axis_origin_(axis_origin),
	mirror_Z_(false)
{}

VirtualCoordinate::VirtualCoordinate(
	numeric::xyzVector< core::Real> axis_x,
	numeric::xyzVector< core::Real> axis_y,
	numeric::xyzVector< core::Real> axis_origin,
	bool mirror_z
):
	axis_x_(axis_x),
	axis_y_(axis_y),
	axis_origin_(axis_origin),
	mirror_Z_(mirror_z)
{}

VirtualCoordinate &
VirtualCoordinate::operator=( VirtualCoordinate const & src ) {
	axis_x_ = src.axis_x_;
	axis_y_ = src.axis_y_;
	axis_origin_ = src.axis_origin_;
	mirror_Z_ = src.mirror_Z_;
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

// @details accessor functions
numeric::xyzVector< core::Real> const &
VirtualCoordinate::get_x() const
{
	return axis_x_;
}

numeric::xyzVector< core::Real> const &
VirtualCoordinate::get_y() const
{
	return axis_y_;
}

numeric::xyzVector< core::Real> const &
VirtualCoordinate::get_origin() const
{
	return axis_origin_;
}


bool
VirtualCoordinate::get_mirror_z() const {
	return mirror_Z_;
}

void
VirtualCoordinate::set_mirror_z( bool val ) {
	mirror_Z_ = val;
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
	debug_assert( coords.size() >= 3 );
	utility::vector1< std::string> split ( utility::string_split( coords[ coord_start  ], ',' ) );
	debug_assert( split.size() == 3 );
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

bool
operator==(
	VirtualCoordinate const & a,
	VirtualCoordinate const & b
) {
	return
		(a.axis_x_ == b.axis_x_) &&
		(a.axis_y_ == b.axis_y_) &&
		(a.axis_origin_ == b.axis_origin_) &&
		(a.mirror_Z_ == b.mirror_Z_);
}

bool
operator!=(
	VirtualCoordinate const & a,
	VirtualCoordinate const & b
) {
	return !(a == b);
}


} // symmetry
} // conformation
} // core
