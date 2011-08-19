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
/// @file   core/conformation/symmetry/VirtualCoordinates.cc
/// @author Ingemar Andre

// Unit headers
#include <core/conformation/symmetry/VirtualCoordinates.hh>

// C++ headers
#include <iostream>

// Utility header
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>

namespace core {
namespace conformation {
namespace symmetry {

VirtualCoordinates::VirtualCoordinates( VirtualCoordinates const & src )
{
    axis_x_ = src.axis_x_;
    axis_y_ = src.axis_y_;
    axis_origin_ = src.axis_origin_;
}

void
VirtualCoordinates::add_coordinate_from_string(
										std::vector< std::string > coords,
                    core::Size coord_start )
{
		assert( coords.size() >= 3 );
		std::vector< std::string> split ( utility::string_split( coords[ coord_start -1 ], ',' ) );
		assert( split.size() == 3 );
    Vector x( ( static_cast<core::Real>( std::atof( split[0].c_str() ) ) ),
                                             ( static_cast<core::Real>( std::atof( split[1].c_str() ) ) ),
                                             ( static_cast<core::Real>( std::atof( split[2].c_str() ) ) ) );
		split = utility::string_split( coords[ coord_start ], ',' );
		Vector y( ( static_cast<core::Real>( std::atof( split[0].c_str() ) ) ),
                                             ( static_cast<core::Real>( std::atof( split[1].c_str() ) ) ),
                                             ( static_cast<core::Real>( std::atof( split[2].c_str() ) ) ) );
		Vector origin(0,0,0);
		if ( coords.size() == 4 ) {
			split = utility::string_split( coords[ coord_start +1 ], ',' );
			origin = Vector( ( static_cast<core::Real>( std::atof( split[0].c_str() ) ) ),
                                             ( static_cast<core::Real>( std::atof( split[1].c_str() ) ) ),
                                             ( static_cast<core::Real>( std::atof( split[2].c_str() ) ) ) );
		}
		push_back_x( x );
		push_back_y( y );
		push_back_origin( origin );
}

} // symmetry
} // conformation
} // core
