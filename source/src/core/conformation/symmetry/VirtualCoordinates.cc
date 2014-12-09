// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
		utility::vector1< std::string> split ( utility::string_split( coords[ coord_start -1 ], ',' ) );
		assert( split.size() == 3 );
    Vector x( ( static_cast<core::Real>( std::atof( split[1].c_str() ) ) ),
                                             ( static_cast<core::Real>( std::atof( split[2].c_str() ) ) ),
                                             ( static_cast<core::Real>( std::atof( split[3].c_str() ) ) ) );
		split = utility::string_split( coords[ coord_start ], ',' );
		Vector y( ( static_cast<core::Real>( std::atof( split[1].c_str() ) ) ),
                                             ( static_cast<core::Real>( std::atof( split[2].c_str() ) ) ),
                                             ( static_cast<core::Real>( std::atof( split[3].c_str() ) ) ) );
		Vector origin(0,0,0);
		if ( coords.size() == 4 ) {
			split = utility::string_split( coords[ coord_start +2 ], ',' );
			origin = Vector( ( static_cast<core::Real>( std::atof( split[1].c_str() ) ) ),
                                             ( static_cast<core::Real>( std::atof( split[2].c_str() ) ) ),
                                             ( static_cast<core::Real>( std::atof( split[3].c_str() ) ) ) );
		}
		push_back_x( x );
		push_back_y( y );
		push_back_origin( origin );
}

} // symmetry
} // conformation
} // core
