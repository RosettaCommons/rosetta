// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/membrane/hull.hh
/// @brief utility functions to compute 2D convex and concave hulls from a slice in a pose
/// @author Julia Koehler Leman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_core_membrane_hull_hh
#define INCLUDED_core_membrane_hull_hh

// Package headers
#include <core/id/AtomID_Map.hh>
#include <core/pose/Pose.hh>

// utility headers
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <numeric/xyzVector.hh>
#include <numeric/util.hh>
#include <utility/vector1.hh>

// C++ headers
#include <cstdlib>

namespace core {
namespace membrane {


////////////////////////////////////////////////////////////////////////////////
// top-level helper functions

/// @brief compute the concave shell from the pose, computes along a slice and
///   takes the atoms within a shell, not just the outermost ones
core::id::AtomID_Map< bool > concave_shell( core::pose::Pose & pose, core::Real min_z, core::Real max_z, core::Real incr_z, core::Real shell_radius, core::Real dist_cutoff );

/// @brief computes 2D convex hull: points in a slice on z-dimension are
///   flattened into xy plane; the output is three lists of points for
///   inside, outside and boundary
utility::vector1< utility::vector1< core::Size > > convex_hull( std::map< core::Size, core::Vector > const & coords, core::Real min_z, core::Real max_z );

/// @brief computes 2D concave hull from three lists of points for inside, outside and boundary
utility::vector1< utility::vector1< core::Size > > concave_hull( std::map< core::Size, core::Vector > const & coords, utility::vector1< utility::vector1< core::Size > > const & point_lists, core::Real dist_cutoff=5.0 );

/// @brief adds a shell to boundary points which are taken from the inside point list
utility::vector1< utility::vector1< core::Size > > add_shell( std::map< core::Size, core::Vector > const & coords, utility::vector1< utility::vector1< core::Size > > const & point_lists, core::Real shell_radius );

////////////////////////////////////////////////////////////////////////////////
// medium-level helper functions

/// @brief compute distances in point list from line p1p2
///   which points: i,o,b for inside, outside or boundary
utility::vector1< core::Real > get_distances( std::map< core::Size, core::Vector > const & coords, core::Size p1, core::Size p2, utility::vector1< core::Size > const & points, bool clock=false );

/// @brief compute enclosing angles in point list from line p1p2
///   which points: i,o,b for inside, outside or boundary
utility::vector1< core::Real > get_angles( std::map< core::Size, core::Vector > const & coords, core::Size p1, core::Size p2, utility::vector1< core::Size > const & points, bool clock=false );

/// @brief find point id that is closest to line connecting p1p2
///   which points: i,o,b for inside, outside or boundary
core::Size find_closest( std::map< core::Size, core::Vector > const & coords, core::Size p1, core::Size p2, utility::vector1< core::Size > const & inside_points, bool clock=false );

/// @brief find point id that is farthest from line connecting p1p2
///   which points: i,o,b for inside, outside or boundary
core::Size find_farthest( std::map< core::Size, core::Vector > const & coords, core::Size p1, core::Size p2, utility::vector1< core::Size > const & outside_points, bool clock=false );

/// @brief checks whether the point q is inside the polygon or not
///   i.e. compute the number of crossing points with polygon
///   0 outside
///   1 inside
///   2 on boundary
std::string inside_polygon( std::map< core::Size, core::Vector > const & coords, utility::vector1< core::Size > const & polygon, core::Vector const & q );

/// @brief returns the points that are within the triangle between points p1, p2, q
utility::vector1< core::Size > points_in_triangle( std::map< core::Size, core::Vector > const & coords, core::Size p1, core::Size p2, core::Size q, utility::vector1< core::Size > const & points );

////////////////////////////////////////////////////////////////////////////////
// bottom-level helper functions

/// @brief distance of point q from line between points p1 and p2
core::Real distance_from_line2D( core::Vector const & p1, core::Vector const & p2, core::Vector const & q );

/// @brief is point q on line segment between points p1 and p2?
bool on_segment( core::Vector const & p1, core::Vector const & p2, core::Vector const & q );

/// @brief q is inside xy boundaries of p1p2
bool inside_boundaries( core::Vector const & p1, core::Vector const & p2, core::Vector const & q );

/// @brief is q clockwise from p1p2?
bool clockwise( core::Vector const & p1, core::Vector const & p2, core::Vector const & q );

/// @brief check whether q is to the left of p1p2
bool to_left( core::Vector const & p1, core::Vector const & p2, core::Vector const & q );

/// @brief check whether line p1p2 intersects with (q to infinity); q is the test point
bool intersect( core::Vector const & p1, core::Vector const & p2, core::Vector const & q );

/// @brief find the sum of angles p1q and p2q; p1, p2, q are points, not vectors
core::Real enclosing_angles( core::Vector const & p1, core::Vector const & p2, core::Vector const & q );

} //core
} //membrane


#endif //core/membrane_hull_hh

