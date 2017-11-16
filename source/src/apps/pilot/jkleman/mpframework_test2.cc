// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/apps/pilot/jkleman/mpframework_test.cc
/// @brief  testing random small stuff in MPframework
/// @author  JKLeman (julia.koehler1982@gmail.com)

// Unit Headers
#include <devel/init.hh>

// Package Headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/util.hh>

#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>

#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/util.hh>

// Utility Headers
#include <core/types.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <utility/pointer/owning_ptr.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/conversions.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/util.hh>

// C++ headers
#include <iostream>
#include <cstdlib>

static basic::Tracer TR( "apps.pilot.jkleman.mpframework_test2" );

using namespace core;
using namespace numeric;
using namespace core::id;

////////////////////////////////////////////////////////////////////////////////
// top-level helper functions

/// @brief compute the concave shell from the pose, computes along a slice and
///   takes the atoms within a shell, not just the outermost ones
AtomID_Map< bool > concave_shell( core::pose::Pose & pose, core::Real min_z, core::Real max_z, core::Real incr_z, core::Real shell_radius );

/// @brief computes 2D convex hull: points in a slice on z-dimension are
///   flattened into xy plane; the output is three lists of points for
///   inside, outside and boundary
utility::vector1< utility::vector1< core::Size > > convex_hull( std::map< core::Size, core::Vector > coords, core::Real min_z, core::Real max_z );

/// @brief computes 2D concave hull from three lists of points for inside, outside and boundary
utility::vector1< utility::vector1< core::Size > > concave_hull( std::map< core::Size, core::Vector > coords, utility::vector1< utility::vector1< core::Size > > point_lists, core::Real dist_cutoff=4.0 );

/// @brief adds a shell to boundary points which are taken from the inside point list
utility::vector1< utility::vector1< core::Size > > add_shell( std::map< core::Size, core::Vector > coords, utility::vector1< utility::vector1< core::Size > > point_lists, core::Real shell_radius );

////////////////////////////////////////////////////////////////////////////////
// medium-level helper functions

/// @brief compute distances in point list from line p1p2
///   which points: i,o,b for inside, outside or boundary
utility::vector1< core::Real > get_distances( std::map< core::Size, core::Vector > coords, core::Size p1, core::Size p2, utility::vector1< core::Size > points, bool clock=false );

/// @brief compute enclosing angles in point list from line p1p2
///   which points: i,o,b for inside, outside or boundary
utility::vector1< core::Real > get_angles( std::map< core::Size, core::Vector > coords, core::Size p1, core::Size p2, utility::vector1< core::Size > points, bool clock=false );

/// @brief find point id that is closest to line connecting p1p2
///   which points: i,o,b for inside, outside or boundary
core::Size find_closest( std::map< core::Size, core::Vector > coords, core::Size p1, core::Size p2, utility::vector1< core::Size > inside_points, bool clock=false );

/// @brief find point id that is farthest from line connecting p1p2
///   which points: i,o,b for inside, outside or boundary
core::Size find_farthest( std::map< core::Size, core::Vector > coords, core::Size p1, core::Size p2, utility::vector1< core::Size > outside_points, bool clock=false );

/// @brief checks whether the point q is inside the polygon or not
///   i.e. compute the number of crossing points with polygon
///   0 outside
///   1 inside
///   2 on boundary
std::string inside_polygon( std::map< core::Size, core::Vector > coords, utility::vector1< core::Size > polygon, core::Vector q );

/// @brief returns the points that are within the triangle between points p1, p2, q
utility::vector1< core::Size > points_in_triangle( std::map< core::Size, core::Vector > coords, core::Size p1, core::Size p2, core::Size q, utility::vector1< core::Size > points );

////////////////////////////////////////////////////////////////////////////////
// bottom-level helper functions

/// @brief distance of point q from line between points p1 and p2
core::Real distance_from_line2D( core::Vector p1, core::Vector p2, core::Vector q );

/// @brief is point q on line segment between points p1 and p2?
bool on_segment( core::Vector p1, core::Vector p2, core::Vector q );

/// @brief is q clockwise from p1p2?
bool clockwise( core::Vector p1, core::Vector p2, core::Vector q );

/// @brief check whether q is to the left of p1p2
bool to_left( core::Vector p1, core::Vector p2, core::Vector q );

/// @brief check whether line p1p2 intersects with (q to infinity); q is the test point
bool intersect( core::Vector p1, core::Vector p2, core::Vector q );

/// @brief find the sum of angles p1q and p2q; p1, p2, q are points, not vectors
core::Real enclosing_angles( core::Vector p1, core::Vector p2, core::Vector q );



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// top-level helper functions

////////////////////////////////////////////////////////////////////////////////
/// @brief adds a shell to boundary points which are taken from the inside point list
utility::vector1< utility::vector1< core::Size > > add_shell( std::map< core::Size, core::Vector > coords, utility::vector1< utility::vector1< core::Size > > point_lists, core::Real shell_radius ) {

	// flatten out point lists
	utility::vector1< core::Size > inside = point_lists[ 1 ];
	utility::vector1< core::Size > outside = point_lists[ 2 ];
	utility::vector1< core::Size > boundary = point_lists[ 3 ];

	// initialize empty list
	utility::vector1< core::Size > add_to_boundary;

	// go through points in boundary
	for ( core::Size i = 1; i <= boundary.size(); ++i ) {

		// go through points in inside
		for ( core::Size j = 1; j <= inside.size(); ++j ) {

			// compute distance between point in boundary and in inside
			core::Real dist = coords[ i ].distance( coords[ j ] );

			// if distance is smaller than cutfoff, add to boundary=shell
			if ( dist < shell_radius && ! boundary.has_value( inside[ j ] ) && ! add_to_boundary.has_value( inside[ j ] ) ) {
				add_to_boundary.push_back( inside[ j ] );
			}
		}
	}

	// go through points in add_to_boundary and add them
	for ( core::Size i = 1; i <= add_to_boundary.size(); ++i ) {
		boundary.push_back( add_to_boundary[ i ] );
		inside.erase( std::remove( inside.begin(), inside.end(), add_to_boundary[ i ] ), inside.end() );
	}

	// combine vectors in one big vector to return
	utility::vector1< utility::vector1< core::Size > > new_lists;
	new_lists.push_back( inside );
	new_lists.push_back( outside );
	new_lists.push_back( boundary );

	return new_lists;

}// add shell

////////////////////////////////////////////////////////////////////////////////
/// @brief computes 2D concave hull from three lists of points for inside, outside and boundary
utility::vector1< utility::vector1< core::Size > > concave_hull( std::map< core::Size, core::Vector > coords, utility::vector1< utility::vector1< core::Size > > point_lists, core::Real dist_cutoff ) {

	// flatten out point lists
	utility::vector1< core::Size > inside = point_lists[ 1 ];
	utility::vector1< core::Size > outside = point_lists[ 2 ];
	utility::vector1< core::Size > boundary = point_lists[ 3 ];

	core::Size iter = 0;
	core::Size i = 0;

	// go through iterations
	while ( iter <= 1000 ) {

		iter++;

		// if counter goes over, start again
		if ( i == boundary.size() ) {
			i = 0;
		}

		// get points
		core::Size p1;
		i == 1 ? p1 = boundary.back() : p1 = boundary[ i-1 ];
		core::Size p2 = boundary[ i ];

		// compute distance between point and previous one
		core::Real dist = coords[ p1 ].distance( coords[ p2 ] );

		// if distance larger than cutoff
		if ( dist >= dist_cutoff ) {

			// find closest point in terms of angles
			core::Size min_point = find_closest( coords, p1, p2, inside, true );

			// add this point to the boundary
			if ( boundary.has_value( min_point ) ) {
				i++;
			} else {
				boundary.insert( boundary.begin()+i, min_point );
				i--;
			}

			// compute points in triangle
			utility::vector1< core::Size > move_to_outside = points_in_triangle( coords, p1, p2, min_point, inside );

			// move points in triangle from inside to outside
			for ( core::Size j = 1; i <= move_to_outside.size(); ++j ) {
				outside.push_back( move_to_outside[ j ] );
				inside.erase( std::remove( inside.begin(), inside.end(), move_to_outside[j] ), inside.end() );
			}
		} else {
			// go to next distance
			i++;
		}
	}

	// combine vectors in one big vector to return
	utility::vector1< utility::vector1< core::Size > > new_lists;
	new_lists.push_back( inside );
	new_lists.push_back( outside );
	new_lists.push_back( boundary );

	return new_lists;

}// concave hull

////////////////////////////////////////////////////////////////////////////////
/// @brief computes 2D convex hull: points in a slice on z-dimension are
///   flattened into xy plane; the output is three lists of points for
///   inside, outside and boundary
utility::vector1< utility::vector1< core::Size > > convex_hull( std::map< core::Size, core::Vector > coords, core::Real min_z, core::Real max_z ) {

	// initialize
	utility::vector1< core::Size > inside;
	utility::vector1< core::Size > outside;
	utility::vector1< core::Size > boundary;
	utility::vector1< core::Real > tmp_x;
	utility::vector1< core::Real > tmp_y;

	// go through map, only put points between min and max into the outside list
	for ( core::Size i = 1; i <= coords.size(); ++i ) {
		if ( coords[ i ].z() >= min_z && coords[ i ].z() <= max_z ) {
			outside.push_back( i );
			tmp_x.push_back( coords[ i ].x() );
			tmp_y.push_back( coords[ i ].y() );
		}
	}

	// find points with minimum and maximum x and y
	core::Size minx = tmp_x.index( numeric::min( tmp_x ) );
	core::Size miny = tmp_y.index( numeric::min( tmp_y ) );
	core::Size maxx = tmp_x.index( numeric::max( tmp_x ) );
	core::Size maxy = tmp_y.index( numeric::max( tmp_y ) );

	// get convex hull first using the quickhull method:
	// connect points of smallest and largest x and y, add those points to the boundary
	// order here matters!!!
	boundary.push_back( minx );
	boundary.push_back( miny );
	boundary.push_back( maxx );
	boundary.push_back( maxy );

	// remove those points in the outside point list
	outside.erase( std::remove( outside.begin(), outside.end(), minx ), outside.end() );
	outside.erase( std::remove( outside.begin(), outside.end(), miny ), outside.end() );
	outside.erase( std::remove( outside.begin(), outside.end(), maxx ), outside.end() );
	outside.erase( std::remove( outside.begin(), outside.end(), maxy ), outside.end() );

	// all points within that boundary can be ignored, figure that out through ray casting:
	// odd number of crossings is within polygon, unless point is on boundary
	// go through all points outside
	utility::vector1< core::Size > move_to_inside;

	for ( core::Size i = 1; i <= outside.size(); ++i ) {

		// if point is inside polygon( i.e. boundary ), add to list of points to move
		if ( inside_polygon( coords, boundary, coords[outside[i]] ) == "i" ) {
			move_to_inside.push_back( outside[i] );
		}
	}

	// move points from outside to inside
	for ( core::Size i = 1; i <= move_to_inside.size(); ++i ) {
		inside.push_back( move_to_inside[ i ] );
		outside.erase( std::remove( outside.begin(), outside.end(), move_to_inside[i] ), outside.end() );
	}

	// while there are points in the outside list
	core::Size i = 0;
	while ( outside.size() > 0 ) {

		// get two points in the boundary vector
		core::Size p1 = boundary[ i-1 ];
		core::Size p2 = boundary[ i ];

		// find point in outside list that is farthest away in clockwise direction
		// from the line segment connecting p1 and p2
		core::Size max_point = find_farthest( coords, p1, p2, outside );

		// add this point to the boundary
		if ( ! boundary.has_value( max_point ) ) {
			boundary.insert( boundary.begin()+i, max_point );
			outside.erase( std::remove( outside.begin(), outside.end(), max_point ), outside.end() );
			i--;
		} else {
			i++;
		}

		// compute points in triangle
		utility::vector1< core::Size > move_to_inside = points_in_triangle( coords, p1, p2, max_point, outside );

		// move points in triangle from outside to inside
		for ( core::Size j = 1; j <= move_to_inside.size(); ++j ) {
			inside.push_back( move_to_inside[ j ] );
			outside.erase( std::remove( outside.begin(), outside.end(), move_to_inside[j] ), outside.end() );
		}
	}

	// combine vectors in one big vector to return
	utility::vector1< utility::vector1< core::Size > > point_lists;
	point_lists.push_back( inside );
	point_lists.push_back( outside );
	point_lists.push_back( boundary );

	return point_lists;

}// convex hull

////////////////////////////////////////////////////////////////////////////////
/// @brief compute the concave shell from the pose, computes along a slice and
///   takes the atoms within a shell, not just the outermost ones from the concave hull
AtomID_Map< bool > concave_shell( core::pose::Pose & pose, core::Real min_z, core::Real max_z, core::Real incr_z, core::Real shell_radius ) {

	core::Size i = 0;
	std::map< core::Size, core::Vector > coords;
	AtomID_Map< bool > shell;

	// create coordinates map; go through residues
	for ( core::Size r = 1; r <= pose.total_residue(); ++r ) {

		// skip for membrane residue
		if ( pose.conformation().is_membrane() ) {
			if ( pose.residue( r ).name3() == "MEM" ) {
				continue;
			}
		}

		// go through atoms
		for ( core::Size a = 1; a <= pose.residue( r ).natoms(); ++a ) {
			i++;
			AtomID aid = AtomID( a, r );

			// set coordinates in coordinates map
			coords.insert( std::pair< core::Size, core::Vector >( i, pose.xyz( aid ) ) );

			// set default values in AtomID_map for later
			shell.set( aid, false );
		}
	}

	// go through slices
	for ( core::Real i = min_z; i <= max_z; i += incr_z ) {

		// compute convex hull from coordinates
		utility::vector1< utility::vector1< core::Size > > convex_lists = convex_hull( coords, i, i+incr_z );

		// compute concave hull from coordinates
		utility::vector1< utility::vector1< core::Size > > concave_lists = concave_hull( coords, convex_lists );

		// add shell around concave hull
		utility::vector1< utility::vector1< core::Size > > shell_lists = add_shell( coords, concave_lists, shell_radius );

		// get information from point lists
		// the points are identifiers that correspond to the i's in the coords map above
		utility::vector1< core::Size > inside = shell_lists[1];
		utility::vector1< core::Size > outside = shell_lists[2];
		utility::vector1< core::Size > boundary = shell_lists[3];

		// add information in the boundary list to AtomID_Map: if any of the atoms
		// in a residue is part of the boundary, set the entire residue as part of the shell
		core::Size id = 0;
		for ( core::Size r = 1; r <= shell.n_residue(); ++r ) {

			// initialize flag
			bool in_boundary( false );

			// go through atoms to see whether any of the atoms in that residue
			// are part of the boundary
			for ( core::Size a = 1; a <= shell.n_atom( r ); ++a ) {
				id += 1;
				AtomID aid = AtomID( a, r );

				// residue is in boundary
				if ( shell.get( aid ) == boundary[id] ) {
					in_boundary = true;
				}
			}

			// if any of the atoms of that residue is part of the boundary
			if ( in_boundary == true ) {

				// go through atoms again and set all atoms for that residue to true
				// in AtomIDMap
				for ( core::Size a = 1; a <= shell.n_atom( r ); ++a ) {
					AtomID aid = AtomID( a, r );
					shell.set( aid, true );
				}// atoms
			}// in boundary
		}// residue
	}// slices

	return shell;

}// concave shell



////////////////////////////////////////////////////////////////////////////////
// medium-level helper functions
////////////////////////////////////////////////////////////////////////////////
/// @brief compute distances in point list from line p1p2
utility::vector1< core::Real > get_distances( std::map< core::Size, core::Vector > coords, core::Size p1, core::Size p2, utility::vector1< core::Size > points, bool clock ) {

	// initialize empty vector
	utility::vector1< core::Real > distances;

	// iterate over point list
	for ( core::Size i = 1; i <= points.size(); ++i ) {

		core::Size q = points[ i ];

		// keep going if we are looking at p1 or p2
		if ( q == p1 || q == p2 ) {
			continue;
		}

		// find the point with the largest distance from it in clockwise direction (i.e. outside)
		// if point is in counter-clockwise direction, then ignore for now
		if ( clockwise( coords[ p1 ], coords[ p2 ], coords[ q ] ) == clock ) {
			continue;
		}

		// add distance value into output vector
		distances.push_back( distance_from_line2D( coords[ p1 ], coords[ p2 ], coords[ q ] ) );

	}
	return distances;

}// get distances

////////////////////////////////////////////////////////////////////////////////
/// @brief compute enclosing angles in point list from line p1p2
utility::vector1< core::Real > get_angles( std::map< core::Size, core::Vector > coords, core::Size p1, core::Size p2, utility::vector1< core::Size > points, bool clock ) {

	// initialize empty vector
	utility::vector1< core::Real > angles;

	// iterate over point list
	for ( core::Size i = 1; i <= points.size(); ++i ) {

		core::Size q = points[ i ];

		// keep going if we are looking at p1 or p2
		if ( q == p1 || q == p2 ) {
			continue;
		}

		// find the point with the smalles enclosing angles in clockwise direction (i.e. inside)
		// if point is in counter-clockwise direction, then ignore for now
		if ( clockwise( coords[ p1 ], coords[ p2 ], coords[ q ] ) == clock ) {
			continue;
		}

		// add angle value into output vector
		angles.push_back( enclosing_angles( coords[ p1 ], coords[ p2 ], coords[ q ] ) );

	}
	return angles;

}// get enclosing angles

////////////////////////////////////////////////////////////////////////////////
/// @brief find point in AtomIDMap that is closest to line connecting p1p2
core::Size find_closest( std::map< core::Size, core::Vector > coords, core::Size p1, core::Size p2, utility::vector1< core::Size > inside_points, bool clock ) {

	utility::vector1< core::Size > angles = get_angles( coords, p1, p2, inside_points, clock );

	// if list is empty, quit
	if ( angles.size() == 0 ) {
		return p1;
	}

	// find point with smallest combined angle
	core::Real min_angle = numeric::min( angles );
	core::Size min_point = inside_points[ angles.index( min_angle ) ];

	// overwrite min point if the angle is larger than 90 degrees
	if ( min_angle > 130.0 ) {
		min_point = p1;
	}

	return min_point;

} // find point id for closest point

////////////////////////////////////////////////////////////////////////////////
/// @brief find point in AtomIDMap that is closest to line connecting p1p2
core::Size find_farthest( std::map< core::Size, core::Vector > coords, core::Size p1, core::Size p2, utility::vector1< core::Size > outside_points, bool clock ) {

	utility::vector1< core::Size > distances = get_distances( coords, p1, p2, outside_points, clock );

	// if list is empty, quit
	if ( distances.size() == 0 ) {
		return p1;
	}

	// find point with maximal distance
	core::Real max_dist = numeric::max( distances );
	core::Size max_point = outside_points[ distances.index( max_dist ) ];

	return max_point;

} // find point id for farthest point

////////////////////////////////////////////////////////////////////////////////
/// @brief checks whether the point q is inside the polygon or not
///   i.e. compute the number of crossing points with polygon
///   0 outside
///   1 inside
///   2 on boundary
std::string inside_polygon( std::map< core::Size, core::Vector > coords, utility::vector1< core::Size > polygon, core::Vector q ) {

	// initialize counts
	core::SSize cnt = 0;

	// go through point pairs on polygon
	for ( core::Size i = 1; i <= polygon.size(); ++i ) {

		// get points
		core::Size p1;
		i == 1 ? p1 = polygon.back() : p1 = polygon[ i-1 ];
		core::Size p2 = polygon[ i ];

		// if (q to infinity) intersects with p1p2, increase counter
		if ( intersect( coords[ p1 ], coords[ p2 ], q ) == true ) {
			cnt++;
		}

		// if q lies on segment between two points
		if ( on_segment( coords[ p1 ], coords[ p2 ], q ) ) {
			cnt--;
			break;
		}
	}

	// if odd number of counts, lies inside
	if ( cnt % 2 == 1 ) {
		return "i";
	} else if ( cnt == -1 ) {
		return "b";
	}
	return "o";

}// inside polygon

////////////////////////////////////////////////////////////////////////////////
/// @brief returns the points that are within the triangle between points p1, p2, q
utility::vector1< core::Size > points_in_triangle( std::map< core::Size, core::Vector > coords, core::Size p1, core::Size p2, core::Size q, utility::vector1< core::Size > points ) {

	utility::vector1< core::Size > in_triangle;

	// create triangle polygon
	utility::vector1< core::Size > triangle;
	triangle.push_back( p1 );
	triangle.push_back( p2 );
	triangle.push_back( q );

	// go through point list
	for ( core::Size i = 1; i <= points.size(); ++i ) {

		// if p1, p2, or q in point list, keep going
		if ( points[ i ] == p1 ) { continue; }
		if ( points[ i ] == p2 ) { continue; }
		if ( points[ i ] == q ) { continue; }

		// if counter-clockwise, ignore
		if ( clockwise( coords[ p1 ], coords[ p2 ], coords[ points[ i ] ] ) == false ) {
			continue;
		}

		// if point is within triangle, append to list to return
		if ( inside_polygon( coords, triangle, coords[ points[ i ] ] ) == "i" ) {
			in_triangle.push_back( points[ i ] );
		}
	}

	return in_triangle;

}// points in triangle


////////////////////////////////////////////////////////////////////////////////
// bottom-level helper functions

/// @brief distance of point q from line between points p1 and p2
core::Real distance_from_line2D( core::Vector p1, core::Vector p2, core::Vector q ) {

	core::Real distance = std::abs( (p2.y()-p1.y())*q.x() - (p2.x()-p1.x())*q.y() + p2.x()*p1.y() - p2.y()*p1.x() ) / sqrt( pow((p2.y()-p1.y()), 2) + pow((p2.x()-p1.x()), 2) );
	return distance;

}// distance from line2D

////////////////////////////////////////////////////////////////////////////////
/// @brief is point q on line segment between points p1 and p2?
bool on_segment( core::Vector p1, core::Vector p2, core::Vector q ) {

	core::Real slope;

	// division by zero catch
	if ( p1.x() == p2.x() ) {
		slope = 0;
	} else {
		// compute slope of line between p1 and p2
		slope = ( p1.y() - p2.y() ) / ( p1.x() - p2.x() );
	}

	// point lies on line
	if ( q.y() == slope * ( q.x() - p1.x() ) + p1.y() ) {

		// point lies on line within boundaries, i.e. on segment
		if ( q.x() >= min( p1.x(), p2.x() ) && q.x() <= max( p1.x(), p2.x() ) &&
				q.y() >= min( p1.y(), p2.y() ) && q.y() <= max( p1.y(), p2.y() ) ) {
			return true;
		}
	}

	// point does not lie on segment
	return false;

}// on segment

////////////////////////////////////////////////////////////////////////////////
/// @brief is q clockwise from p1p2?
bool clockwise( core::Vector p1, core::Vector p2, core::Vector q ) {

	// get cross-product in 2D, if 3rd dimension is zero
	core::Real val = ( p2.x()-p1.x() )*( q.y()-p2.y() ) - ( p2.y()-p1.y() )*( q.x()-p2.x() );

	if ( val > 0 ) {
		return false;
	} else {
		return true;
	}
}// clockwise

////////////////////////////////////////////////////////////////////////////////
/// @brief check whether q is to the left of p1p2
bool to_left( core::Vector p1, core::Vector p2, core::Vector q ) {

	// check whether q is clockwise from p1p2
	if ( p1.y() <= p2.y() && clockwise( p1, p2, q ) == false ) {
		return true;
	} else if ( p1.y() >= p2.y() && clockwise( p1, p2, q) == true ) {
		return true;
	}

	return false;

}// to left

////////////////////////////////////////////////////////////////////////////////
/// @brief check whether line p1p2 intersects with (q to infinity); q is the test point
bool intersect( core::Vector p1, core::Vector p2, core::Vector q ) {

	// if slope is parallel, lines don't intersect
	if ( p1.x() == p2.x() ) {
		return false;
	}

	// compute the slope of p1p2
	core::Real slope = ( p1.y() - p2.y() )/( p1.x() - p2.x() );

	// check whether q test point is within y-boundaries of p1p2
	if ( q.y() >= min( p1.y(), p2.y() ) && q.y() <= max( p1.y() ,p2.y() ) ) {

		// check whether q is to the left of p1p2
		if ( to_left( p1, p2, q ) ) {
			return true;
		}
	}

	// else: they don't intersect
	return false;

}// intersect

////////////////////////////////////////////////////////////////////////////////
/// @brief find the sum of angles p1q and p2q; p1, p2, q are points, not vectors
core::Real enclosing_angles( core::Vector p1, core::Vector p2, core::Vector q ) {

	using namespace numeric;

	core::Vector a = q - p1;
	core::Vector b = q - p2;
	core::Vector c = p2 - p1;
	core::Vector mc = p1 - p2;

	core::Real angle_ac = numeric::conversions::degrees( angle_of( a, c ) );
	core::Real angle_bc = numeric::conversions::degrees( angle_of( b, mc ) );

	return angle_ac+angle_bc;

}// enclosing angles


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//typedef utility::pointer::shared_ptr< MPframeworkTestMover > MPframeworkTestMoverOP;

/// @brief Main method
int
main( int argc, char * argv [] )
{
	try {

		devel::init(argc, argv);

		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "4DXW_opm_ABCD.pdb" , core::import_pose::PDB_file);

		core::Real min_z = -15.0;
		core::Real max_z = 15.0;
		core::Real incr_z = 5.0;
		core::Real shell_radius = 5.0;

		AtomID_Map< bool > shell = concave_shell( pose, min_z, max_z, incr_z, shell_radius );


		return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

