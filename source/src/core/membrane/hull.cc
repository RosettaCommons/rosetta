// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/membrane/hull.cc
/// @brief utility functions to compute 2D convex and concave hulls from a slice in a pose
/// @author Julia Koehler Leman (julia.koehler1982@gmail.com)

#include <core/membrane/hull.hh>

// Package headers
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>

// utility headers
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <numeric/conversions.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/util.hh>
#include <utility/vector1.hh>

// C++ headers
#include <cstdlib>
#include <cmath>

static THREAD_LOCAL basic::Tracer TR( "core.membrane.hull" );


namespace core {
namespace membrane {

	
////////////////////////////////////////////////////////////////////////////////
// top-level helper functions

////////////////////////////////////////////////////////////////////////////////
/// @brief compute the concave shell from the pose, computes along a slice and
///			takes the atoms within a shell, not just the outermost ones from the concave hull
core::id::AtomID_Map< bool > concave_shell( core::pose::Pose & pose, core::Real min_z, core::Real max_z, core::Real incr_z, core::Real shell_radius, core::Real dist_cutoff ) {
	
	using namespace core::id;
	
	TR << "concave shell" << std::endl;
	
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
		
		// set coordinates in coordinates map and set default values in AtomID_map for later
		AtomID aid = AtomID( 2, r );
		
		// only insert if atom in membrane
		coords.insert( std::pair< core::Size, core::Vector >( r, pose.conformation().xyz( aid ) ) );
		shell.set( aid, false );
	}
	
	// go through slices
	for ( core::Real i = min_z; i < max_z; i += incr_z ) {
		
		// compute convex hull from coordinates
		utility::vector1< utility::vector1< core::Size > > convex_lists = convex_hull( coords, i, i+incr_z );
		
		// compute concave hull from coordinates
		utility::vector1< utility::vector1< core::Size > > concave_lists = concave_hull( coords, convex_lists, dist_cutoff );

		// add shell around concave hull
		utility::vector1< utility::vector1< core::Size > > shell_lists = add_shell( coords, concave_lists, shell_radius );

		// get information from point lists
		// the points are identifiers that correspond to the i's in the coords map above
//		utility::vector1< core::Size > inside = convex_lists[1];
//		utility::vector1< core::Size > outside = convex_lists[2];
//		utility::vector1< core::Size > boundary = convex_lists[3];

//		utility::vector1< core::Size > inside = concave_lists[1];
//		utility::vector1< core::Size > outside = concave_lists[2];
//		utility::vector1< core::Size > boundary = concave_lists[3];

		utility::vector1< core::Size > inside = shell_lists[1];
		utility::vector1< core::Size > outside = shell_lists[2];
		utility::vector1< core::Size > boundary = shell_lists[3];

		// check lists
//		TR << "after adding shell:" << std::endl;
//		TR << "inside " << std::endl;
//		for ( core::Size k = 1; k <= inside.size(); ++k ) {
//			TR << inside[k] << " ";
//		}
//		TR << std::endl;
//		TR << "outside " << std::endl;
//		for ( core::Size k = 1; k <= outside.size(); ++k ) {
//			TR << outside[k] << " ";
//		}
//		TR << std::endl;
//		TR << "boundary " << std::endl;
//		for ( core::Size k = 1; k <= boundary.size(); ++k ) {
//			TR << boundary[k] << " ";
//		}
//		TR << std::endl;

		// add information in the boundary list to AtomID_Map: if any of the atoms
		// in a residue is part of the boundary, set the entire residue as part of the shell
		for ( core::Size r = 1; r <= shell.n_residue(); ++r ) {
			
			// residue is in boundary
			if ( boundary.has_value( r ) ) {
				AtomID aid = AtomID( 2, r );
				if ( shell.has( aid ) ) {
					shell.set( aid, true );
				}
			}
		}// residue
	}// slices
	
	return shell;
	
}// concave shell

////////////////////////////////////////////////////////////////////////////////
/// @brief adds a shell to boundary points which are taken from the inside point list
utility::vector1< utility::vector1< core::Size > > add_shell( std::map< core::Size, core::Vector > coords, utility::vector1< utility::vector1< core::Size > > point_lists, core::Real shell_radius ) {
	
	TR << "adding shell" << std::endl;
	
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
			core::Vector a = coords[ boundary[ i ] ];
			core::Vector b = coords[ inside[ j ] ];
			a.z() = 0;
			b.z() = 0;
			core::Real dist = a.distance( b );
			
			// if distance is smaller than cutoff, add to boundary=shell
			if ( dist < shell_radius && boundary.has_value( inside[ j ] ) == false && add_to_boundary.has_value( inside[ j ] ) == false ) {
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
	
	TR << "concave hull" << std::endl;
	
	// flatten out point lists
	utility::vector1< core::Size > inside = point_lists[ 1 ];
	utility::vector1< core::Size > outside = point_lists[ 2 ];
	utility::vector1< core::Size > boundary = point_lists[ 3 ];
	
	core::SSize i = 1;
	core::Size iter = 0;
	
	// go through iterations
	while ( iter <= 1000 ) {
		
		iter++;
		
		// wrap around vector
		if ( static_cast< core::Size >( i ) > boundary.size() ) {
			i -= boundary.size();
		}
		if ( i <= 0 ) {
			i += boundary.size();
		}
		
		// get points
		core::Size p1;
		i == 1 ? p1 = boundary.back() : p1 = boundary[ i-1 ];
		core::Size p2 = boundary[ i ];
		
		core::Vector c1 = coords[ p1 ];
		core::Vector c2 = coords[ p2 ];
		
		// compute distance between point and previous one
		core::Real dist = sqrt( pow((c1.x()-c2.x()), 2) + pow((c1.y()-c2.y()), 2) );
		
		// if distance larger than cutoff
		if ( dist >= dist_cutoff ) {
			
			// find closest point in terms of angles
			core::Size min_point = find_closest( coords, p1, p2, inside, true );
			
			// add this point to the boundary
			if ( boundary.has_value( min_point ) ) {
				i++;
			}
			else {
				boundary.insert( boundary.begin()+i-1, min_point );
				i--;
			}
			
			// compute points in triangle
			utility::vector1< core::Size > move_to_outside = points_in_triangle( coords, p1, p2, min_point, inside );
			
			// move points in triangle from inside to outside
			for ( core::Size j = 1; j <= move_to_outside.size(); ++j ) {
				outside.push_back( move_to_outside[ j ] );
				inside.erase( std::remove( inside.begin(), inside.end(), move_to_outside[j] ), inside.end() );
			}
		}
		
		// go to next distance
		else {
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
///			flattened into xy plane; the output is three lists of points for
///			inside, outside and boundary
utility::vector1< utility::vector1< core::Size > > convex_hull( std::map< core::Size, core::Vector > coords, core::Real min_z, core::Real max_z ) {
	
	TR << "convex hull" << std::endl;
	
	// initialize
	utility::vector1< core::Size > inside;
	utility::vector1< core::Size > outside;
	utility::vector1< core::Size > boundary;
	utility::vector1< core::Real > tmp_x;
	utility::vector1< core::Real > tmp_y;
	
	// go through map, only put points between min and max into the outside list
	for ( core::Size i = 1; i <= coords.size(); ++i ) {
		if ( coords[ i ].z() >= min_z && coords[ i ].z() <= max_z ) {
			outside.push_back( coords.find( i )->first );
			tmp_x.push_back( coords.find( i )->second.x() );
			tmp_y.push_back( coords.find( i )->second.y() );
		}
	}
	
	// find points with minimum and maximum x and y
	core::Size minx = outside[ tmp_x.index( numeric::min( tmp_x ) ) ];
	core::Size miny = outside[ tmp_y.index( numeric::min( tmp_y ) ) ];
	core::Size maxx = outside[ tmp_x.index( numeric::max( tmp_x ) ) ];
	core::Size maxy = outside[ tmp_y.index( numeric::max( tmp_y ) ) ];

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
	core::SSize i = 1;
	core::Size length = outside.size();
	core::Size iter = 0;
	
	while( length > 0 && iter <= 10 ) {
		
		// wrap around vector
		if ( static_cast< core::Size >( i ) > boundary.size() ) {
			i -= boundary.size();
		}
		if ( i <= 0 ) {
			i += boundary.size();
		}
		
		// get two points in the boundary vector
		core::Size p1;
		i == 1 ? p1 = boundary.back() : p1 = boundary[ static_cast< core::Size >( i-1 ) ];
		core::Size p2 = boundary[ static_cast< core::Size >( i ) ];

		// find point in outside list that is farthest away in clockwise direction
		// from the line segment connecting p1 and p2
		core::Size max_point = find_farthest( coords, p1, p2, outside );

		// add this point to the boundary
		if ( ! boundary.has_value( max_point ) ) {
			boundary.insert( boundary.begin()+i-1, max_point );
			outside.erase( std::remove( outside.begin(), outside.end(), max_point ), outside.end() );
			i--;
		}
		else {
			i++;
		}
		
		// compute points in triangle
		utility::vector1< core::Size > move_to_inside = points_in_triangle( coords, p1, p2, max_point, outside );
		
		// move points in triangle from outside to inside
		for ( core::Size j = 1; j <= move_to_inside.size(); ++j ) {
			inside.push_back( move_to_inside[ j ] );
			outside.erase( std::remove( outside.begin(), outside.end(), move_to_inside[ j ] ), outside.end() );
		}
		
		length = outside.size();
		
		// workaround for that one time, where one of the points is exactly on the line
		if ( length == 1 ) {
			iter++;
		}
	}

	outside.push_back( 0 );
	
	// combine vectors in one big vector to return
	utility::vector1< utility::vector1< core::Size > > point_lists;
	point_lists.push_back( inside );
	point_lists.push_back( outside );
	point_lists.push_back( boundary );
	
	return point_lists;
	
}// convex hull

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
			distances.push_back( 0.0 );
		}
		
		// find the point with the largest distance from it in clockwise direction (i.e. outside)
		// if point is in counter-clockwise direction, then ignore for now
		else if ( clockwise( coords[ p1 ], coords[ p2 ], coords[ q ] ) == clock ) {
			distances.push_back( -1000 );
		}
		// add distance value into output vector
		else {
			distances.push_back( distance_from_line2D( coords[ p1 ], coords[ p2 ], coords[ q ] ) );
		}
	}
	return distances;
	
}// get distances

////////////////////////////////////////////////////////////////////////////////
/// @brief compute enclosing angles in point list from line p1p2
utility::vector1< core::Real > get_angles( std::map< core::Size, core::Vector > coords, core::Size p1, core::Size p2, utility::vector1< core::Size > points, bool clock ) {
	
//	TR << "get angles" << std::endl;
	
	// initialize empty vector
	utility::vector1< core::Real > angles;
	
	// iterate over point list
	for ( core::Size i = 1; i <= points.size(); ++i ) {
		
		core::Size q = points[ i ];
		
		// keep going if we are looking at p1 or p2
		if ( q == p1 || q == p2 ) {
			angles.push_back( 9999 );
		}
		
		// find the point with the smalles enclosing angles in clockwise direction (i.e. inside)
		// if point is in counter-clockwise direction, then ignore for now
		else if ( clockwise( coords[ p1 ], coords[ p2 ], coords[ q ] ) == clock ) {
			angles.push_back( 9999 );
		}
		// add angle value into output vector
		else {
			angles.push_back( enclosing_angles( coords[ p1 ], coords[ p2 ], coords[ q ] ) );
		}
	}
	return angles;
	
}// get enclosing angles

////////////////////////////////////////////////////////////////////////////////
/// @brief find point in AtomIDMap that is closest to line connecting p1p2
core::Size find_closest( std::map< core::Size, core::Vector > coords, core::Size p1, core::Size p2, utility::vector1< core::Size > inside_points, bool clock ) {
	
	utility::vector1< core::Real > angles = get_angles( coords, p1, p2, inside_points, clock );
	
	// if list is empty, quit
	if ( angles.size() == 0 || numeric::min( angles ) == 9999 ) {
		return p1;
	}
	
	// find point with smallest combined angle
	core::Real min_angle = numeric::min( angles );
	core::Size min_point = inside_points[ angles.index( min_angle ) ];
	
	// overwrite min point if the angle is larger than X degrees
	if ( min_angle > 130.0 ) {
		min_point = p1;
	}
	
	return min_point;
	
} // find point id for closest point

////////////////////////////////////////////////////////////////////////////////
/// @brief find point in AtomIDMap that is farthest from line connecting p1p2
core::Size find_farthest( std::map< core::Size, core::Vector > coords, core::Size p1, core::Size p2, utility::vector1< core::Size > outside_points, bool clock ) {
	
	// check length of input vector
	if ( outside_points.size() == 0 ) {
		return p1;
	}
	
	// write the distances in a vector
	utility::vector1< core::Real > distances = get_distances( coords, p1, p2, outside_points, clock );
	
	// if list is empty, quit
	if ( distances.size() == 0 || numeric::max( distances ) == -1000 ) {
		return p1;
	}

	// find point with maximal distance
	core::Real max_dist = numeric::max( distances );
	core::Size max_point = outside_points[ distances.index( max_dist ) ];
	
	return max_point;
	
} // find point id for farthest point

////////////////////////////////////////////////////////////////////////////////
/// @brief checks whether the point q is inside the polygon or not
///			i.e. compute the number of crossing points with polygon
///			0 outside
///			1 inside
///			2 on boundary
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
		if ( core::membrane::intersect( coords[ p1 ], coords[ p2 ], q ) == true ) {
			cnt += 1;
		}
		
		// if q lies on segment between two points
		if ( core::membrane::on_segment( coords[ p1 ], coords[ p2 ], q )
				|| p1 == p2 || coords[ p1 ] == q || coords[ p2 ] == q ) {
			cnt = -1;
			break;
		}
	}

	// if odd number of counts, lies inside
	if ( cnt % 2 == 1 ) {
		return "i";
	}
	else if ( cnt == -1 ) {
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
//		if ( clockwise( coords[ p1 ], coords[ p2 ], coords[ points[ i ] ] ) == false ) {
//			continue;
//		}
		
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
	
	// checks
	if ( p1 == p2 || p1 == q || p2 == q ) {
		return 0.0;
	}
	
	core::Real numer = std::abs( (p2.y()-p1.y())*q.x() - (p2.x()-p1.x())*q.y() + p2.x()*p1.y() - p2.y()*p1.x() );
	core::Real denom = sqrt( pow((p2.y()-p1.y()), 2) + pow((p2.x()-p1.x()), 2) );
	core::Real dist = numer / denom;
	
	return dist;
	
}// distance from line2D

////////////////////////////////////////////////////////////////////////////////
/// @brief is point q on line segment between points p1 and p2?
bool on_segment( core::Vector p1, core::Vector p2, core::Vector q ) {
	
	bool value( false );
	
	using namespace numeric;
	core::Real slope;
	
	// division by zero catch
	if ( p1.x() == p2.x() ){
		
		if ( inside_boundaries( p1, p2, q ) == true ) {
			value = true;
		}
	}
	// compute slope of line between p1 and p2
	else {
		slope = ( p1.y() - p2.y() ) / ( p1.x() - p2.x() );

		if ( ( q.y() == p1.y() + slope * ( q.x() - p1.x() ) )
			&& inside_boundaries( p1, p2, q ) == true ) {
			value = true;
		}
	}
	
	return value;
	
}// on segment

////////////////////////////////////////////////////////////////////////////////
/// @brief q is inside xy boundaries of p1p2
bool inside_boundaries( core::Vector p1, core::Vector p2, core::Vector q ) {
	
	using namespace numeric;

	if ( ( q.x() >= min( p1.x(), p2.x() ) ) && ( q.x() <= max( p1.x(), p2.x() ) )
		&& ( q.y() >= min( p1.y(), p2.y() ) ) && ( q.y() <= max( p1.y(), p2.y() ) ) ) {
		return true;
	}
	return false;
	
}// inside xy boundaries

////////////////////////////////////////////////////////////////////////////////
/// @brief is q clockwise from p1p2?
bool clockwise( core::Vector p1, core::Vector p2, core::Vector q ) {
	
	// get cross-product in 2D, if 3rd dimension is zero
	core::Real val = ( p2.x()-p1.x() )*( q.y()-p2.y() ) - ( p2.y()-p1.y() )*( q.x()-p2.x() );
	
	if ( val > 0 ) {
		return false;
	}
	else {
		return true;
	}
}// clockwise

////////////////////////////////////////////////////////////////////////////////
/// @brief check whether q is to the left of p1p2
bool to_left( core::Vector p1, core::Vector p2, core::Vector q ) {
	
	// check whether q is clockwise from p1p2
	if ( p1.y() <= p2.y() && clockwise( p1, p2, q ) == false ) {
		return true;
	}
	else if ( p1.y() >= p2.y() && clockwise( p1, p2, q) == true ) {
		return true;
	}
	
	return false;
	
}// to left

////////////////////////////////////////////////////////////////////////////////
/// @brief check whether line p1p2 intersects with (q to infinity); q is the test point
bool intersect( core::Vector p1, core::Vector p2, core::Vector q ) {

	using namespace numeric;

	// if slope is parallel, lines don't intersect
	if ( p1.y() == p2.y() ) {
		return false;
	}
	
	// check whether q test point is within y-boundaries of p1p2
	if ( q.y() >= min( p1.y(), p2.y() ) && q.y() <= max( p1.y() ,p2.y() ) ) {
		
		// check whether q is to the left of p1p2
		if ( to_left( p1, p2, q ) == true ) {
			return true;
		}
	}
	
	// else: they don't intersect
	return false;
	
}// intersect

////////////////////////////////////////////////////////////////////////////////
/// @brief find the sum of angles p1q and p2q in 2D; p1, p2, q are points, not vectors
core::Real enclosing_angles( core::Vector p1, core::Vector p2, core::Vector q ) {
	
	using namespace numeric;
	
	core::Vector a = q - p1;
	core::Vector b = q - p2;
	core::Vector c = p2 - p1;
	core::Vector mc = p1 - p2;
	
	// this is in 2D
	a.z() = 0;
	b.z() = 0;
	c.z() = 0;
	mc.z() = 0;
	
	core::Real angle_ac = numeric::conversions::degrees( angle_of( a, c ) );
	core::Real angle_bc = numeric::conversions::degrees( angle_of( b, mc ) );
	
	return angle_ac+angle_bc;
	
}// enclosing angles


} //core
} //membrane


