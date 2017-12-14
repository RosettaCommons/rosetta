// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Brandon Frenz
/// @author Frank DiMaio

#include <protocols/loop_grower/DensSkeleton.hh>

#include <iostream>
#include <fstream>
#include <queue>

#include <protocols/moves/Mover.hh>

#include <ObjexxFCL/format.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.hh>

#include <utility/vector1.hh>

#include <basic/Tracer.hh>

//possibily duplicate includes here
#include <basic/basic.hh>
#include <basic/database/open.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/option.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/electron_density/util.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>

#include <boost/unordered/unordered_map.hpp>

static basic::Tracer TR( "protocols.loop_grower.DensSkeleton" );

namespace protocols {
namespace loop_grower {

using namespace core;
using namespace protocols::moves;

core::Real
DensSkeleton::shortest_path_bfs( numeric::xyzVector< core::Real >const & start, numeric::xyzVector< core::Real >const & end, core::Real max_length,
	core::Size max_que, bool hit_max, core::Size max_grid){

	utility::vector1< utility::vector1< numeric::xyzVector< int > > > paths = breadth_first_connect( start, end, max_length, max_que, hit_max, max_grid );
	if ( hit_max == true ) {
		return -1;
	} else if ( paths.size() == 0 ) {
		return -2;
	} else {
		return shortest_path(paths);
	}
	return 0;

}

utility::vector1< utility::vector1< numeric::xyzVector< int > > >
DensSkeleton::breadth_first_connect(numeric::xyzVector< core::Real >const & start_cart, numeric::xyzVector< core::Real >const & end_cart, core::Real max_length, core::Size maxque,
	bool & hit_max, core::Size max_grid ){

	numeric::xyzVector< core::Real > start_point = find_closest_occupied_point(start_cart, max_grid);
	numeric::xyzVector< core::Real > end_point = find_closest_occupied_point(end_cart, max_grid);
	if ( (start_point[0] == 1e5 && start_point[1] == 1e5 && start_point[2] == 1e5) || (end_point[0] == 1e5 && end_point[1] == 1e5 && end_point[2] == 1e5 ) ) {
		hit_max = true;
		TR << " no nearby occupied points could be found. Perhaps try increasing your max_grid distance from " << max_grid << std::endl;
		utility::vector1< utility::vector1< numeric::xyzVector< int > > > empty_vect;
		return empty_vect;
	}

	utility::vector1< utility::vector1< numeric::xyzVector< int > > > paths = breadth_first_connect( start_point, end_point, max_length, hit_max, maxque );
	return paths;

}

utility::vector1< utility::vector1< numeric::xyzVector< int > > >
DensSkeleton::breadth_first_connect(numeric::xyzVector< int >const & start_point, numeric::xyzVector< int >const & end_point, core::Real max_length, bool & hit_max, core::Size max_que ){

	ObjexxFCL::FArray3D< float > const & densdata = density_.get_data();
	std::queue<utility::vector1<numeric::xyzVector< int > > > q;
	utility::vector1< numeric::xyzVector< int > > points;
	utility::vector1< utility::vector1< numeric::xyzVector< int > > > completed_paths;
	points.push_back(start_point);
	q.push(points);
	core::Size level = points.size();
	utility::vector1< utility::vector1< numeric::xyzVector < int > > > sampled_paths;
	sampled_paths.push_back(points);

	///core::Size cntr = 0;
	while ( !q.empty() ) {
		utility::vector1< numeric::xyzVector< int > > path = q.front();
		//if( path.size() == 8 ){
		//    //write_points_to_MRC(path,"path"+utility::to_string(cntr)+".mrc");
		//    write_points_to_pdb(path,"level"+utility::to_string(cntr)+".txt");
		//    cntr++;
		//}
		if ( level != path.size() ) {
			TR << "level: " << path.size() << " que: " << q.size() << std::endl;
			level = path.size();
		}
		q.pop();
		numeric::xyzVector< int > current_point = path.back();
		std::map< core::Size, numeric::xyzVector< int > > neighbors = assign_neighbors(current_point);
		//clean points
		utility::vector1< numeric::xyzVector< int > > clean_path;
		for ( Size i=1; i<=path.size(); i++ ) {
			numeric::xyzVector< int > c_point = path[i];
			if ( !breaks_skeleton(path,neighbors) ) {
				TR << " breaks skeleton " << std::endl;
			} else {
				clean_path.push_back(c_point);
			}

		}
		path = clean_path;
		core::Size total_neighbors = 27;
		for ( Size i=1; i<=total_neighbors; i++ ) {
			//skip values already in the path. Including the current point
			if ( path.has_value(neighbors[i]) ) continue;
			if ( densdata(neighbors[i][0],neighbors[i][1],neighbors[i][2] ) > 1e-10 ) {
				utility::vector1< numeric::xyzVector< int > > next_path = path;
				next_path.push_back(neighbors[i]);

				//if the point is the last point is the end point we are done in which case store the path
				if ( neighbors[i] != end_point ) {
					//abort paths that are too long
					core::Real distance = path_length(next_path,end_point);
					if ( distance < max_length ) {
						//if we've sampled this path already skip it
						utility::vector1< numeric::xyzVector< int > > sorted_path = next_path;
						std::sort(sorted_path.begin(), sorted_path.end()); //, sort_xyz);
						if ( sampled_paths.has_value( sorted_path ) ) continue;
						sampled_paths.push_back(next_path);
						q.push(next_path);
						if ( q.size() > max_que ) {
							TR << " max que exceeded returning empty vector " << std::endl;
							hit_max = true;
							return completed_paths;
						}
					}
				} else {
					completed_paths.push_back(next_path);
				}
			}
		}
	}
	return completed_paths;

}

core::Real
DensSkeleton::shortest_path(utility::vector1< utility::vector1< numeric::xyzVector< int > > > paths ){
	core::Real shortest_path = 0;
	for ( Size i=1; i<=paths.size(); i++ ) {
		core::Real length = 0;
		for ( Size point=1; point<=paths[i].size()-1; point++ ) {
			numeric::xyzVector< core::Real > coord1;
			numeric::xyzVector< core::Real > coord2;
			density_.idx2cart(paths[i][point],coord1);
			density_.idx2cart(paths[i][point+1],coord2);
			length+=coord1.distance(coord2);
		}
		if ( length < shortest_path || shortest_path == 0 ) {
			shortest_path = length;
		}
	}
	return shortest_path;
}

numeric::xyzVector< int >
DensSkeleton::find_closest_occupied_point(numeric::xyzVector< core::Real >const & coord, core::Size max_grid){

	numeric::xyzVector< core::Real > initial_point;
	density_.cart2idx(coord,initial_point);
	core::Real min_distance = -1;
	numeric::xyzVector< int > closest_point(1e5,1e5,1e5);
	ObjexxFCL::FArray3D< float > const & densdata = density_.get_data();
	for ( int x=std::floor(initial_point[0]-max_grid); x<=std::ceil(initial_point[0]+max_grid); x++ ) {
		for ( int y=std::floor(initial_point[1]-max_grid); y<=std::ceil(initial_point[1]+max_grid); y++ ) {
			for ( int z=std::floor(initial_point[2]-max_grid); z<=std::ceil(initial_point[2]+max_grid); z++ ) {
				//skip any points out of bounds
				if ( ( x < 1 || x > densdata.u1() ) || ( y < 1 || y > densdata.u2() ) || ( z < 1 || z > densdata.u3() ) ) continue;
				//skip empty points
				if ( densdata(x,y,z) < 1e-10 ) continue;
				core::Real distance = initial_point.distance( numeric::xyzVector< core::Real >(x,y,z) );
				if ( distance < min_distance or min_distance == -1 ) {
					min_distance = distance;
					closest_point = numeric::xyzVector< int >(x,y,z);
				}
			}
		}
	}
	//no occupied point found return impossible value as signal
	if ( min_distance == -1 ) {
		return numeric::xyzVector< int >(1e5,1e5,1e5);
	}
	return closest_point;
}

std::map< core::Size, numeric::xyzVector< int > >
DensSkeleton::assign_neighbors(numeric::xyzVector< int > point){
	std::map< core::Size, numeric::xyzVector< int > > neighbor_map;
	//z = -1
	neighbor_map[1] = numeric::xyzVector< int >(point[0]-1,point[1]+1,point[2]-1);
	neighbor_map[2] = numeric::xyzVector< int >(point[0],point[1]+1,point[2]-1);
	neighbor_map[3] = numeric::xyzVector< int >(point[0]+1,point[1]+1,point[2]-1);

	neighbor_map[4] = numeric::xyzVector< int >(point[0]-1,point[1],point[2]-1);
	neighbor_map[5] = numeric::xyzVector< int >(point[0],point[1],point[2]-1);
	neighbor_map[6] = numeric::xyzVector< int >(point[0]+1,point[1],point[2]-1);

	neighbor_map[7] = numeric::xyzVector< int >(point[0]-1,point[1]-1,point[2]-1);
	neighbor_map[8] = numeric::xyzVector< int >(point[0],point[1]-1,point[2]-1);
	neighbor_map[9] = numeric::xyzVector< int >(point[0]+1,point[1]-1,point[2]-1);

	//z = 0
	neighbor_map[10] = numeric::xyzVector< int >(point[0]-1,point[1]+1,point[2]);
	neighbor_map[11] = numeric::xyzVector< int >(point[0],point[1]+1,point[2]);
	neighbor_map[12] = numeric::xyzVector< int >(point[0]+1,point[1]+1,point[2]);

	neighbor_map[13] = numeric::xyzVector< int >(point[0]+1,point[1],point[2]);
	neighbor_map[14] = numeric::xyzVector< int >(point[0],point[1],point[2]);
	neighbor_map[15] = numeric::xyzVector< int >(point[0]-1,point[1],point[2]);

	neighbor_map[16] = numeric::xyzVector< int >(point[0]+1,point[1]-1,point[2]);
	neighbor_map[17] = numeric::xyzVector< int >(point[0],point[1]-1,point[2]);
	neighbor_map[18] = numeric::xyzVector< int >(point[0]-1,point[1]-1,point[2]);

	//z = +1
	neighbor_map[19] = numeric::xyzVector< int >(point[0]+1,point[1]+1,point[2]+1);
	neighbor_map[20] = numeric::xyzVector< int >(point[0],point[1]+1,point[2]+1);
	neighbor_map[21] = numeric::xyzVector< int >(point[0]-1,point[1]+1,point[2]+1);

	neighbor_map[22] = numeric::xyzVector< int >(point[0]-1,point[1],point[2]+1);
	neighbor_map[23] = numeric::xyzVector< int >(point[0]-1,point[1],point[2]+1);
	neighbor_map[24] = numeric::xyzVector< int >(point[0]-1,point[1],point[2]+1);

	neighbor_map[25] = numeric::xyzVector< int >(point[0]-1,point[1]-1,point[2]+1);
	neighbor_map[26] = numeric::xyzVector< int >(point[0],point[1]-1,point[2]+1);
	neighbor_map[27] = numeric::xyzVector< int >(point[0]+1,point[1]-1,point[2]+1);

	return neighbor_map;

}

bool
DensSkeleton::breaks_skeleton( utility::vector1< numeric::xyzVector< int > > accepted_points, std::map< core::Size, numeric::xyzVector< int > > neighbors ){

	core::Size total_neighbors = 27;
	utility::vector1< numeric::xyzVector< int > > connected_neighbors;
	for ( Size i=1; i<=total_neighbors; i++ ) {
		if ( !accepted_points.has_value(neighbors[i]) or i == 14 ) {
			continue;
		}
		connected_neighbors.push_back(neighbors[i]);
		//bool done = false;
		core::Size checked_points = 1;
		while ( checked_points != connected_neighbors.size() && checked_points == 0 ) {
			std::map< core::Size, numeric::xyzVector< int > > new_neighbors = assign_neighbors(connected_neighbors[checked_points]);
			for ( Size ii=1; ii<=total_neighbors; ii++ ) {
				if ( accepted_points.has_value(new_neighbors[ii]) ) connected_neighbors.push_back(new_neighbors[i]);
			}
			checked_points+=1;
		}
	}
	bool breaks_skeleton = false;
	for ( core::Size i=1; i<=accepted_points.size(); i++ ) {
		if ( !connected_neighbors.has_value(accepted_points[i]) ) {
			breaks_skeleton = true;
			break;
		}
	}
	return breaks_skeleton;
}

core::Real
DensSkeleton::path_length(utility::vector1< numeric::xyzVector< int > > path, numeric::xyzVector< int > end_point ){

	core::Real distance = 0;
	for ( Size i=1; i<=path.size()-1; i++ ) {
		numeric::xyzVector< core::Real > coord1;
		numeric::xyzVector< core::Real > coord2;
		core::scoring::electron_density::getDensityMap().idx2cart(path[i],coord1);
		core::scoring::electron_density::getDensityMap().idx2cart(path[i+1],coord2);
		distance+=coord1.distance(coord2);
	}

	numeric::xyzVector< core::Real > coord1;
	numeric::xyzVector< core::Real > coord2;
	core::scoring::electron_density::getDensityMap().idx2cart(path.back(),coord1);
	core::scoring::electron_density::getDensityMap().idx2cart(end_point,coord2);

	distance += coord1.distance(coord2);
	return distance;
}

core::Real
DensSkeleton::matchAtomFast( core::Size resid, core::Size atomid, core::conformation::Residue & res, core::pose::Pose& pose ){
	return density_.matchAtomFast(resid, atomid, res, pose, nullptr);
}

} //loop_grower
} //protocols
