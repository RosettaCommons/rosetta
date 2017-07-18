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


#ifndef INCLUDED_protocols_loop_grower_DensSkeleton_hh
#define INCLUDED_protocols_loop_grower_DensSkeleton_hh

#include <protocols/loop_grower/DensSkeleton.fwd.hh>

#include <iostream>
#include <fstream>

#include <protocols/moves/Mover.hh>

#include <ObjexxFCL/format.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
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
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/option.hh>

#include <core/conformation/symmetry/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/electron_density/util.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>

#include <boost/unordered/unordered_map.hpp>



namespace protocols {
namespace loop_grower {

class DensSkeleton {
public:
	DensSkeleton():
		density_loaded_(false)
	{}

	DensSkeleton(std::string filename, core::Real mapreso, core::Real gridspace){
		density_loaded_ = true;
		density_.readMRCandResize(filename,mapreso,gridspace);
	};

	//load the density map
	void
	load_density(std::string filename, core::Real mapreso, core::Real gridspace){
		density_loaded_ = true;
		density_.readMRCandResize(filename,mapreso,gridspace);
	}

	//This function uses a breadth first search to connect two points through the density map. This can blow up very quickly so it must be capped.
	//If max_que is exceeded it will return an empty vector.
	utility::vector1< utility::vector1< numeric::xyzVector < int > > >
	breadth_first_connect(numeric::xyzVector< int > const & start_point, numeric::xyzVector< int >const & end_point, core::Real max_length, bool & hit_max, core::Size max_que = 1e5 );

	//this function takes cartesian coordinates and finds the nearest occupied points in the density before calling the breadth first connect
	utility::vector1< utility::vector1< numeric::xyzVector< int > > >
	breadth_first_connect(numeric::xyzVector< core::Real > const & start_cart, numeric::xyzVector< core::Real > const & end_cart, core::Real max_length, core::Size maxque,
		bool & hit_max, core::Size max_grid );


	//Takes a vector of paths (vector of points) and returns the length of the shorest one
	core::Real
	shortest_path(utility::vector1< utility::vector1< numeric::xyzVector< int > > > paths );

	//takes a cartesian coordinate and finds the closest non-zero grid point
	numeric::xyzVector < int >
	find_closest_occupied_point(numeric::xyzVector< core::Real >const & coord, core::Size max_grid);

	//uses the breadth first search to find the closest way to connect two points(cartesian) through the density (only uses occupied points)
	core::Real
	shortest_path_bfs( numeric::xyzVector< core::Real >const & start, numeric::xyzVector< core::Real >const & end, core::Real max_length, core::Size max_que, bool hit_max, core::Size max_grid);

	//returns whether or not the density has been loaded
	bool
	has_density(){ return density_loaded_; }

	//takes a point and creates a map of all the neighbors to be used in the skeletonization calculations
	std::map< core::Size, numeric::xyzVector< int > >
	assign_neighbors(numeric::xyzVector< int > point);

	//does removing the point break the skeleton
	bool
	breaks_skeleton( utility::vector1< numeric::xyzVector< int > > accepted_points, std::map< core::Size, numeric::xyzVector< int > > neighbors );

	//finds the length of the path plus the distance between the last point and the target end.
	core::Real
	path_length(utility::vector1< numeric::xyzVector< int > > path, numeric::xyzVector< int > end_point );

	//returns a pointer to the density
	inline core::scoring::electron_density::ElectronDensity const & get_density() const { return density_; }

	//Calls the matchAtomFast function using the density_ object
	core::Real
	matchAtomFast( core::Size resid, core::Size atomid, core::conformation::Residue & res, core::pose::Pose& pose );


private:
	bool density_loaded_;
	core::scoring::electron_density::ElectronDensity density_;
};



} //loop_grower
} //protocols

#endif
