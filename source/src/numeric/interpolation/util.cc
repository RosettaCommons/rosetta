// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/numeric/interpolation/util.cc
/// @author Sam DeLuca
/// @author Stephanie Hirst

#include <numeric/interpolation/util.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <iostream>

#include <boost/algorithm/string.hpp>

namespace numeric {
namespace interpolation {

/// @details read in a file, read out a spline.
///The file should be tab separated, and have two lines
///The first field of one line should be "x_axis", the next fields should be the x values of the points for the spline
///The first field of the other line should be "y_axis", the next fields should be the y values for the points on the spline
///For an example, see "scoring/constraints/epr_distance_potential.histogram"
spline::SplineGenerator spline_from_file(std::string const &  filename,platform::Real const & bin_size)
{
	utility::io::izstream potential_file;
	potential_file.open(filename.c_str());

	utility::vector1<platform::Real> bins_vect;
	utility::vector1<platform::Real> potential_vect;

	platform::Real lower_bound_x(0.0);
	platform::Real lower_bound_y(0.0);
	platform::Real upper_bound_x(0.0);
	platform::Real upper_bound_y(0.0);


	std::string line;
	utility::vector1< std::pair<std::string,spline::LinearFunction> > boundary_functions;

	while ( getline(potential_file,line) )
			{
		if ( line.size() == 0 ) {
			continue;
		}
		boost::trim(line);
		utility::vector1<std::string>  split_fields(utility::string_split(line,'\t'));
		if ( split_fields[1] == "x_axis" ) {
			for ( platform::Size count = 2; count <= split_fields.size(); ++count ) {
				bins_vect.push_back(utility::from_string(split_fields[count],platform::Real(0.0)));
			}

			utility::vector1<platform::Real>::const_iterator lower_bound_x_iterator(bins_vect.begin());
			utility::vector1<platform::Real>::const_reverse_iterator upper_bound_x_iterator(bins_vect.rbegin());
			lower_bound_x = ((*lower_bound_x_iterator) - (0.5*bin_size));
			upper_bound_x = ((*upper_bound_x_iterator) + (0.5*bin_size));
		} else if ( split_fields[1] == "y_axis" ) {
			for ( platform::Size count = 2; count <= split_fields.size(); ++count ) {
				potential_vect.push_back(utility::from_string(split_fields[count],platform::Real(0.0)));
			}

			utility::vector1<platform::Real>::const_iterator lower_bound_y_iterator(potential_vect.begin());
			utility::vector1<platform::Real>::const_reverse_iterator upper_bound_y_iterator(potential_vect.end());
			lower_bound_y = *lower_bound_y_iterator;
			upper_bound_y = *upper_bound_y_iterator;
		}
		if ( split_fields[1] == "lb_function" || split_fields[1] == "ub_function" ) {
			platform::Real cutoff= utility::from_string(split_fields[2],platform::Real(0.0));
			platform::Real slope = utility::from_string(split_fields[3],platform::Real(0.0));
			platform::Real intercept = utility::from_string(split_fields[4],platform::Real(0.0));
			boundary_functions.push_back(std::make_pair(split_fields[1],spline::LinearFunction(cutoff,slope,intercept)));
		}

	}
	potential_file.close();
	if ( bins_vect.size() != potential_vect.size() ) {
		utility_exit_with_message("the x and y axis lines in the spline file"+filename+"do not have the same number of entries");
	}

	numeric::interpolation::spline::SplineGenerator spline(lower_bound_x,lower_bound_y,0 /*lbdy*/,upper_bound_x,upper_bound_y,0/*ubdy*/);
	//std::cout << "entering for loop in util.cc..." << std::endl;
	for ( platform::Size index = 1; index <= bins_vect.size(); ++index ) {
		spline.add_known_value(bins_vect[index],potential_vect[index]);

		//std::cout << "bins_vect.size() is:  " << bins_vect.size() << std::endl;
		//std::cout << "current index is:  " << index << " current bins_vect value is:  " << bins_vect[index]
		// << " current potential_vect_value is:  " << potential_vect[index] << std::endl; }
	}


	for ( platform::Size index = 1; index <= boundary_functions.size(); ++index ) {
		std::string tag = boundary_functions[index].first;
		spline::LinearFunction function = boundary_functions[index].second;
		spline.add_boundary_function(tag,function.cutoff,function.slope,function.intercept);
	}

	return spline;


}
}
}
