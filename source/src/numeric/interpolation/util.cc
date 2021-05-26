// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/numeric/interpolation/util.cc
/// @author Sam DeLuca
/// @author Stephanie Hirst

#include <numeric/interpolation/util.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <iostream>
#include <string>
#include <tuple>

#include <boost/algorithm/string/trim.hpp>                  // for trim
#include <utility/vector1.hh>


namespace numeric {
namespace interpolation {


/// @details Easier make a spline wrapper from code
spline::SplineGenerator make_spline(
	utility::vector1<platform::Real> const & bins_vect,
	utility::vector1<platform::Real> const & potential_vect,
	platform::Real const bin_size,
	utility::vector1<std::tuple<std::string, platform::Real, platform::Real, platform::Real>> const & boundary_functions
) {

	if ( bins_vect.size() != potential_vect.size() ) {
		std::stringstream ss;
		ss << "The x and y axis lines in the given spline do not have the same number of entries\n";
		ss << "x_axis: " << bins_vect << "\n";
		ss << "y_axis: " << potential_vect << "\n";
		utility_exit_with_message(ss.str());
	}

	if ( bins_vect.empty() || potential_vect.empty() ) {
		std::stringstream ss;
		ss << "The x and/or y axis lines in the given spline are empty\n";
		ss << "x_axis: " << bins_vect << "\n";
		ss << "y_axis: " << potential_vect << "\n";
		utility_exit_with_message(ss.str());
	}

	platform::Real lower_bound_y(potential_vect.front());
	platform::Real upper_bound_y(potential_vect.back());
	platform::Real lower_bound_x(bins_vect.front() - (0.5*bin_size));
	platform::Real upper_bound_x(bins_vect.back() + (0.5*bin_size));

	numeric::interpolation::spline::SplineGenerator spline(lower_bound_x,lower_bound_y,0 /*lbdy*/,upper_bound_x,upper_bound_y,0/*ubdy*/);

	for ( platform::Size index = 1; index <= bins_vect.size(); ++index ) {
		spline.add_known_value(bins_vect[index], potential_vect[index]);
	}

	for ( auto const & current_boundary_fn : boundary_functions ) {
		spline.add_boundary_function(
			std::get<0>(current_boundary_fn),  // TAG
			std::get<1>(current_boundary_fn),  // cutoff
			std::get<2>(current_boundary_fn),  // slope
			std::get<3>(current_boundary_fn));  // intercept
	}
	return spline;
}


/// @details read in a file, read out a spline.
///The file should be tab separated, and have two lines
///The first field of one line should be "x_axis", the next fields should be the x values of the points for the spline
///The first field of the other line should be "y_axis", the next fields should be the y values for the points on the spline
///For an example, see "scoring/constraints/epr_distance_potential.histogram"
spline::SplineGenerator spline_from_file(std::string const &  filename, platform::Real const bin_size)
{
	utility::io::izstream potential_file;
	potential_file.open(filename.c_str());

	utility::vector1<platform::Real> bins_vect;
	utility::vector1<platform::Real> potential_vect;
	utility::vector1<std::tuple<std::string, platform::Real, platform::Real, platform::Real>> boundary_functions;

	std::string line;

	while ( getline(potential_file, line) ) {
		boost::trim(line);
		if ( ( line.size() == 0 ) || ( line[ 0 ] == '#' ) ) {
			continue;
		}
		utility::vector1<std::string>  split_fields(utility::string_split(line,'\t'));
		if ( split_fields[1] == "x_axis" ) {
			for ( platform::Size count = 2; count <= split_fields.size(); ++count ) {
				bins_vect.push_back(utility::from_string(split_fields[count],platform::Real(0.0)));
			}
		} else if ( split_fields[1] == "y_axis" ) {
			for ( platform::Size count = 2; count <= split_fields.size(); ++count ) {
				potential_vect.push_back(utility::from_string(split_fields[count],platform::Real(0.0)));
			}
		}
		if ( split_fields[1] == "lb_function" || split_fields[1] == "ub_function" ) {
			platform::Real cutoff= utility::from_string(split_fields[2],platform::Real(0.0));
			platform::Real slope = utility::from_string(split_fields[3],platform::Real(0.0));
			platform::Real intercept = utility::from_string(split_fields[4],platform::Real(0.0));
			boundary_functions.push_back(std::make_tuple(split_fields[1], cutoff, slope, intercept));
		}
	}
	potential_file.close();

	if ( bins_vect.size() != potential_vect.size() ) {
		std::stringstream ss;
		ss << "The x and y axis lines in the given spline do not have the same number of entries\n";
		ss << "Filename: " << filename << "\n";
		ss << "x_axis: " << bins_vect << "\n";
		ss << "y_axis: " << potential_vect << "\n";
		utility_exit_with_message(ss.str());
	}

	return make_spline(bins_vect, potential_vect, bin_size, boundary_functions);
}

/// @brief A little helper enum to determine what I'm reading
/// right now in the spline_from_stream() function.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
enum class SplineFromStreamMode {
	none = 0,
	x_axis,
	y_axis,
	lb_or_ub_function
};

/// @brief Initialize a SplineGenerator from a std::istream.  Note: only reads from the currenet line of the istream.
/// @details Assumes that the stream contains only the following:
/// "x_axis" followed by N values -- mandatory.
/// "y_axis" followed by N values -- mandatory.  N must match "xaxis".
/// "lb_function" [cutoff] [slope] [intercept] -- optional.
/// "ub_function" [cutoff] [slope] [intercept] -- optional.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
spline::SplineGenerator
spline_from_stream(
	std::istream & iss,
	platform::Real const bin_size
) {
	std::string curline;
	std::getline( iss, curline );
	iss.unget(); //Go back one character, to before the carriage return.
	std::istringstream ss( curline );

	std::string const errmsg( "Could not parse spline definition \"" + curline + "\"." );
	SplineFromStreamMode mode = SplineFromStreamMode::none;
	utility::vector1<platform::Real > xpoints, ypoints;
	utility::vector1<std::tuple<std::string, platform::Real, platform::Real, platform::Real>> boundary_functions;

	do {
		std::string buffer;
		ss >> buffer;
		if ( buffer == "x_axis" ) {
			mode = SplineFromStreamMode::x_axis;
			continue;
		} else if ( buffer == "y_axis" ) {
			mode = SplineFromStreamMode::y_axis;
			continue;
		} else if ( buffer == "lb_function" || buffer == "ub_function" ) {
			//mode = SplineFromStreamMode::lb_or_ub_function; //Do not delete.  This bothers the clang-analysis warnings (since the stored value is not read), but should be uncommented if mode is ever needed before the end of the scope.
			platform::Real cutoff, slope, intercept;
			ss >> cutoff;
			runtime_assert_string_msg( !(ss.bad() || ss.fail() || ss.eof()), errmsg + "  Failed when trying to get cutoff for " + buffer + "." );
			ss >> slope;
			runtime_assert_string_msg( !(ss.bad() || ss.fail() || ss.eof()), errmsg + "  Failed when trying to get slope for " + buffer + "." );
			ss >> intercept;
			runtime_assert_string_msg( !(ss.bad() || ss.fail()), errmsg + "  Failed when trying to get intercept for " + buffer + "." );
			boundary_functions.push_back( std::make_tuple( buffer, cutoff, slope, intercept ) );
			mode = SplineFromStreamMode::none;
			continue;
		}
		//Parse entries in x- or y-axis lines.
		runtime_assert_string_msg( mode == SplineFromStreamMode::x_axis || mode == SplineFromStreamMode::y_axis, errmsg + "  Unrecognized entry \"" + buffer + "\"." );
		std::istringstream ssentry( buffer );
		platform::Real val;
		ssentry >> val;
		runtime_assert_string_msg( !(ssentry.bad() || ssentry.fail()), errmsg + "  Could not parse \"" + buffer + "\" as a floating-point number." );
		if ( mode == SplineFromStreamMode::x_axis ) {
			xpoints.push_back(val);
		} else {
			debug_assert( mode == SplineFromStreamMode::y_axis ); //Should be guaranteed true.
			ypoints.push_back(val);
		}
	} while( !ss.eof() );

	runtime_assert_string_msg( xpoints.size() > 0, errmsg + "  No x_axis values were found." );
	runtime_assert_string_msg( ypoints.size() > 0, errmsg + "  No y_axis values were found." );
	runtime_assert_string_msg( xpoints.size() == ypoints.size(), errmsg + "  Mismatch in number of x_axis and y_axis values." );

	return make_spline(xpoints, ypoints, bin_size, boundary_functions);
}

}
}
