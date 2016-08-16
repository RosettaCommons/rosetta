// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/numeric/XyzVector.string.hh
/// @author Sam DeLuca

#ifndef INCLUDED_numeric_XyzVector_string_hh
#define INCLUDED_numeric_XyzVector_string_hh

#include <numeric/xyzVector.hh>

#include <numeric/types.hh>

#include <sstream>
#include <string>
namespace numeric
{

template <typename T>
std::string truncate_and_serialize_xyz_vector(xyzVector<T> vector, Real precision)
{
	std::stringstream vector_stream;

	//set up proper formatting
	vector_stream.setf(std::ios::fixed,std::ios::floatfield);
	vector_stream.precision(static_cast<std::streamsize>(precision));

	//the output of this is not meant to
	vector_stream << vector.x() << vector.y() << vector.z();
	return vector_stream.str();
}

}
#endif /* XYZVECTOR_STRING_HH_ */
