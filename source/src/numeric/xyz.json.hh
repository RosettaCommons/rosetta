// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/numeric/xyz.json.hh
/// @author Sam DeLuca

#ifndef INCLUDED_numeric_xyz_json_hh
#define INCLUDED_numeric_xyz_json_hh

#include <numeric/xyzVector.hh>

#include <utility/tools/make_vector.hh>
#include <utility/json_spirit/json_spirit_value.h>


namespace numeric {


/// @brief Convert vector to a json_spirit Value
/// @note Format is a list in the form [x,y,z]
template<typename T>
inline
utility::json_spirit::Value serialize(xyzVector<T> coords)
{
	utility::json_spirit::Value x(coords.x());
	utility::json_spirit::Value y(coords.y());
	utility::json_spirit::Value z(coords.z());

	return utility::json_spirit::Value(utility::tools::make_vector(x,y,z));
}

template<typename T>
inline
xyzVector<T> deserialize(utility::json_spirit::mArray data)
{
	xyzVector<T> coords;
	coords.x(data[0].get_value<T>());
	coords.y(data[1].get_value<T>());
	coords.z(data[2].get_value<T>());
	return coords;
}


}


#endif /* XYZ_JSON_HH_ */
