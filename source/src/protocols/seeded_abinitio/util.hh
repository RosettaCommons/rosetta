// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/RosettaScripts/util.hh
/// @brief Utility functions useful in Seeded Abinitio.
/// @author Alex Ford (fordas@uw.edu)
//
#ifndef INCLUDED_protocols_seeded_abinitio_util_hh
#define INCLUDED_protocols_seeded_abinitio_util_hh

// Unit headers

// Project Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.fwd.hh>

// C++ headers
#include <string>
#include <set>

#include <numeric/types.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>

#include <ObjexxFCL/FArray.all.hh>

#include <utility/vector1.hh>
#include <basic/Tracer.hh>

namespace protocols {
namespace seeded_abinitio {


/* @brief Calculates superposition transform from coords to ref_coords.
*
* Modifies ref_coords and coords, moving coords into superposition.
*/
void
superposition_transform(
	utility::vector1< numeric::xyzVector< numeric::Real > > & init_coords,
	utility::vector1< numeric::xyzVector< numeric::Real > > & ref_coords,
	numeric::xyzMatrix< numeric::Real > & rotation,
	numeric::xyzVector< numeric::Real > & toCenter,
	numeric::xyzVector< numeric::Real > & toFitCenter);

template <typename T>
void
vector_vector_to_FArray2(
	utility::vector1< numeric::xyzVector< T > > & from,
	ObjexxFCL::FArray2D< T > & to);


}
}


#endif
