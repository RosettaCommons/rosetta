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
#include <numeric/model_quality/rms.hh>

#include <utility/vector1.hh>
#include <basic/Tracer.hh>

#include <protocols/seeded_abinitio/util.hh>
#include <protocols/toolbox/superimpose.hh>


namespace protocols {
namespace seeded_abinitio {

static basic::Tracer TR( "protocols.seeded_abinitio.util" );

void
superposition_transform(
	utility::vector1< numeric::xyzVector< numeric::Real > > & init_coords,
	utility::vector1< numeric::xyzVector< numeric::Real > > & ref_coords,
	numeric::xyzMatrix< numeric::Real > & rotation,
	numeric::xyzVector< numeric::Real > & toCenter,
	numeric::xyzVector< numeric::Real > & toFitCenter)
{
	core::Size count = init_coords.size();
	runtime_assert(count == ref_coords.size());

	ObjexxFCL::FArray2D< numeric::Real > init_fa( 3, init_coords.size());
	ObjexxFCL::FArray2D< numeric::Real > ref_fa( 3, ref_coords.size());

	vector_vector_to_FArray2(init_coords, init_fa);
	vector_vector_to_FArray2(ref_coords, ref_fa);

	ObjexxFCL::FArray1D< numeric::Real > weights_fa(count, 1);

	protocols::toolbox::superposition_transform(
		count,
		weights_fa,
		ref_fa,
		init_fa,
		rotation,
		toCenter,
		toFitCenter);
}

template <typename T>
void
vector_vector_to_FArray2(
	utility::vector1< numeric::xyzVector< T > > & from,
	ObjexxFCL::FArray2D< T > & to)
{
	int count = from.size();

	/*
	runtime_assert(to.I1() == 3);
	runtime_assert(to.I2() == count);
	*/

	for ( int i = 1; i <= count; i++ ) {
		for ( int j = 1; j <= 3; j++ ) {
			to(j,i) = from[i](j);
		}
	}
}

}
}
