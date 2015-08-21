// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   SurfaceParameters.hh
/// @brief
/// @author Robin A Thottungal (rathottungal@gmail.com)
/// @author Michael Pacella (mpacella88@gmail.com)
#ifndef INCLUDED_protocols_surface_docking_SurfaceParameters_hh
#define INCLUDED_protocols_surface_docking_SurfaceParameters_hh

// Unit headers
#include <protocols/surface_docking/SurfaceParameters.fwd.hh>

// Project headers
#include <core/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// Numeric Headers
#include <numeric/xyzVector.hh>

namespace protocols {
namespace surface_docking {

class SurfaceParameters: public utility::pointer::ReferenceCount  {

public:

	SurfaceParameters( numeric::xyzVector< core::Real > SURFA0,
		numeric::xyzVector< core::Real > SURFA1, numeric::xyzVector< core::Real > SURFA2);

	SurfaceParameters( SurfaceParameters const & src );

	~SurfaceParameters();

	SurfaceParametersOP clone() const;

	core::Vector slide_axis()
	{
		return slide_axis_;
	}

	core::Vector vecAB()
	{
		return vecAB_;
	}

	core::Vector vecAC()
	{
		return vecAC_;
	}

	void set_slide_axis(core::Vector slide_axis_in)
	{
		slide_axis_ = slide_axis_in;
	}
private:
	// Private, never called!
	SurfaceParameters() {}

	// Holds the xyz co-ordinates that form the AB and AC vector
	numeric::xyzVector< core::Real > SURFA0_;
	numeric::xyzVector< core::Real > SURFA1_;
	numeric::xyzVector< core::Real > SURFA2_;

	// AB & AC Vector
	core::Vector vecAB_;
	core::Vector vecAC_;

	// Vector along which the slide into surface needs to occur
	core::Vector slide_axis_;
};

} // namespace surface_docking
} // namespace protocols


#endif
