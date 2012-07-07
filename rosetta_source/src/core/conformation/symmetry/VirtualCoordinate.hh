// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @brief  Symmetry data container
/// @file   core/conformation/symmetry/SymmData.hh
/// @author Ingemar Andre


#ifndef INCLUDED_core_conformation_symmetry_VirtualCoordinate_hh
#define INCLUDED_core_conformation_symmetry_VirtualCoordinate_hh

// Utility headers
#include <core/conformation/symmetry/VirtualCoordinate.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <core/types.hh>
#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>
#include <vector>

namespace core {
namespace conformation {
namespace symmetry {

class VirtualCoordinate {

	public:

	VirtualCoordinate();

	/// @brief copy constructor
	VirtualCoordinate( VirtualCoordinate const & src );

	VirtualCoordinate(
		numeric::xyzVector< core::Real> axis_x,
		numeric::xyzVector< core::Real> axis_y,
		numeric::xyzVector< core::Real> axis_origin
	);

	VirtualCoordinate &
  operator=( VirtualCoordinate const & src );

	~VirtualCoordinate();

	// @details accessor functions
	numeric::xyzVector< core::Real> &
	get_x();

	numeric::xyzVector< core::Real> &
	get_y();

	numeric::xyzVector< core::Real> &
	get_origin();

	numeric::xyzVector< core::Real> const &
	get_x() const;

	numeric::xyzVector< core::Real> const &
	get_y() const;

	numeric::xyzVector< core::Real> const &
	get_origin() const;

	void
	add_coordinate_from_string(
		utility::vector1< std::string > coords,
		core::Size coord_start=2
	);


	private:

	numeric::xyzVector< core::Real> axis_x_; // store unit vector for X
	numeric::xyzVector< core::Real> axis_y_; // store unit vector for Y
	numeric::xyzVector< core::Real> axis_origin_; // store origin for coordinate system
};

} // symmetry
} // conformation
} // core
#endif
