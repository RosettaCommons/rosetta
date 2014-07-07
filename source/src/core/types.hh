// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/types.hh
/// @brief  rosetta project type declarations
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_core_types_hh
#define INCLUDED_core_types_hh

// Numeric headers
#include <numeric/xyzVector.fwd.hh>

// Platform headers
#include <platform/types.hh>

// C++ headers
#include <cstddef>
#include <limits>

namespace core {
	// Integer scalars
	typedef platform::Size Size;
	typedef platform::SSize SSize;
	typedef platform::uint uint;

	// Minimum and maximum values of integer scalar type
	//static const Size SZ_MIN = std::numeric_limits<Size>::min();
	const Size SZ_MAX = std::numeric_limits<Size>::max();

#ifndef WIN32
  typedef platform::Real Real;
#else
	// Floating point precision control scalar
  #ifdef ROSETTA_FLOAT // Real == float
	  typedef  float  Real;
	#else // Real == double
	  typedef  double  Real;
	#endif
#endif

	typedef unsigned short ShortSize; // used in the conformation::Atom

	// Floating point scalars
	typedef Real Length;
	typedef Real LengthSquared;
	typedef Real Distance;
	typedef Real DistanceSquared;
	typedef Real Volume;
	typedef Real Angle;
	typedef Real Trig; // Trigonometric values of angles
	typedef Real Mass;
	typedef Real Charge;
	typedef Real Energy;
	typedef Real EnergyDerivative;

	// Double precision slows the packer considerably
	typedef float PackerEnergy;

	// Floating point arrays
	typedef  numeric::xyzVector< Length >  PointPosition;
	typedef  numeric::xyzVector< Length >  Vector;
	typedef  numeric::xyzVector< EnergyDerivative >  EnergyGradient;
}  // namespace core

#endif // INCLUDED_core_types_HH
