// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/packing/PoseBalls.hh
/// @brief
/// @author

#ifndef INCLUDED_core_scoring_packing_Ball_hh
#define INCLUDED_core_scoring_packing_Ball_hh

// Project headers
#include <core/types.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <numeric/xyzVector.hh>

namespace core {
namespace scoring {
namespace packing {


class Ball : public utility::pointer::ReferenceCount {

public:

	Ball( core::Real x, core::Real y, core::Real z  , core::Real r ) : xyz_(x,y,z), r_(r) {}
	Ball( numeric::xyzVector<core::Real> const & xyz, core::Real r ) : xyz_(xyz), r_(r) {}

	inline core::Real                     const & x()      const { return xyz_.x(); }
	inline core::Real                     const & y()      const { return xyz_.y(); }
	inline core::Real                     const & z()      const { return xyz_.z(); }
	inline core::Real                     const & r()      const { return r_; }
	inline core::Real                     const & radius() const { return r_; }
	inline numeric::xyzVector<core::Real> const & xyz()    const { return xyz_; }

	inline core::Real                     & x()      { return xyz_.x(); }
	inline core::Real                     & y()      { return xyz_.y(); }
	inline core::Real                     & z()      { return xyz_.z(); }
	inline core::Real                     & r()      { return r_; }
	inline core::Real                     & radius() { return r_; }
	inline numeric::xyzVector<core::Real> & xyz()    { return xyz_; }

private:

	numeric::xyzVector<core::Real> xyz_;
	core::Real r_;

};

} // namespace packing
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_packing_Balls_HH
