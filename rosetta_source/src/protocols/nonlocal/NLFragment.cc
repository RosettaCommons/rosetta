// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/NLFragment.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <protocols/nonlocal/NLFragment.hh>

// Project headers
#include <core/types.hh>

namespace protocols {
namespace nonlocal {

typedef core::Size Size;
typedef core::Real Real;

NLFragment::NLFragment(Size position,
                       Real phi,
                       Real psi,
                       Real omega,
                       Real x,
                       Real y,
                       Real z)
    : pos_(position), phi_(phi), psi_(psi), omega_(omega),
      x_(x), y_(y), z_(z) {}

Size NLFragment::position() const {
  return pos_;
}

Real NLFragment::phi() const {
  return phi_;
}

Real NLFragment::psi() const {
  return psi_;
}

Real NLFragment::omega() const {
  return omega_;
}

Real NLFragment::x() const {
  return x_;
}

Real NLFragment::y() const {
  return y_;
}

Real NLFragment::z() const {
  return z_;
}

bool NLFragment::operator==(const NLFragment& other) const {
  return
      x() == other.x() &&
      y() == other.y() &&
      z() == other.z() &&
      phi() == other.phi() &&
      psi() == other.psi() &&
      omega() == other.omega() &&
      position() == other.position();
}

bool NLFragment::operator!=(const NLFragment& other) const {
  return !(*this == other);
}

bool NLFragment::operator<(const NLFragment& other) const {
  return position() < other.position();
}

}  // namespace nonlocal
}  // namespace protocols
