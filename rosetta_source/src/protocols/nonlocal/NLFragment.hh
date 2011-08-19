// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/NLFragment.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef PROTOCOLS_NONLOCAL_NLFRAGMENT_HH_
#define PROTOCOLS_NONLOCAL_NLFRAGMENT_HH_

#include <core/types.hh>

namespace protocols {
namespace nonlocal {

/// @class A read-only data structure containing known or predicted information
/// about a single residue. This includes the residue's position within the
/// sequence (1-based), its backbone torsions, and the Cartesian coordinates of
/// its Ca atom.
class NLFragment {
  typedef core::Size Size;
  typedef core::Real Real;

 public:
  /// @brief Constructs a new NLFragment with the given position, torsions, and
  /// Cartesian coordinates
  NLFragment(Size position, Real phi, Real psi, Real omega, Real x, Real y, Real z);

  // Accessors
  Size position() const;
  Real phi() const;
  Real psi() const;
  Real omega() const;
  Real x() const;
  Real y() const;
  Real z() const;

  // Operators
  bool operator==(const NLFragment& other) const;
  bool operator!=(const NLFragment& other) const;
  bool operator< (const NLFragment& other) const;

 private:
  Size pos_;

  // Torsions
  Real phi_;
  Real psi_;
  Real omega_;

  // Coordinates
  Real x_;
  Real y_;
  Real z_;
};

}  // namespace nonlocal
}  // namespace protocols

#endif  // PROTOCOLS_NONLOCAL_NLFRAGMENT_HH_
