// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author


#ifndef INCLUDED_protocols_viewer_triangleIterator_hh
#define INCLUDED_protocols_viewer_triangleIterator_hh


#include <numeric/xyzVector.hh>

#include <vector>

#include <ObjexxFCL/FArray3D.fwd.hh>
#include <ObjexxFCL/FArray4D.fwd.hh>


namespace protocols {
namespace viewer {

class triangleIterator {
 public:
  // Constructors and destructors
  triangleIterator(ObjexxFCL::FArray3D_float const & density, float const & threshold);
  ~triangleIterator();

  // Returns true if there is another triangle.
  bool hasNext() const;

  // Takes three parallel arrays of length three which this methods sets
  // to contain the vertix positions, normals, and velocities, respectively
  // for the next triangle triangle. The precondition is hasNext() is
  // true.
  void next(numeric::xyzVector_float vertices[3], numeric::xyzVector_float normals[3]);

 private:
  typedef std::vector<numeric::xyzVector_float> vecQueue;

  void aquireNextQueue();
  void computeGradient();
  void evalGradient(const numeric::xyzVector_float & pt, numeric::xyzVector_float & gradResult);

  // next indices
  int nextX, nextY, nextZ;

  // grids and dimensions
  float threshold;
  ObjexxFCL::FArray3D_float const * densityPtr;
  ObjexxFCL::FArray4D_float * gradPtr;
  numeric::xyzVector_int size;

  // ques for future triangles
  vecQueue vertQueue;
  vecQueue nrmlQueue;
};

}
}

#endif
