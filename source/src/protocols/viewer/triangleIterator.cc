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

#include <protocols/viewer/triangleIterator.hh>
#include <protocols/viewer/marchingCubes.hh>

// Utility headers
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.hh>
#include <utility/exit.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>

#include <utility/vector1.hh>
#include <ObjexxFCL/FArray4D.hh>


using namespace ObjexxFCL;
using namespace numeric;

namespace protocols {
namespace viewer {


triangleIterator::triangleIterator(FArray3D_float const & density, float const & threshold) {
  densityPtr = &density;
  size = numeric::xyzVector<int>(density.size1(), density.size2(), density.size3());

  gradPtr = new FArray4D_float(3, size[0], size[1], size[2]);
  computeGradient();
  this->threshold = threshold;

  // This will start us off looking for triangles at (0, 0, 0)
  nextX = nextY = 0;
  nextZ = -1;

  aquireNextQueue();
}

triangleIterator::~triangleIterator() {
  delete gradPtr; gradPtr = NULL;
}

bool triangleIterator::hasNext() const {
  return !vertQueue.empty();
}

// Takes three parallel arrays of length three which this methods sets
// to contain the vertix positions, normals, and velocities, respectively
// for the next triangle. The precondition is hasNext() is true.
void triangleIterator::next(xyzVector_float vertices[3], xyzVector_float normals[3]) {
  assert(hasNext());

  for (int v = 0 ; v < 3 ; v++) {
    vertices[v] = vertQueue.back();
    normals[v] = nrmlQueue.back();

    vertQueue.pop_back();
    nrmlQueue.pop_back();
  }

  if (vertQueue.empty())
    aquireNextQueue();
}


void triangleIterator::aquireNextQueue() {
  const FArray3D_float & density(*(this->densityPtr));

  // This function should only be called to fill up the queue.
  assert(vertQueue.size() == 0);

  // We always start at (x, y, z) + (0, 0, 1)
  // so that we don't process the same triangle twice.
  bool starting = true;
  nextZ++;

  // Look for a 2x2x2 cube with triangles.
  for (int x = starting ? nextX : 0 ; x < size[0] -1 ; x++) {
    for (int y = starting ? nextY : 0 ; y < size[1] - 1; y++) {
      for (int z = starting ? nextZ : 0 ; z < size[2] - 1; z++) {
		// Load this baby into a bitfield to see if we have any triangels.
		int bitfield = 0;
		for (int i = 1 ; i <= 8 ; i++) {
		  int vx = x + VERTEX_OFF[i][0];
		  int vy = y + VERTEX_OFF[i][1];
		  int vz = z + VERTEX_OFF[i][2];
		  if (density(vx + 1, vy + 1, vz + 1) > threshold) {
			bitfield |= (1 << (i - 1));
	  }
	}

	// Load up the triangles
	for (int j = 0 ; POLY_CASES[bitfield][j] != 0 ; j++) {
	  for (int k = 0 ; k < 3 ; k++) {
	    // look up the edge
	    int edgeIndex = POLY_CASES[bitfield][j++];

	    // find the adjacent vertices
	    int v0 = EDGE_NGHBRS[edgeIndex][0];
	    int x0 = x + VERTEX_OFF[v0][0] + 1;
	    int y0 = y + VERTEX_OFF[v0][1] + 1;
	    int z0 = z + VERTEX_OFF[v0][2] + 1;

	    int v1 = EDGE_NGHBRS[edgeIndex][1];
	    int x1 = x + VERTEX_OFF[v1][0] + 1;
	    int y1 = y + VERTEX_OFF[v1][1] + 1;
	    int z1 = z + VERTEX_OFF[v1][2] + 1;

	    // look up the level set values
	    float phi0 = density(x0, y0, z0) - threshold;
	    float phi1 = density(x1, y1, z1) - threshold;

	    // compute the distance to front
	    xyzVector_float p0(x0, y0, z0);
	    xyzVector_float p1(x1, y1, z1);
	    float dist = phi0 / (phi0 - phi1);

	    // compute the position, normal and gradient
	    xyzVector_float vert = p0 * (1.0f - dist) + p1 * dist;
	    xyzVector_float nrml;
	    evalGradient(vert, nrml);

	    // stick them in queues
	    vertQueue.push_back(vert);
	    nrmlQueue.push_back(-nrml);
	  }
	}

	// If we managaed to fill up the queue then we're done.
	if (!vertQueue.empty()) {
	  nextX = x;
	  nextY = y;
	  nextZ = z;
	  return;
	}
	starting = false;
      }
      starting = false;
    }
    starting = false;
  }
}

void triangleIterator::computeGradient() {
/*
  FArray4D_float & grad(*gradPtr);
  const FArray3D_float & density(*densityPtr);

  grad = 0.0;
  // this is proportional to the gradient, but not equal because we don't know the grid spacing
  for (int i = 2 ; i < size[0] ; i++) {
    for (int j = 2 ; j < size[1] ; j++) {
      for (int k = 2 ; k < size[2] ; k++) {
		grad(1, i, j, k) = density(i+1, j, k) - density(i-1, j, k);
		grad(2, i, j, k) = density(i, j+1, k) - density(i, j-1, k);
		grad(3, i, j, k) = density(i, j, k+1) - density(i, j, k-1);
      }
    }
  }
*/
}

void triangleIterator::evalGradient(const xyzVector_float & pt, xyzVector_float & gradResult) {
	// get the gradient of the density
	const core::scoring::electron_density::ElectronDensity& edm = core::scoring::electron_density::getDensityMap();
	numeric::xyzVector< core::Real > pt_d = numeric::xyzVector< core::Real >( pt[0] , pt[1] , pt[2] );
	numeric::xyzVector< core::Real > gx = edm.dens_grad( pt_d );
	gradResult = xyzVector_float( gx[0] , gx[1], gx[2] );

/*
  // hopefully this linearly interpolates the gradient (with normalization)

  FArray4D_float & gradient(*gradPtr);

  float x = pt(1);
  float y = pt(2);
  float z = pt(3);
  if (x < 2.0) x = 2.0; else if (x > (size[0] - 1.0)) x = size[0] - 1.0;
  if (y < 2.0) y = 2.0; else if (y > (size[1] - 1.0)) y = size[1] - 1.0;
  if (z < 2.0) z = 2.0; else if (z > (size[2] - 1.0)) z = size[2] - 1.0;

  int x1 = (int) x;
  int y1 = (int) y;
  int z1 = (int) z; // x;
  int x2 = x1 + 1;
  int y2 = y1 + 1;
  int z2 = z1 + 1;

  float coefX2 = x - x1;
  float coefY2 = y - y1;
  float coefZ2 = z - z1;
  float coefX1 = 1.0 - coefX2;
  float coefY1 = 1.0 - coefY2;
  float coefZ1 = 1.0 - coefZ2;

  gradResult(1) = gradResult(2) = gradResult(3) = 0.0;
  for (int dim = 1 ; dim <= 3 ; dim++) {
    gradResult(dim) += gradient(dim, x1, y1, z1) * coefX1 * coefY1 * coefZ1;
    gradResult(dim) += gradient(dim, x2, y1, z1) * coefX2 * coefY1 * coefZ1;
    gradResult(dim) += gradient(dim, x1, y2, z1) * coefX1 * coefY2 * coefZ1;
    gradResult(dim) += gradient(dim, x2, y2, z1) * coefX2 * coefY2 * coefZ1;
    gradResult(dim) += gradient(dim, x1, y1, z2) * coefX1 * coefY1 * coefZ2;
    gradResult(dim) += gradient(dim, x2, y1, z2) * coefX2 * coefY1 * coefZ2;
    gradResult(dim) += gradient(dim, x1, y2, z2) * coefX1 * coefY2 * coefZ2;
    gradResult(dim) += gradient(dim, x2, y2, z2) * coefX2 * coefY2 * coefZ2;
  }
  gradResult.normalize_any();
*/
}


}
}
