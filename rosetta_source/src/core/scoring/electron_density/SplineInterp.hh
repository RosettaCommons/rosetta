// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/electron_density/SplineInterp.hh
/// @brief  3D spline interpolation methods.  Based on implementation by Philippe Thevenaz, see comments below.


//  Fast 3d spline interpolation routines
//
//  Based on implementation by Philippe Thevenaz.
//  See http://bigwww.epfl.ch/thevenaz/interpolation/ for details.

#ifndef INCLUDED_core_scoring_electron_density_SplineInterp_hh
#define INCLUDED_core_scoring_electron_density_SplineInterp_hh

namespace core {
namespace scoring {
namespace electron_density {
namespace SplineInterp {

double interp3(double *Bcoeff, int dims[3], double X[3], int degree);
int grad3(double grad[3], double *Bcoeff, int dims[3], double X[3], int degree);
int compute_coefficients(double *data, int dims[3], int degree);

}
}
}
}

#endif
