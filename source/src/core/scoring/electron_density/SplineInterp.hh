// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/electron_density/SplineInterp.hh
/// @brief  3D spline interpolation methods.  Based on implementation by Philippe Thevenaz, see comments below.


//  Fast 3d spline interpolation routines
//
//  Based on implementation by Philippe Thevenaz.
//  See http://bigwww.epfl.ch/thevenaz/interpolation/ for details.

#ifndef INCLUDED_core_scoring_electron_density_SplineInterp_hh
#define INCLUDED_core_scoring_electron_density_SplineInterp_hh

#include <cmath>

namespace core {
namespace scoring {
namespace electron_density {
namespace SplineInterp {

// 3d --- periodic or mirrored boundaries on all dimensions
double interp3(const double *Bcoeff, const int dims[3], double X[3], bool mirrored=false);
int grad3(double grad[3], const double *Bcoeff, const int dims[3], double X[3], bool mirrored=false);
int compute_coefficients3(double *data, int dims[3], bool mirrored=false);

//
template <class S>
double interp3(const S *Bcoeff, const int dims[3], double X[3], bool mirrored=false)
{
	double wt[3][4];
	double value;
	double sum_k;
	int idx[3][4];
	int i,j,k, dim, dims2;

	// interpolation indexes
	for ( dim=0; dim<3; dim++ ) {
		int pt = (int)floor(X[dim] - (3-1) / 2.0);
		for ( i = 0L; i <= 3; i++ ) {
			idx[dim][i] = pt++;
		}

		// interpolation weights
		double w = X[dim] - (double)idx[dim][1];
		wt[dim][3] = (1.0 / 6.0) * w * w * w;
		wt[dim][0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - wt[dim][3];
		wt[dim][2] = w + wt[dim][0] - 2.0 * wt[dim][3];
		wt[dim][1] = 1.0 - wt[dim][0] - wt[dim][2] - wt[dim][3];

		if ( mirrored ) {
			// _mirror_ boundary conditions
			for ( i = 0L; i <= 3; i++ ) {
				if ( dims[dim] == 1 ) {
					idx[dim][i] = 0;
				} else {
					dims2 = 2*dims[dim]-2;
					if ( idx[dim][i]<0 ) {
						idx[dim][i] = (-idx[dim][i] - dims2 * (-idx[dim][i]/dims2));
					} else {
						idx[dim][i] = (idx[dim][i] - dims2 * (idx[dim][i]/dims2));
					}
					if ( idx[dim][i] >= dims[dim] ) {
						idx[dim][i] = dims2 - idx[dim][i];
					}
				}
			}
		} else {
			// _periodic_ boundary conditions
			for ( i = 0L; i <= 3; i++ ) {
				if ( dims[dim] == 1 ) {
					idx[dim][i] = 0;
				} else {
					idx[dim][i] = idx[dim][i] % dims[dim];
				}
				if ( idx[dim][i] < 0 ) {
					idx[dim][i] += dims[dim];
				}
			}
		}
	}

	value = 0.0;
	for ( i = 0; i <= 3; i++ ) {  // x
		double sum_jk = 0.0;
		for ( j = 0; j <= 3; j++ ) {  // y
			sum_k = 0.0;
			for ( k = 0; k <= 3; k++ ) {  // z
				sum_k += wt[2][k] * S(Bcoeff[idx[0][i]*dims[1]*dims[2] + idx[1][j]*dims[2] + idx[2][k]]);
			}
			sum_jk += wt[1][j] * sum_k;
		}
		value += wt[0][i] * sum_jk;
	}

	return(value);
}

template <class S>
int grad3(double grad[3], const S *Bcoeff, const int dims[3], double X[3], bool mirrored=false)
{
	double wt[3][4];
	double w;
	double sum_k, sum_jk;
	int idx[3][4];
	int i,j,k, pt, dim, gradDim, dims2;

	// compute interpolation indexes
	for ( gradDim=0; gradDim<3; gradDim++ ) {
		for ( dim=0; dim<3; dim++ ) {
			pt = (int)floor(X[dim] - (3-1) / 2.0);
			for ( i = 0L; i <= 3; i++ ) {
				idx[dim][i] = pt++;
			}

			// compute the interpolation weights
			if ( dim == gradDim ) {
				w = X[dim] - (double)idx[dim][1];
				wt[dim][3] = (1.0 / 2.0) * w * w;
				wt[dim][0] = (w - 1.0/2.0) - wt[dim][3];
				wt[dim][2] = 1.0 + wt[dim][0] - 2.0 * wt[dim][3];
				wt[dim][1] = - wt[dim][0] - wt[dim][2] - wt[dim][3];
			} else {
				w = X[dim] - (double)idx[dim][1];
				wt[dim][3] = (1.0 / 6.0) * w * w * w;
				wt[dim][0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - wt[dim][3];
				wt[dim][2] = w + wt[dim][0] - 2.0 * wt[dim][3];
				wt[dim][1] = 1.0 - wt[dim][0] - wt[dim][2] - wt[dim][3];
			}

			if ( mirrored ) {
				// _mirror_ boundary conditions
				for ( i = 0L; i <= 3; i++ ) {
					if ( dims[dim] == 1 ) {
						idx[dim][i] = 0;
					} else {
						dims2 = 2*dims[dim]-2;
						if ( idx[dim][i]<0 ) {
							idx[dim][i] = (-idx[dim][i] - dims2 * (-idx[dim][i]/dims2));
						} else {
							idx[dim][i] = (idx[dim][i] - dims2 * (idx[dim][i]/dims2));
						}
						if ( idx[dim][i] >= dims[dim] ) {
							idx[dim][i] = dims2 - idx[dim][i];
						}
					}
				}
			} else {
				// _periodic_ boundary conditions
				for ( i = 0L; i <= 3; i++ ) {
					if ( dims[dim] == 1 ) {
						idx[dim][i] = 0;
					} else {
						idx[dim][i] = idx[dim][i] % dims[dim];
					}
					if ( idx[dim][i] < 0 ) {
						idx[dim][i] += dims[dim];
					}
				}
			}
		}

		grad[gradDim] = 0.0;
		for ( i = 0; i <= 3; i++ ) {  // x
			sum_jk = 0.0;
			for ( j = 0; j <= 3; j++ ) {  // y
				sum_k = 0.0;
				for ( k = 0; k <= 3; k++ ) {  // z
					sum_k += wt[2][k] * S(Bcoeff[idx[0][i]*dims[1]*dims[2] + idx[1][j]*dims[2] + idx[2][k]]);
				}
				sum_jk += wt[1][j] * sum_k;
			}
			grad[gradDim] += wt[0][i] * sum_jk;
		}
	}

	return(0);
}

// 4d --- periodic boundaries on dimensions 2-4; mirror boundaries on dimension 1
double interp4(const double *Bcoeff, const int dims[4], double X[4]);
int grad4(double grad[4], const double *Bcoeff, const int dims[4], double X[4]);
int compute_coefficients4(double *data, int dims[3]);

}
}
}
}

#endif
