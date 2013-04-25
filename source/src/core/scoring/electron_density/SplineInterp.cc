// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/electron_density/SplineInterp.cc
/// @brief  3D spline interpolation methods.  Based on implementation by Philippe Thevenaz, see comments below.


//  Fast 3d spline interpolation routines
//   --> pretty specific to density map interpolation (3d only, periodic boundaries)
//       so it lives in core/scoring/methods/electron_density
//
//  Based on implementation by Philippe Thevenaz.
//  See http://bigwww.epfl.ch/thevenaz/interpolation/ for details.

#include <core/scoring/electron_density/SplineInterp.hh>

#include <cmath>
#include <float.h>  // for DBL_EPSILON
#include <vector>
#include <iostream>


namespace core {
namespace scoring {
namespace electron_density {
namespace SplineInterp {

void put_line(double* data, int dim, int x1, int x2, double line[], int dims[]) {
	int i, inc;
	double* ptr;

	if (dim == 0) {          // x1 == y, x2 == z
		ptr = &data[x1*dims[2] + x2];
		inc = dims[1]*dims[2];
	} else if (dim == 1) {   // x1 == x, x2 == z
		ptr = &data[x1*dims[1]*dims[2] + x2];
		inc = dims[2];
	} else {                 // x1 == x, x2 == y
		ptr = &data[x1*dims[1]*dims[2] + x2*dims[2]];
		inc = 1;
	}

	for (i=0; i<dims[dim]; i++) {
		*ptr = line[i] ;
		ptr = &ptr[inc];
	}
}


void get_line(double* data, int dim, int x1, int x2, double line[], int dims[]) {
	int i, inc;
	double* ptr;

	if (dim == 0) {          // x1 == y, x2 == z
		ptr = &data[x1*dims[2] + x2];
		inc = dims[1]*dims[2];
	} else if (dim == 1) {   // x1 == x, x2 == z
		ptr = &data[x1*dims[1]*dims[2] + x2];
		inc = dims[2];
	} else {                 // x1 == x, x2 == y
		ptr = &data[x1*dims[1]*dims[2] + x2*dims[2]];
		inc = 1;
	}

	for (i=0; i<dims[dim]; i++) {
		line[i] = *ptr;
		ptr = &ptr[inc];
	}
}

static double	InitialCausalCoefficient
				(
					double	c[],		// coefficients
					long	DataLength,	// number of coefficients
					double	z,			// actual pole
					double	Tolerance	// admissible relative error
				)
{

	double	Sum, zn;
	long	n, Horizon;

	// this initialization corresponds to mirror boundaries
	// modified FPD -- periodic boundaries
	Horizon = DataLength;
	if (Tolerance > 0.0) {
		Horizon = (long)ceil(log(Tolerance) / log(fabs(z)));
	}
	if (Horizon < DataLength) {
		// accelerated loop
		zn = z;
		Sum = c[0];
		for (n = 1L; n < Horizon; n++) {
			Sum += zn * c[DataLength - n];
			zn *= z;
		}
		return(Sum);
	}
	else {
		// full loop
		zn = z;
		Sum = c[0];
		for (n = 1L; n < DataLength; n++) {
			Sum += zn * c[DataLength - n];
			zn *= z;
		}
		return(Sum / (1.0 - zn));
	}
}


static double	InitialAntiCausalCoefficient
				(
					double	c[],		// coefficients
					long	DataLength,	// number of samples or coefficients
					double	z,			// actual pole
					double	Tolerance	// admissible relative error
				)
{
	// this initialization corresponds to mirror boundaries
	// return((z / (z * z - 1.0)) * (z * c[DataLength - 2L] + c[DataLength - 1L]));
	// modified FPD -- periodic boundaries
	double Sum, zn;
	int n, Horizon;

	Horizon = DataLength;
	if (Tolerance > 0.0) {
		Horizon = (long)ceil(log(Tolerance) / log(fabs(z)));
	}
	if (Horizon < DataLength) {
		// accelerated loop
		zn = z;
		Sum = c[DataLength-1];
		for (n = 0L; n < Horizon; n++) {
			Sum += zn * c[n];
			zn *= z;
		}
		return(-z*Sum);
	}
	else {
		// full loop
		zn = z;
		Sum = c[DataLength-1];
		for (n = 0L; n < DataLength-1; n++) {
			Sum += zn * c[n];
			zn *= z;
		}
		return(-z*Sum / (1.0 - zn));
	}
}



void ConvertToInterpolationCoefficients
				(
					double	c[],		// input samples --> output coefficients
					long	DataLength,	// number of samples or coefficients
					double	z[],		// poles
					long	NbPoles,	// number of poles
					double	Tolerance	// admissible relative error
				)
{

	double	Lambda = 1.0;
	long	n, k;

	// special case required by mirror boundaries
	if (DataLength == 1L) {
		return;
	}
	// compute the overall gain
	for (k = 0L; k < NbPoles; k++) {
		Lambda = Lambda * (1.0 - z[k]) * (1.0 - 1.0 / z[k]);
	}
	// apply the gain
	for (n = 0L; n < DataLength; n++) {
		c[n] *= Lambda;
	}
	// loop over all poles
	for (k = 0L; k < NbPoles; k++) {
		// causal initialization
		c[0] = InitialCausalCoefficient(c, DataLength, z[k], Tolerance);
		// causal recursion
		for (n = 1L; n < DataLength; n++) {
			c[n] += z[k] * c[n - 1L];
		}
		// anticausal initialization
		c[DataLength - 1L] = InitialAntiCausalCoefficient(c, DataLength, z[k], Tolerance);
		// anticausal recursion
		for (n = DataLength - 2L; 0 <= n; n--) {
			c[n] = z[k] * (c[n + 1L] - c[n]);
		}
	}
}

//////////////////////////////////

int compute_coefficients(double *data, int dims[3], int degree) {
	//double *line;
	std::vector< double > line;
	double Pole[2];
	int NbPoles;
	int x,y,z;

	// recover the poles from a lookup table
	switch (degree) {
		case 2L:
			NbPoles = 1;
			Pole[0] = sqrt(8.0) - 3.0;
			break;
		case 3L:
			NbPoles = 1;
			Pole[0] = sqrt(3.0) - 2.0;
			break;
		case 4L:
			NbPoles = 2;
			Pole[0] = sqrt(664.0 - sqrt(438976.0)) + sqrt(304.0) - 19.0;
			Pole[1] = sqrt(664.0 + sqrt(438976.0)) - sqrt(304.0) - 19.0;
			break;
		case 5L:
			NbPoles = 2;
			Pole[0] = sqrt(135.0 / 2.0 - sqrt(17745.0 / 4.0)) + sqrt(105.0 / 4.0)
				- 13.0 / 2.0;
			Pole[1] = sqrt(135.0 / 2.0 + sqrt(17745.0 / 4.0)) - sqrt(105.0 / 4.0)
				- 13.0 / 2.0;
			break;
		default:
			std::cerr << "Invalid spline degree\n";
			return(1);
	}


	// convert the image samples into interpolation coefficients
	// in-place separable process, along x
	//line = (double *)malloc((size_t)(dims[0] * sizeof(double)));
	line.resize( dims[0] );
	if ((int)line.size() != dims[0]) { std::cerr << "Row allocation failed\n"; return(1); }
	for (y = 0L; y < dims[1]; y++) {
		for (z = 0L; z < dims[2]; z++) {
			get_line(data, 0, y, z, &line[0], dims);
			ConvertToInterpolationCoefficients(&line[0], dims[0], Pole, NbPoles, DBL_EPSILON);
			put_line(data, 0, y, z, &line[0], dims);
		}
	}
	//free(line);

	// in-place separable process, along y
	//line = (double *)malloc((size_t)(dims[1] * sizeof(double)));
	line.resize( dims[1] );
	if ((int)line.size() != dims[1]) { std::cerr << "Row allocation failed\n"; return(1); }
	for (x = 0L; x < dims[0]; x++) {
		for (z = 0L; z < dims[2]; z++) {
			get_line(data, 1, x, z, &line[0], dims);
			ConvertToInterpolationCoefficients(&line[0], dims[1], Pole, NbPoles, DBL_EPSILON);
			put_line(data, 1, x, z, &line[0], dims);
		}
	}
	//free(line);

	// in-place separable process, along z
	//line = (double *)malloc((size_t)(dims[2] * sizeof(double)));
	line.resize( dims[2] );
	if ((int)line.size() != dims[2]) { std::cerr << "Row allocation failed\n"; return(1); }
	for (x = 0L; x < dims[0]; x++) {
		for (y = 0L; y < dims[1]; y++) {
			get_line(data, 2, x, y, &line[0], dims);
			ConvertToInterpolationCoefficients(&line[0], dims[2], Pole, NbPoles, DBL_EPSILON);
			put_line(data, 2, x, y, &line[0], dims);
		}
	}
	//free(line);

	return(0);
}


int grad3(double grad[3], double *Bcoeff, int dims[3], double X[3], int degree) {
	double wt[3][6];
	double w, w2, w4, t, t0, t1;
	double sum_k, sum_jk;
	int idx[3][6];
	int i,j,k, pt, dim, gradDim;

	//double wtP[6], wtM[6];

	// compute interpolation indexes
	for (gradDim=0; gradDim<3; gradDim++) {
		for (dim=0; dim<3; dim++) {
			pt = (int)floor(X[dim] - (degree-1) / 2.0);
			for (i = 0L; i <= degree; i++)
				idx[dim][i] = pt++;

			// compute the interpolation weights
			if (dim == gradDim) {
				switch (degree) {
					// to do: other orders
					case 3L:
						w = X[dim] - (double)idx[dim][1];
						wt[dim][3] = (1.0 / 2.0) * w * w;
						wt[dim][0] = (w - 1.0/2.0) - wt[dim][3];
						wt[dim][2] = 1.0 + wt[dim][0] - 2.0 * wt[dim][3];
						wt[dim][1] = - wt[dim][0] - wt[dim][2] - wt[dim][3];
						break;
					default:
						std::cerr << "Invalid spline degree\n";
						return(1);
				}
			} else {
				switch (degree) {
					case 2L:
						w = X[dim] - (double)idx[dim][1];
						wt[dim][1] = 3.0 / 4.0 - w * w;
						wt[dim][2] = (1.0 / 2.0) * (w - wt[dim][1] + 1.0);
						wt[dim][0] = 1.0 - wt[dim][1] - wt[dim][2];
						break;
					case 3L:
						w = X[dim] - (double)idx[dim][1];
						wt[dim][3] = (1.0 / 6.0) * w * w * w;
						wt[dim][0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - wt[dim][3];
						wt[dim][2] = w + wt[dim][0] - 2.0 * wt[dim][3];
						wt[dim][1] = 1.0 - wt[dim][0] - wt[dim][2] - wt[dim][3];
						break;
					case 4L:
						w = X[dim] - (double)idx[dim][2];
						w2 = w * w;
						t = (1.0 / 6.0) * w2;
						wt[dim][0] = 1.0 / 2.0 - w;
						wt[dim][0] *= wt[dim][0];
						wt[dim][0] *= (1.0 / 24.0) * wt[dim][0];
						t0 = w * (t - 11.0 / 24.0);
						t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
						wt[dim][1] = t1 + t0;
						wt[dim][3] = t1 - t0;
						wt[dim][4] = wt[dim][0] + t0 + (1.0 / 2.0) * w;
						wt[dim][2] = 1.0 - wt[dim][0] - wt[dim][1] - wt[dim][3] - wt[dim][4];
						break;
					case 5L:
						w = X[dim] - (double)idx[dim][2];
						w2 = w * w;
						wt[dim][5] = (1.0 / 120.0) * w * w2 * w2;
						w2 -= w;
						w4 = w2 * w2;
						w -= 1.0 / 2.0;
						t = w2 * (w2 - 3.0);
						wt[dim][0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - wt[dim][5];
						t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
						t1 = (-1.0 / 12.0) * w * (t + 4.0);
						wt[dim][2] = t0 + t1;
						wt[dim][3] = t0 - t1;
						t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
						t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
						wt[dim][1] = t0 + t1;
						wt[dim][4] = t0 - t1;
						break;
					default:
						std::cerr << "Invalid spline degree\n";
						return(1);
				}
			}

			// _periodic_ boundary conditions
			for (i = 0L; i <= degree; i++) {
				if (dims[dim] == 1)
					idx[dim][i] = 0;
				else
					idx[dim][i] = idx[dim][i] % dims[dim];
				if (idx[dim][i] < 0)
					idx[dim][i] += dims[dim];
			}
		}

		grad[gradDim] = 0.0;
		for (i = 0; i <= degree; i++) {  // x
			sum_jk = 0.0;
			for (j = 0; j <= degree; j++) {  // y
				sum_k = 0.0;
				for (k = 0; k <= degree; k++) {  // z
					sum_k += wt[2][k] * Bcoeff[idx[0][i]*dims[1]*dims[2] + idx[1][j]*dims[2] + idx[2][k]];
				}
				sum_jk += wt[1][j] * sum_k;
			}
			grad[gradDim] += wt[0][i] * sum_jk;
		}
		// fprintf(stderr, "---------------------------\n");
	}

	return(0);
}


double interp3(double *Bcoeff, int dims[3], double X[3], int degree) {
	double wt[3][6];
	double value;
	double w, w2, w4, t, t0, t1;
	double sum_k, sum_jk;
	int idx[3][6];
	int i,j,k, pt, dim;

	// interpolation indexes
	for (dim=0; dim<3; dim++) {
		pt = (int)floor(X[dim] - (degree-1) / 2.0);
		for (i = 0L; i <= degree; i++)
			idx[dim][i] = pt++;

		// interpolation weights
		switch (degree) {
			case 2L:
				w = X[dim] - (double)idx[dim][1];
				wt[dim][1] = 3.0 / 4.0 - w * w;
				wt[dim][2] = (1.0 / 2.0) * (w - wt[dim][1] + 1.0);
				wt[dim][0] = 1.0 - wt[dim][1] - wt[dim][2];
				break;
			case 3L:
				w = X[dim] - (double)idx[dim][1];
				wt[dim][3] = (1.0 / 6.0) * w * w * w;
				wt[dim][0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - wt[dim][3];
				wt[dim][2] = w + wt[dim][0] - 2.0 * wt[dim][3];
				wt[dim][1] = 1.0 - wt[dim][0] - wt[dim][2] - wt[dim][3];
				break;
			case 4L:
				w = X[dim] - (double)idx[dim][2];
				w2 = w * w;
				t = (1.0 / 6.0) * w2;
				wt[dim][0] = 1.0 / 2.0 - w;
				wt[dim][0] *= wt[dim][0];
				wt[dim][0] *= (1.0 / 24.0) * wt[dim][0];
				t0 = w * (t - 11.0 / 24.0);
				t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
				wt[dim][1] = t1 + t0;
				wt[dim][3] = t1 - t0;
				wt[dim][4] = wt[dim][0] + t0 + (1.0 / 2.0) * w;
				wt[dim][2] = 1.0 - wt[dim][0] - wt[dim][1] - wt[dim][3] - wt[dim][4];
				break;
			case 5L:
				w = X[dim] - (double)idx[dim][2];
				w2 = w * w;
				wt[dim][5] = (1.0 / 120.0) * w * w2 * w2;
				w2 -= w;
				w4 = w2 * w2;
				w -= 1.0 / 2.0;
				t = w2 * (w2 - 3.0);
				wt[dim][0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - wt[dim][5];
				t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
				t1 = (-1.0 / 12.0) * w * (t + 4.0);
				wt[dim][2] = t0 + t1;
				wt[dim][3] = t0 - t1;
				t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
				t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
				wt[dim][1] = t0 + t1;
				wt[dim][4] = t0 - t1;
				break;
			default:
				std::cerr << "Invalid spline degree\n";
				return(0.0);
		}

		// _periodic_ boundary conditions
		for (i = 0L; i <= degree; i++) {
			if (dims[dim] == 1)
				idx[dim][i] = 0;
			else
				idx[dim][i] = idx[dim][i] % dims[dim];
			if (idx[dim][i] < 0)
				idx[dim][i] += dims[dim];
		}
		// mirror boundary conditions
		// 		dims2[dim] = 2*dims[dim]-2;
		// 		for (k = 0L; k <= degree; k++) {
		// 			idx[dim][k] = (dims[dim] == 1L) ? (0L) : ((idx[dim][k] < 0L) ?
		// 				(-idx[dim][k] - dims2[dim] * ((-idx[dim][k]) / dims2[dim]))
		// 				: (idx[dim][k] - dims2[dim] * (idx[dim][k] / dims2[dim])));
		// 			if (dims[dim] <= idx[dim][k]) {
		// 				idx[dim][k] = dims2[dim] - idx[dim][k];
		// 			}
		// 		}
	}

	value = 0.0;
	for (i = 0; i <= degree; i++) {  // x
		sum_jk = 0.0;
		for (j = 0; j <= degree; j++) {  // y
			sum_k = 0.0;
			for (k = 0; k <= degree; k++) {  // z
				sum_k += wt[2][k] * Bcoeff[idx[0][i]*dims[1]*dims[2] + idx[1][j]*dims[2] + idx[2][k]];
			}
			sum_jk += wt[1][j] * sum_k;
		}
		value += wt[0][i] * sum_jk;
	}

	return(value);
}

}
}
}
}

