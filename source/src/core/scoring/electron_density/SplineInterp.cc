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

void put_line3(double* data, int dim, int x1, int x2, double line[], int dims[]) {
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


void get_line3(double* data, int dim, int x1, int x2, double line[], int dims[]) {
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

void put_line4(double* data, int dim, int x1, int x2, int x3, double line[], int dims[]) {
	int i, inc;
	double* ptr;

	if (dim == 0) {          // x1 == y, x2 == z, x3 == w
		ptr = &data[x1*dims[2]*dims[3] + x2*dims[3] + x3];
		inc = dims[1]*dims[2]*dims[3];
	} else if (dim == 1) {   // x1 == x, x2 == z, x3 == w
		ptr = &data[x1*dims[1]*dims[2]*dims[3] + x2*dims[3]+ x3];
		inc = dims[2]*dims[3];
	} else if (dim == 2) {   // x1 == x, x2 == y, x3 == w
		ptr = &data[x1*dims[1]*dims[2]*dims[3] + x2*dims[2]*dims[3]+ x3];
		inc = dims[3];
	} else {                 // x1 == x, x2 == y, x3 == z
		ptr = &data[x1*dims[1]*dims[2]*dims[3] + x2*dims[2]*dims[3] + x3*dims[3]];
		inc = 1;
	}

	for (i=0; i<dims[dim]; i++) {
		*ptr = line[i] ;
		ptr = &ptr[inc];
	}
}


void get_line4(double* data, int dim, int x1, int x2, int x3, double line[], int dims[]) {
	int i, inc;
	double* ptr;

	if (dim == 0) {          // x1 == y, x2 == z, x3 == w
		ptr = &data[x1*dims[2]*dims[3] + x2*dims[3] + x3];
		inc = dims[1]*dims[2]*dims[3];
	} else if (dim == 1) {   // x1 == x, x2 == z, x3 == w
		ptr = &data[x1*dims[1]*dims[2]*dims[3] + x2*dims[3]+ x3];
		inc = dims[2]*dims[3];
	} else if (dim == 2) {   // x1 == x, x2 == y, x3 == w
		ptr = &data[x1*dims[1]*dims[2]*dims[3] + x2*dims[2]*dims[3]+ x3];
		inc = dims[3];
	} else {                 // x1 == x, x2 == y, x3 == z
		ptr = &data[x1*dims[1]*dims[2]*dims[3] + x2*dims[2]*dims[3] + x3*dims[3]];
		inc = 1;
	}

	for (i=0; i<dims[dim]; i++) {
		line[i] = *ptr;
		ptr = &ptr[inc];
	}
}

static double	InitialCausalCoefficient (
		double c[],       // coefficients
		long DataLength,  // number of coefficients
		double z,         // actual pole
		double Tolerance, // admissible relative error
		bool Mirrored     // mirror boundary?
) {
	double	Sum, zn, iz, z2n;
	long	n, Horizon;

	// this initialization corresponds to mirror boundaries
	// modified FPD -- periodic boundaries
	Horizon = DataLength;
	if (Tolerance > 0.0) {
		Horizon = (long)ceil(log(Tolerance) / log(fabs(z)));
	}

	if (!Mirrored) {
		if (Horizon < DataLength) {
			// accelerated loop
			zn = z;
			Sum = c[0];
			for (n = 1L; n < Horizon; n++) {
				Sum += zn * c[DataLength - n];
				zn *= z;
			}
			return(Sum);
		} else {
			// full loop
			zn = z;
			Sum = c[0];
			for (n = 1L; n < DataLength; n++) {
				Sum += zn * c[DataLength - n];
				zn *= z;
			}
			return(Sum / (1.0 - zn));
		}
	} else {
		if (Horizon < DataLength) {
			// accelerated loop
			zn = z;
			Sum = c[0];
			for (n = 1L; n < Horizon; n++) {
				Sum += zn * c[n];
				zn *= z;
			}
			return(Sum);
		}
		else {
			// full loop
			zn = z;
			iz = 1.0 / z;
			z2n = std::pow(z, (double)(DataLength - 1));
			Sum = c[0] + z2n * c[DataLength - 1];
			z2n *= z2n * iz;
			for (n = 1L; n <= DataLength - 2; n++) {
				Sum += (zn + z2n) * c[n];
				zn *= z;
				z2n *= iz;
			}
			return(Sum / (1.0 - zn * zn));
		}
	}
}


static double	InitialAntiCausalCoefficient (
		double c[],       // coefficients
		long DataLength,  // number of coefficients
		double z,         // actual pole
		double Tolerance, // admissible relative error
		bool Mirrored     // mirror boundary?
) {
	double Sum, zn;
	int n, Horizon;

	Horizon = DataLength;
	if (Tolerance > 0.0) {
		Horizon = (long)ceil(log(Tolerance) / log(fabs(z)));
	}

	if (!Mirrored) {
		if (Horizon < DataLength) {
			zn = z;
			Sum = c[DataLength-1];
			for (n = 0L; n < Horizon; n++) {
				Sum += zn * c[n];
				zn *= z;
			}
			return(-z*Sum);
		}
		else {
			zn = z;
			Sum = c[DataLength-1];
			for (n = 0L; n < DataLength-1; n++) {
				Sum += zn * c[n];
				zn *= z;
			}
			return(-z*Sum / (1.0 - zn));
		}
	} else {
		return((z / (z * z - 1.0)) * (z * c[DataLength - 2] + c[DataLength - 1]));
	}
}



void ConvertToInterpolationCoefficients (
		double c[],		    // input samples --> output coefficients
		long   DataLength,// number of samples or coefficients
		double z[],       // poles
		long   NbPoles,   // number of poles
		double Tolerance, // admissible relative error
		bool   Mirrored   // mirror boundary?
) {
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
		c[0] = InitialCausalCoefficient(c, DataLength, z[k], Tolerance, Mirrored);
		// causal recursion
		for (n = 1L; n < DataLength; n++) {
			c[n] += z[k] * c[n - 1L];
		}
		// anticausal initialization
		c[DataLength - 1L] = InitialAntiCausalCoefficient(c, DataLength, z[k], Tolerance, Mirrored);
		// anticausal recursion
		for (n = DataLength - 2L; 0 <= n; n--) {
			c[n] = z[k] * (c[n + 1L] - c[n]);
		}
	}
}

//////////////////////////////////

int compute_coefficients3(double *data, int dims[3]) {
	std::vector< double > line;
	double Pole[2];
	int NbPoles;
	int x,y,z;

	NbPoles = 1;
	Pole[0] = sqrt(3.0) - 2.0;

	// convert the image samples into interpolation coefficients
	// in-place separable process, along x
	//line = (double *)malloc((size_t)(dims[0] * sizeof(double)));
	line.resize( dims[0] );
	if ((int)line.size() != dims[0]) { std::cerr << "Row allocation failed\n"; return(1); }
	for (y = 0L; y < dims[1]; y++) {
		for (z = 0L; z < dims[2]; z++) {
			get_line3(data, 0, y, z, &line[0], dims);
			ConvertToInterpolationCoefficients(&line[0], dims[0], Pole, NbPoles, DBL_EPSILON, false);
			put_line3(data, 0, y, z, &line[0], dims);
		}
	}

	// in-place separable process, along y
	//line = (double *)malloc((size_t)(dims[1] * sizeof(double)));
	line.resize( dims[1] );
	if ((int)line.size() != dims[1]) { std::cerr << "Row allocation failed\n"; return(1); }
	for (x = 0L; x < dims[0]; x++) {
		for (z = 0L; z < dims[2]; z++) {
			get_line3(data, 1, x, z, &line[0], dims);
			ConvertToInterpolationCoefficients(&line[0], dims[1], Pole, NbPoles, DBL_EPSILON, false);
			put_line3(data, 1, x, z, &line[0], dims);
		}
	}

	// in-place separable process, along z
	//line = (double *)malloc((size_t)(dims[2] * sizeof(double)));
	line.resize( dims[2] );
	if ((int)line.size() != dims[2]) { std::cerr << "Row allocation failed\n"; return(1); }
	for (x = 0L; x < dims[0]; x++) {
		for (y = 0L; y < dims[1]; y++) {
			get_line3(data, 2, x, y, &line[0], dims);
			ConvertToInterpolationCoefficients(&line[0], dims[2], Pole, NbPoles, DBL_EPSILON, false);
			put_line3(data, 2, x, y, &line[0], dims);
		}
	}

	return(0);
}


int grad3(double grad[3], double *Bcoeff, int dims[3], double X[3]) {
	double wt[3][4];
	double w, w2, w4, t, t0, t1;
	double sum_k, sum_jk;
	int idx[3][4];
	int i,j,k, pt, dim, gradDim;

	// compute interpolation indexes
	for (gradDim=0; gradDim<3; gradDim++) {
		for (dim=0; dim<3; dim++) {
			pt = (int)floor(X[dim] - (3-1) / 2.0);
			for (i = 0L; i <= 3; i++)
				idx[dim][i] = pt++;

			// compute the interpolation weights
			if (dim == gradDim) {
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

			// _periodic_ boundary conditions
			for (i = 0L; i <= 3; i++) {
				if (dims[dim] == 1)
					idx[dim][i] = 0;
				else
					idx[dim][i] = idx[dim][i] % dims[dim];
				if (idx[dim][i] < 0)
					idx[dim][i] += dims[dim];
			}
		}

		grad[gradDim] = 0.0;
		for (i = 0; i <= 3; i++) {  // x
			sum_jk = 0.0;
			for (j = 0; j <= 3; j++) {  // y
				sum_k = 0.0;
				for (k = 0; k <= 3; k++) {  // z
					sum_k += wt[2][k] * Bcoeff[idx[0][i]*dims[1]*dims[2] + idx[1][j]*dims[2] + idx[2][k]];
				}
				sum_jk += wt[1][j] * sum_k;
			}
			grad[gradDim] += wt[0][i] * sum_jk;
		}
	}

	return(0);
}


double interp3(double *Bcoeff, int dims[3], double X[3]) {
	double wt[3][4];
	double value;
	double w, w2, w4, t, t0, t1;
	double sum_k, sum_jk;
	int idx[3][4];
	int i,j,k, pt, dim;

	// interpolation indexes
	for (dim=0; dim<3; dim++) {
		pt = (int)floor(X[dim] - (3-1) / 2.0);
		for (i = 0L; i <= 3; i++)
			idx[dim][i] = pt++;

		// interpolation weights
		w = X[dim] - (double)idx[dim][1];
		wt[dim][3] = (1.0 / 6.0) * w * w * w;
		wt[dim][0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - wt[dim][3];
		wt[dim][2] = w + wt[dim][0] - 2.0 * wt[dim][3];
		wt[dim][1] = 1.0 - wt[dim][0] - wt[dim][2] - wt[dim][3];

		// _periodic_ boundary conditions
		for (i = 0L; i <= 3; i++) {
			if (dims[dim] == 1)
				idx[dim][i] = 0;
			else
				idx[dim][i] = idx[dim][i] % dims[dim];
			if (idx[dim][i] < 0)
				idx[dim][i] += dims[dim];
		}
	}

	value = 0.0;
	for (i = 0; i <= 3; i++) {  // x
		sum_jk = 0.0;
		for (j = 0; j <= 3; j++) {  // y
			sum_k = 0.0;
			for (k = 0; k <= 3; k++) {  // z
				sum_k += wt[2][k] * Bcoeff[idx[0][i]*dims[1]*dims[2] + idx[1][j]*dims[2] + idx[2][k]];
			}
			sum_jk += wt[1][j] * sum_k;
		}
		value += wt[0][i] * sum_jk;
	}

	return(value);
}


int compute_coefficients4(double *data, int dims[4]) {
	std::vector< double > line;
	double Pole[2];
	int NbPoles;
	int x,y,z,w;

	NbPoles = 1;
	Pole[0] = sqrt(3.0) - 2.0;

	// convert the image samples into interpolation coefficients
	// in-place separable process, along x
	line.resize( dims[0] );
	if ((int)line.size() != dims[0]) { std::cerr << "Row allocation failed\n"; return(1); }
	for (y = 0L; y < dims[1]; y++) {
		for (z = 0L; z < dims[2]; z++) {
			for (w = 0L; w < dims[3]; w++) {
				get_line4(data, 0, y, z, w, &line[0], dims);
				ConvertToInterpolationCoefficients(&line[0], dims[0], Pole, NbPoles, DBL_EPSILON, true);
				put_line4(data, 0, y, z, w, &line[0], dims);
			}
		}
	}

	// in-place separable process, along y
	line.resize( dims[1] );
	if ((int)line.size() != dims[1]) { std::cerr << "Row allocation failed\n"; return(1); }
	for (x = 0L; x < dims[0]; x++) {
		for (z = 0L; z < dims[2]; z++) {
			for (w = 0L; w < dims[3]; w++) {
				get_line4(data, 1, x, z, w, &line[0], dims);
				ConvertToInterpolationCoefficients(&line[0], dims[1], Pole, NbPoles, DBL_EPSILON, false);
				put_line4(data, 1, x, z, w, &line[0], dims);
			}
		}
	}

	// in-place separable process, along z
	line.resize( dims[2] );
	if ((int)line.size() != dims[2]) { std::cerr << "Row allocation failed\n"; return(1); }
	for (x = 0L; x < dims[0]; x++) {
		for (y = 0L; y < dims[1]; y++) {
			for (w = 0L; w < dims[3]; w++) {
				get_line4(data, 2, x, y, w, &line[0], dims);
				ConvertToInterpolationCoefficients(&line[0], dims[2], Pole, NbPoles, DBL_EPSILON, false);
				put_line4(data, 2, x, y, w, &line[0], dims);
			}
		}
	}

	// in-place separable process, along w
	line.resize( dims[3] );
	if ((int)line.size() != dims[3]) { std::cerr << "Row allocation failed\n"; return(1); }
	for (x = 0L; x < dims[0]; x++) {
		for (y = 0L; y < dims[1]; y++) {
			for (z = 0L; z < dims[2]; z++) {
				get_line4(data, 3, x, y, z, &line[0], dims);
				ConvertToInterpolationCoefficients(&line[0], dims[3], Pole, NbPoles, DBL_EPSILON, false);
				put_line4(data, 3, x, y, z, &line[0], dims);
			}
		}
	}

	return(0);
}


int grad4(double grad[4], double *Bcoeff, int dims[4], double X[4]) {
	double wt[4][4];
	double w, w2, w4, t, t0, t1;
	double sum_l, sum_kl, sum_jkl;
	int idx[4][4];
	int i,j,k,l, pt, dim, gradDim;

	// compute interpolation indexes
	for (gradDim=0; gradDim<4; gradDim++) {
		for (dim=0; dim<4; dim++) {
			pt = (int)floor(X[dim] - (3-1) / 2.0);
			for (i = 0L; i <= 3; i++)
				idx[dim][i] = pt++;

			// compute the interpolation weights
			if (dim == gradDim) {
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

			if (dim > 0) {
				// _periodic_ boundary conditions
				for (i = 0L; i <= 3; i++) {
					if (dims[dim] == 1)
						idx[dim][i] = 0;
					else
						idx[dim][i] = idx[dim][i] % dims[dim];
					if (idx[dim][i] < 0)
						idx[dim][i] += dims[dim];
				}
			} else {
				// flat boundary
				for (i = 0L; i <= 3; i++) {
					if (dims[dim] == 1)
						idx[dim][i] = 0;
					if (idx[dim][i] < 0)
						idx[dim][i] = 0;
					if (idx[dim][i] > dims[dim]-1)
						idx[dim][i] = dims[dim]-1;
				}
			}
		}

		grad[gradDim] = 0.0;
		for (i = 0; i <= 3; i++) {  // x
			sum_jkl = 0.0;
			for (j = 0; j <= 3; j++) {  // y
				sum_kl = 0.0;
				for (k = 0; k <= 3; k++) {  // z
					sum_l = 0;
					for (l = 0; l <= 3; l++) {  // w
						sum_l += wt[3][l] * Bcoeff[idx[0][i]*dims[1]*dims[2]*dims[3] + idx[1][j]*dims[2]*dims[3] + idx[2][k]*dims[3] + idx[3][l]];
					}
					sum_kl += wt[2][k] * sum_l;
				}
				sum_jkl += wt[1][j] * sum_kl;
			}
			grad[gradDim] += wt[0][i] * sum_jkl;
		}

	}

	return(0);
}


double interp4(double *Bcoeff, int dims[4], double X[4]) {
	double wt[4][4];
	double value;
	double w, w2, w4, t, t0, t1;
	double sum_l, sum_kl, sum_jkl;
	int idx[4][4];
	int i,j,k,l, pt, dim;

	// interpolation indexes
	for (dim=0; dim<4; dim++) {
		pt = (int)floor(X[dim] - (3-1) / 2.0);
		for (i = 0L; i <= 3; i++)
			idx[dim][i] = pt++;

		// interpolation weights
		w = X[dim] - (double)idx[dim][1];
		wt[dim][3] = (1.0 / 6.0) * w * w * w;
		wt[dim][0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - wt[dim][3];
		wt[dim][2] = w + wt[dim][0] - 2.0 * wt[dim][3];
		wt[dim][1] = 1.0 - wt[dim][0] - wt[dim][2] - wt[dim][3];

		if (dim > 0) {
			// _periodic_ boundary conditions
			for (i = 0L; i <= 3; i++) {
				if (dims[dim] == 1)
					idx[dim][i] = 0;
				else
					idx[dim][i] = idx[dim][i] % dims[dim];
				if (idx[dim][i] < 0)
					idx[dim][i] += dims[dim];
			}
		} else {
			for (i = 0L; i <= 3; i++) {
				if (idx[dim][i] < 0)
					idx[dim][i] = 0;
				if (idx[dim][i] > dims[dim]-1)
					idx[dim][i] = dims[dim]-1;
			}
		}
	}

	value = 0.0;
	for (i = 0; i <= 3; i++) {  // x
		sum_jkl = 0.0;
		for (j = 0; j <= 3; j++) {  // y
			sum_kl = 0.0;
			for (k = 0; k <= 3; k++) {  // z
				sum_l = 0;
				for (l = 0; l <= 3; l++) {  // w
					sum_l += wt[3][l] * Bcoeff[idx[0][i]*dims[1]*dims[2]*dims[3] + idx[1][j]*dims[2]*dims[3] + idx[2][k]*dims[3] + idx[3][l]];
				}
				sum_kl += wt[2][k] * sum_l;
			}
			sum_jkl += wt[1][j] * sum_kl;
		}
		value += wt[0][i] * sum_jkl;
	}

	return(value);
}

}
}
}
}

