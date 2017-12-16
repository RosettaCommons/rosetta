// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/kalngyk/cubic.cc
/// @author YK Ng & SC Li (kalngyk@gmail.com)

#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <numeric/NumericTraits.hh>
#include <core/types.hh>
#include <protocols/cluster/calibur/cubic.hh>

namespace protocols {
namespace cluster {
namespace calibur {

#define PI 3.14159265358979323846

//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -


void cubic_roots1(double a2, double a1, double a0, double * z) {
	double p_over_3 = (3*a1 - a2*a2)/9.;
	double q_over_2 = (9*a2*a1 - 27*a0 - 2*a2*a2*a2) /54.;
	if ( p_over_3 >= 0 ) {
		// can proceed to obtain root using sinh(C)
		z[0] = z[1] = z[2] = 0;
		return;
	} else { // p is negative
		p_over_3 = -p_over_3;
		double C = q_over_2 / sqrt(p_over_3 * p_over_3 * p_over_3);
		if ( C >= -1 && C <= 1 ) {
			// Then, cos(3*theta) = C. This gives 3 solutions, respectively
			// cos(acos(C)/3), cos((acos(C)+2pi)/3), cos((acos(C)+4pi)/3)
			double theta = acos(C);
			double a2_over_3 = a2/3.;
			z[0] = z[1] = z[2] = 2 * sqrt(p_over_3);
			z[0] = z[0] * cos( theta    / 3.) - a2_over_3;
			z[1] = z[1] * cos((theta + 2.*PI) / 3.) - a2_over_3;
			z[2] = z[2] * cos((theta + 4.*PI) / 3.) - a2_over_3;
		} else { // only one root is real
			// can proceed to obtain root using cosh(C)
			z[0] = z[1] = z[2] = 0;
		}
	}
}

/**
* This is faster but leaves no option for obtaining imaginary roots.
* (After "-O" in GCC, the speed gain is not much, only ~4%)
*/
void cubic_roots2(double a2, double a1, double a0, double* z) {
	double Q = (3*a1 - a2*a2)/9.;
	double R = (9*a2*a1 - 27*a0 - 2*a2*a2*a2) /54.;
	double Q3 = Q*Q*Q;
	if ( Q3 + R*R < 0 ) { // condition for 3 real roots
		// Q3+R2<0 implies Q3<0, and sqrt(-Q3) is real
		double negsqrtQ3 = sqrt(-Q3);
		double theta = acos(R/negsqrtQ3);
		double a2_over_3 = a2/3.;
		z[0] = z[1] = z[2] = 2 * sqrt(-Q);
		z[0] = z[0] * cos( theta    /3.) - a2_over_3;
		z[1] = z[1] * cos((theta + 2.*PI) /3.) - a2_over_3;
		z[2] = z[2] * cos((theta + 4.*PI) /3.) - a2_over_3;
	} else { // only one root is real
		z[0] = z[1] = z[2] = 0; // give up...
	}
}


}
}
}
