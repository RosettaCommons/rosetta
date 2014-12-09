// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/numeric/interpolation/spline_functions.cc
/// @brief  Interpolation with cubic splines
/// @author Will Sheffler
///

#include <numeric/interpolation/spline/spline_functions.hh>

#include <iostream>

namespace numeric {
namespace interpolation {
namespace spline {



//Given arrays x[1..n] and y[1..n] containing a tabulated function,
//i.e., y i = f(xi),with x1 < x2 <...<xN , and given values yp1 and
//ypn for the first derivative of the interpolating function at
//points 1 and n, respectively, this routine returns an array y2[1..n]
//that contains the second derivatives of the interpolating function
//at the tabulated points xi.Ifyp1 and/or ypn are equal to 1  10 30
//or larger, the routine is signaled to set the corresponding
//boundary condition for a natural spline, with zero second
//derivative on that boundary.
utility::vector1<Real>
spline_second_derivative(
 utility::vector1<Real> const & x,
 utility::vector1<Real> const & y,
 Real yp1,
 Real ypn
)
{

	Real p,qn,sig,un;

	// using namespace std;
	// cerr << "spline second derivative 1 yp1 " << yp1 << " ypn " << ypn<< endl;
	// cerr << "                      x:" ;
	// for(int ii=1; ii<=(int)x.size();ii++)
	// 	cerr << ' ' << x[ii] ;
	// cerr << endl;
	// cerr << "                      y:" ;
	// for(int ii=1; ii<=(int)y.size();ii++)
	// 	cerr << ' ' << y[ii] ;
	// cerr << endl;
	// std::cerr << "starting spline_second_derivative" << std::endl;

	int n = x.size();
	utility::vector1<Real> y2(n);
	utility::vector1<Real> u(n-1);

	// cout << "spline second derivative 2 "<< n <<  endl;
	//The lower boundary condition is set either to be  nat-ural  y2[1]=u[1]=0.0;
	if (yp1 > 0.99e30) {
		//or else to have a specified first derivative.
		y2[1] = u[1]=0.0;
	} else {
		y2[1] = -0.5;
		u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
	}
	// cout << "spline second derivative 3 " << u[1] << endl;
	//This is the decomposition loop of the tridiagonal al-gorithm.
	//y2 and u are used for tem-porary storage of the decomposed factors.

	for (int i=2; i<=n-1; i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
		// cout << i << ' ' << y2[i] << ' ' << u[i] << endl;
	}
	// cout << "spline second derivative 4 "<< endl;
	//The upper boundary condition is set either to be  natural  qn=un=0.0;
	if (ypn > 0.99e30) {
		qn = un = 0.0;
	} else { //or else to have a specified first derivative.
		qn=0.5;
		un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
	}
	// cout << "spline second derivative 5 "<< endl;
	//This is the backsubstitution loop of the tridiagonal algorithm
	y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
	// cout << "spline second derivative 6 "<< endl;
	for( int k=n-1; k>=1; k--) {
		// cout << k << ' ' << y2[k] << ' ' << u[k] << endl;
		y2[k]=y2[k]*y2[k+1]+u[k];
	}
	// cout << "spline second derivative 7 "<< endl;
	return y2;
}

// !!!!!!!!!! THIS COULD BE MADE MUCH FASTER IF XA ARE IN ORDER
//                  AND IT WAS MODIFIED NOT TO LOOK UP "BIN"
//Given the arrays xa[1..n] and ya[1..n], which tabulate a function
//(with the xai s in order), and given the array y2a[1..n], which
//is the output from spline above, and given a value of x, this
//routine returns a cubic-spline interpolated value y.
void spline_interpolate(
 utility::vector1<Real> const & xa,
 utility::vector1<Real> const & ya,
 utility::vector1<Real> const & y2a,
 Real x, Real & y, Real & dy
)
{
	int klo,khi,k;
	Real h,b,a;

	//We will find the right place in the table by means of bisection.
	//This is optimal if sequential calls to this routine are at
	//random values of x. If sequential calls are in order, and
	//closely spaced, one would do better to store previous values
	//of klo and khi and test if they remain appropriate on the next call.
	int n = xa.size();
	klo=1;
	khi=n;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x)
			khi=k;
		else
			klo=k;
	}
	//klo and khi now bracket the input value of x.
	h=xa[khi]-xa[klo];
	// if (h == 0.0)
		// cout << "spline interpolate: Bad xa input to routine splint" << endl;
	//The xa s must be dis-tinct.
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	//Cubic spline polynomial is now evaluated.
	y  = a*ya[klo]+b*ya[khi] + ( (a*a*a-a)*y2a[klo] + (b*b*b-b)*y2a[khi] ) *(h*h)/6.0;
	dy = (ya[khi]-ya[klo])/h + ( (3*b*b-1)*y2a[khi] - (3*a*a-1)*y2a[klo] ) *  h  /6.0;
}

} // end namespace spline
} // end namespace interpolation
} // end namespace numeric
