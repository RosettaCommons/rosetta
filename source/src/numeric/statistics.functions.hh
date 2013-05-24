// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/statistics.functions.hh
/// @brief  a collection of various functions to compute statistics. feel free to add your own
/// @author Florian Richter (floric@u.washington.edu), sep 08


#ifndef INCLUDED_numeric_statistics_functions_hh
#define INCLUDED_numeric_statistics_functions_hh


// Platform headers
#include <platform/types.hh>

#include <numeric/types.hh>
#include <utility/vector1.hh>

// C++ headers
#include <algorithm>
#include <cmath>


namespace numeric {
namespace statistics {

/// @brief mean value of an input vector
template< class Iterator, typename T>
inline
T
mean( Iterator first, Iterator last, T )
{
	T mean = T(0);
	int size(0);
	for( ; first != last; ++first){
		mean += *first;
		size++;
		//	std::cout << "D2H " << *first << " " << size << std::endl;
	}
	return mean / size;
}


template< class Iterator, typename T>
inline
T
std_dev_with_provided_mean( Iterator first, Iterator last, T mean )
{
	T std_dev2 = T(0);
	int size(0);
	for( ; first != last; ++first){
		T meandev = *first - mean;
		std_dev2 += (meandev * meandev);
		size++;
	}
  std_dev2 /= size;
	return sqrt( std_dev2 );
}

template< class Iterator, typename T>
inline
T
std_dev( Iterator first, Iterator last, T )
{
	T meanval = mean( first, last, *first );
	return std_dev_with_provided_mean( first, last, meanval );
}

template< class T >
T errfc( T x, double tol=1e-12);

template< class T >
T errf( T x, double tol=1e-12);

///@brief error function (provided since it is not available on all platforms' math.h)
///       implemented using a Taylor series expansion
template< class T >
T errf( T x, double tol/*=1e-12*/) {
	static const double two_sqrtpi = 1.128379167095512574;        // 2/sqrt(pi)
	if (fabs(x) > 2.2) {
		return 1.0 - errfc(x);     // use continued fraction when fabs(x) > 2.2
	}
	double sum=x, term=x, xsqr=x*x;
	int j= 1;
	do {
		term*= xsqr/j;
		sum-= term/(2*j+1);
		++j;
		term*= xsqr/j;
		sum+= term/(2*j+1);
		++j;
	} while (fabs(term/sum) > tol);
	return two_sqrtpi*sum;
}

/// @brief complementary error function (provided since it is not available on
/// all platforms' math.h)
///       implemented as continued fraction
template< class T >
T errfc( T x, double tol/*=1e-12*/) {
	static const double one_sqrtpi = 0.564189583547756287;        // 1/sqrt(pi)
	if (fabs(x) < 2.2) {
		return 1.0 - errf(x);  // use series when fabs(x) < 2.2
	}
	if (x<0) {               // continued fraction only valid for x>0
		return 2.0 - errfc(-x);
	}
	double a=1,b=x;        //last two convergent numerators
	double c=x,d=x*x+0.5;  //last two convergent denominators
	double q1,q2= b/d;     //last two convergents (a/c and b/d)
	double n= 1.0, t;
	do {
		t = a*n+b*x;
		a = b;
		b = t;
		t = c*n+d*x;
		c = d;
		d = t;
		n += 0.5;
		q1 = q2;
		q2 = b/d;
	} while (fabs(q1-q2)/q2 > tol);
	return one_sqrtpi*exp(-x*x)*q2;
}

/// @brief Returns the Kullback-Leibler divergence (aka relative entropy)
/// between two discrete probability distributions.
numeric::Real kl_divergence(
	utility::vector1< numeric::Real > const & prior,
	utility::vector1< numeric::Real > const & posterior
);


numeric::Real 
corrcoef(
					utility::vector1< numeric::Real > const & vec1,
					utility::vector1< numeric::Real > const & vec2);

numeric::Real 
corrcoef_with_provided_mean_and_std_dev(
																				 utility::vector1< numeric::Real > const & vec1,
																				 numeric::Real m1,
																				 numeric::Real sd1,
																				 utility::vector1< numeric::Real > const & vec2,
																				 numeric::Real m2,
																				 numeric::Real sd2);

numeric::Real 
cov(
					utility::vector1< numeric::Real > const & vec1,
					utility::vector1< numeric::Real > const & vec2);

numeric::Real 
cov_with_provided_mean(
																				 utility::vector1< numeric::Real > const & vec1,
																				 numeric::Real m1,
																				 utility::vector1< numeric::Real > const & vec2,
																				 numeric::Real m2);


} // namespace statistics
} // namespace numeric


#endif // INCLUDED_numeric_statistics_functions_HH
