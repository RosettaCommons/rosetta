// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/statistics/functions.hh
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
#include <complex>


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
	for ( ; first != last; ++first ) {
		mean += *first;
		size++;
		// std::cout << "D2H " << *first << " " << size << std::endl;
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
	for ( ; first != last; ++first ) {
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

//FPD erf and erfc with real & imaginary arguments
//    code borrowed from http://ab-initio.mit.edu/Faddeeva

// compute w(z) = exp(-z^2) erfc(-iz) [ Faddeeva / scaled complex error func ]
std::complex<double> w(std::complex<double> z,double relerr=0);
double w_im(double x); // special-case code for Im[w(x)] of real x

// compute erfcx(z) = exp(z^2) erfc(z)
std::complex<double> errfcx(std::complex<double> z, double relerr=0);
double errfcx(double x);

// compute erf(z), the error function of complex arguments
std::complex<double> errf(std::complex<double> z, double relerr=0);
double errf(double x);

// compute erfi(z) = -i erf(iz), the imaginary error function
std::complex<double> errfi(std::complex<double> z, double relerr=0);
double errfi(double x);

// compute erfc(z) = 1 - erf(z), the complementary error function
std::complex<double> errfc(std::complex<double> z, double relerr=0);
double errfc(double x);

// compute Dawson(z) = sqrt(pi)/2  *  exp(-z^2) * erfi(z)
std::complex<double> Dawson(std::complex<double> z, double relerr=0);
double Dawson(double x); // special case for real x

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
