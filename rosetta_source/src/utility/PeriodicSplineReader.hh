// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/PeriodicSplineReader.hh
/// @brief  Class for reading in Splines and Calculating values of a function, first order derivatives
/// @author Kristian Kaufmann


#ifndef INCLUDED_utility_PeriodicSplineReader_hh
#define INCLUDED_utility_PeriodicSplineReader_hh


// Unit headers
#include <utility/PeriodicSplineReader.fwd.hh>

// Utility headers
#include <utility/io/izstream.hh>

// C++ headers
#include <iostream>
#include <vector>


namespace utility {


/// @brief Class: PeriodicSplineReader
class PeriodicSplineReader
{


public: // Creation


	/// @brief Default constructor
	inline
	PeriodicSplineReader() :
		start_( 0.0 ),
		delta_( 0.0 )
	{}


	/// @brief Copy constructor
	inline
	PeriodicSplineReader( PeriodicSplineReader const & spline ) :
		start_( spline.start_ ),
		delta_( spline.delta_ ),
		values_( spline.values_ ),
		dsecox_( spline.dsecox_ )
	{}


	/// @brief Destructor
	inline
	//virtual // Don't need this if no class hierarchy
	~PeriodicSplineReader()
	{}


public: // Inspectors


	/// @brief return value at certain x
	inline
	double
	F( double const x )
	{
		int dim = int(values_.size());
		//determine i with m_Start+(i-1)*m_Delta < x < m_Start+i*m_Delta
		// for the correct supporting points
		int    i(0);
		while  (start_+i*delta_<x) ++i;
		if(!i){
			while  (start_+i*delta_>x) --i;
			++i;
		}

		//see Numerical recipes in C++, pages 116-118
		double delta_akt=start_+i*delta_-x;
		// the following while moves the point into range making
		// the assumption that the spline is periodic
		// the modulus by dim takes care of the rest
		while (i<0) i+=dim;
		double y0(values_.at(size_t((i+dim-1)%dim)));
		double y1(values_.at(size_t(i%dim)));
		double ddy0(dsecox_.at(size_t((i+dim-1)%dim)));
		double ddy1(dsecox_.at(size_t(i%dim)));
		double a(delta_akt/delta_);
		double b( 1 - a);
		double c(( a*a*a - a) * (delta_*delta_) / 6);
		double d(( b*b*b - b) * ( delta_*delta_) / 6);
		double y( a * y0 + b * y1 + c * ddy0 + d * ddy1);
		return y;
	}


	/// @brief return derivative at certain x
	inline
	double dF( double const x )
	{
		int dim = int(values_.size());
		//determine i with start_+(i-1)*delta_ < x < start_+i*delta_
		// for the correct supporting points
		int    i(0);
		while ( start_+i*delta_<x) ++i;
		if(!i){
			while ( start_+i*delta_>x) --i;
			++i;
		}
		//see Numerical recipes in C++, pages 116-118
		double delta_akt=start_+i*delta_-x;
		// the following while moves the point into range making
		// the assumption that the spline is periodic
		// the modulus by dim takes care of the rest.
		while (i<0) i+=dim;
		double y0(values_.at(size_t((i+dim-1)%dim)));
		double y1(values_.at(size_t(i%dim)));
		double ddy0(dsecox_.at(size_t((i+dim-1)%dim)));
		double ddy1(dsecox_.at(size_t(i%dim)));
		double a(delta_akt/delta_);
		double b( 1 - a);

		return ( y1-y0 )/delta_ - ( 3*a*a-1 )/6*delta_*ddy0
		 + ( 3*b*b-1 )/6 * delta_ * ddy1;
	}


public: // Modifiers


	/// @brief Read from a stream
	inline
	int
	read( utility::io::izstream &iunit )
	{
		// read parameters
		int dimension;
		double in_data;
		iunit >> start_; ///start value of the spline
		iunit >> delta_; /// distance between points
		iunit >> dimension; /// number of points
		int i=1;
		// following while block reads in the values of the function
		while (!(iunit.fail()) && i<dimension+1) {
			iunit >> in_data;
			if( !(iunit.fail()) && !(iunit.eof())) {
				values_.push_back(in_data);
				++i;
			}
		}
		i=1;
		// following block reads in the second order derivatives of the function
		while (!(iunit.fail()) && i<dimension+1) {
			iunit >> in_data;
			if( !(iunit.fail()) && !(iunit.eof())) {
				dsecox_.push_back(in_data);
				++i;
			}
		}
		if (iunit.fail()) {
			std::cout << "read spline failed." << std::endl;
			return -1;
		}
		iunit.clear();
		return 1;
	}


private: // Fields


	/// @brief Gives the arguments as a sequence of equidistant points
	double start_, delta_;

	/// @brief f(x)
	std::vector< double > values_;

	/// @brief Second order derivatives
	std::vector< double > dsecox_;


}; // Class PeriodicSplineReader


} // namespace utility


#endif // INCLUDED_utility_PeriodicSplineReader_HH
