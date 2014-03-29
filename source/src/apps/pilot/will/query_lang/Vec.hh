// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   demo/utility/Vec.hh
/// @brief  a few vector<Real> helpers for use in true_false.cc
/// @author Will Sheffler (willsheffler@gmail.com)
/// @date   Thu Aug  9 19:49:23 2007
///

#ifndef demo_value_utility_Vec_HH
#define demo_value_utility_Vec_HH

#include <vector>
#include <string>
#include <iostream>

#include <core/types.hh>
#include <utility/query/types.hh>

using namespace utility::query;
using core::Real;

typedef std::vector<Real> Vec;

struct VecContains : public Converter1Param<Vec,int,bool>{
	inline bool convert( Vec v, int elem) {
		for(Vec::iterator i = v.begin(); i != v.end(); ++i)
			if( *i == elem ) return true;
		return false;
	}
	std::string description() { return "Vec contains"; }
};

struct VecMaxElement : public Converter<Vec,Real> {
	inline Real convert( Vec v ) {
		if( v.size() == 0 ) return -1234;
		Real m = v[0];
		for( Size i=0; i < v.size(); ++i )
			if( v[i] > m ) m = v[i];
		return m;
	}
	std::string description() { return "max element of vector"; }
	// NC_Vec_Real::OP clone() { return new VecMaxElement; }
};

struct VecMinElement : public Converter<Vec,Real> {
	inline Real convert( Vec v ) {
		if( v.size() == 0 ) return -1234;
		Real m = v[0];
		for( Size i=0; i < v.size(); ++i )
			if( v[i] < m ) m = v[i];
		return m;
	}
	std::string description() { return "min element of vector"; }
	// NC_Vec_Real::OP clone() { return new VecMinElement; }
};

struct VecLen : public Converter<Vec,Real> {
	inline Real convert( Vec v ) { return v.size(); }
	std::string description() { return "lenght of vector"; }
	// NC_Vec_Real::OP clone() { return new VecLen; }
};

struct VecMean : public Converter<Vec,Real> {
	inline Real convert(Vec v ) {
		float m = 0;
		for( Size i=0; i < v.size(); ++i )
			m += v[i];
		return m / (Real)v.size();
	}
	std::string description() { return "mean of vector elements"; }
	// NC_Vec_float::OP clone() { return new VecMean; }
};

std::ostream &
operator<<( std::ostream & os, Vec const & v ) {
	VecMean vecmean;
	os << "Vec( ";
	for(Vec::const_iterator i = v.begin(); i != v.end(); ++i )
		os << *i << ", ";
	float m = vecmean.convert(v);
	os << ") mean " << m;
	return os;
}

Vec rand_Vec(Real len=-1, int MAX = 10) {
	Vec v;
	if( len < 0 ) len = rand()%9+1;
	for(Real i=0; i<len; ++i) {
		v.push_back( rand()%MAX );
	}
	return v;
}


#endif
