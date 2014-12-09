// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/Filter.hh
/// @brief  True/False query class which handles basic boolean logic
/// @author Will Sheffler (willsheffler@gmail.com)
/// @date   Thu Aug  10 19:49:23 2007
///

#ifndef INCLUDED_utility_query_types_HH
#define INCLUDED_utility_query_types_HH

#include <string>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>

namespace utility {
namespace query {

#define C_TN_OP utility::pointer::owning_ptr< Converter<T,N> >
#define C1P_TAN_OP utility::pointer::owning_ptr< Converter1Param<T,A,N> >
#define C1P_TAB_OP utility::pointer::owning_ptr< Converter1Param<T,A,bool> >
#define C1P_TAR_OP utility::pointer::owning_ptr< Converter1Param<T,A,Real> >

template<class T,class N>
struct Converter : public utility::pointer::ReferenceCount
{
	// typedef utility::pointer::owning_ptr<Converter<T,N> > OP;
	virtual N convert( T arg ) = 0;
	virtual std::string description() = 0;
	// virtual C_TN_OP clone() = 0;
};

template<class T, class A, class N>
struct Converter1Param : public utility::pointer::ReferenceCount
{
	// typedef utility::pointer::owning_ptr<Converter<T,N> > OP;
	virtual N convert( T arg, A arg2 ) = 0;
	virtual std::string description() = 0;
	// virtual C_TN_OP clone() = 0;
};


template<class T,class N>
struct ImplicitConverter : public Converter<T,N> {
	// typedef utility::pointer::owning_ptr<Converter<T,N> > OP;
	inline N convert( T arg ) { return arg; }
	std::string description() { return "(implicit conversion)"; }
	// C_TN_OP clone() { return new ImplicitConverter; }
};




} // end namespace query
} // end namespace utility




#endif
