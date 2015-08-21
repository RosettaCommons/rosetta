// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/unti/OStream.hh
/// @brief  base ostream like class for mini
/// @author Sergey Lyskov

#ifndef INCLUDED_utility_stream_util_hh
#define INCLUDED_utility_stream_util_hh

#include <map>
#include <list>
#include <ostream>

//Auto Headers
#include <utility/vector1.fwd.hh>
#include <vector>


namespace utility {

/// -------------------------------------------------------------------
/// Predefined functions for Tracer IO
/// We originally moved them out of name spaces so we can use them right away - without specifying it.
/// Returned to utility namespace, as some compilers did not find the operator otherwise! -- rhiju

template <typename T1, typename T2>
std::ostream & operator <<(std::ostream & os, std::pair<T1, T2> const & v); // declare here to allow output of vector1< std::pair >.

/// @brief Output function for std::vector object.
template <class T>
std::ostream & operator <<( std::ostream & os, std::vector<T> const & v)
{
	os << "[";
	for ( size_t i = 0; i < v.size(); i++ ) {
		os << v[i];
		if ( i < v.size()-1 ) os << ", ";
	}
	os << "]";
	return os;
}


/// @brief Output function for utility::vector1 object.
template <class T>
std::ostream & operator <<(std::ostream & os, utility::vector1<T> const & v) {
	os << "[";
	for ( size_t i=1; i<=v.size(); i++ ) {
		os << v[i];
		if ( i < v.size() ) os << ", ";
	}
	os << "]";
	return os;
}

/// @brief Output function for std::pair object.
template <typename T1, typename T2>
std::ostream & operator <<(std::ostream & os, std::pair<T1, T2> const & v) {
	os << "[";
	os << v.first;
	os << ", ";
	os << v.second;
	os << "]";
	return os;
}


/// @brief Output function for std::map object.
template <typename T1, typename T2>
std::ostream & operator <<(std::ostream & os, std::map<T1, T2> const & m) {
	typedef typename std::map<T1, T2>::const_iterator ConstIterator;
	ConstIterator p;

	os << "{";

	for ( p=m.begin(); p!=m.end(); p++ ) {
		os << p->first << ":" << p->second << ", ";
	}

	os << "}";
	return os;
}


/// @brief Output function for std::list object.
template <typename T>
std::ostream & operator <<(std::ostream & os, std::list<T> const & l) {
	typedef typename std::list<T>::const_iterator ConstIterator;
	ConstIterator p;

	os << "[[";

	for ( p=l.begin(); p!=l.end(); p++ ) {
		os << *p << ", ";
	}

	os << "]]";
	return os;
}

} // namespace utility

#endif // INCLUDED_utility_stream_util_hh
