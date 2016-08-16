// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/stream_util.hh
/// @brief  Implemention of ostream operator << for various common types
/// @author Sergey Lyskov

#ifndef INCLUDED_utility_stream_util_hh
#define INCLUDED_utility_stream_util_hh

#include <utility/stream_util.fwd.hh>

#include <utility/type_traits.hh>

#include <utility/vectorL.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <ostream>


/// -------------------------------------------------------------------
/// Predefined functions for Tracer IO
/// We originally moved them out of name spaces so we can use them right away - without specifying it.
/// Returned to utility namespace, as some compilers did not find the operator otherwise! -- rhiju


#ifdef CXX11

namespace utility {

/// @brief Output function for utility::vectorL object.
template <platform::SSize L, class T, typename std::enable_if< utility::has_insertion_operator_s<T>::value >::type *>
std::ostream & operator <<(std::ostream & os, utility::vectorL<L, T> const & v) {
	os << "[";
	if ( v.size() ) {
		for ( size_t i=v.l(); i<=v.u(); ++i ) {
			os << v[i];
			if ( i < v.u() ) os << ", ";
		}
	}
	os << "]";
	return os;
}

} // namespace utility


namespace std { // inserting operator for ::std types in to std namespace

/// @brief Output function for std::pair object.
template <typename T1, typename T2, typename std::enable_if< utility::has_insertion_operator_s<T1>::value  and  utility::has_insertion_operator_s<T2>::value >::type *>
std::ostream & operator <<(std::ostream & os, std::pair<T1, T2> const & v) {
	os << "(";
	os << v.first;
	os << ", ";
	os << v.second;
	os << ")";
	return os;
}


/// @brief Output function for std::map object.
template <typename T1, typename T2, typename std::enable_if< utility::has_insertion_operator_s<T1>::value  and  utility::has_insertion_operator_s<T2>::value >::type *>
std::ostream & operator <<(std::ostream & os, std::map<T1, T2> const & m) {
	typedef typename std::map<T1, T2>::const_iterator ConstIterator;
	ConstIterator p;

	os << "{";

	for ( p=m.begin(); p!=m.end(); ++p ) {
		os << p->first << ":" << p->second << ", ";
	}

	os << "}";
	return os;
}


/// @brief Output function for std::list object.
template <typename T, typename std::enable_if< utility::has_insertion_operator_s<T>::value >::type *>
std::ostream & operator <<(std::ostream & os, std::list<T> const & l) {
	typedef typename std::list<T>::const_iterator ConstIterator;
	ConstIterator p;

	os << "[[";

	for ( p=l.begin(); p!=l.end(); ++p ) {
		os << *p << ", ";
	}

	os << "]]";
	return os;
}

/// @brief Output function for std::set object.
template <typename T, typename std::enable_if< utility::has_insertion_operator_s<T>::value >::type *>
std::ostream & operator <<(std::ostream & os, std::set<T> const & s) {
	typedef typename std::set<T>::const_iterator ConstIterator;
	ConstIterator p;

	os << "[";

	for ( p=s.begin(); p!=s.end(); ++p ) {
		os << *p << ", ";
	}

	os << "]";
	return os;
}

/// @brief Output function for std::vector object.
template <class T, typename std::enable_if< utility::has_insertion_operator_s<T>::value >::type * >
std::ostream & operator <<( std::ostream & os, std::vector<T> const & v)
{
	os << "[";
	for ( size_t i = 0; i < v.size(); ++i ) {
		os << v[i];
		if ( i < v.size()-1 ) os << ", ";
	}
	os << "]";
	return os;
}

} // namespace std

#else // CXX11

// declare functions first to allow output of vector1< std types > and vice/versa
namespace utility {

/// @brief Output function for utility::vector1 object.
template <platform::SSize L, class T>
std::ostream & operator <<(std::ostream & os, utility::vectorL<L, T> const & v);

}

namespace std {

template <typename T1, typename T2>
std::ostream & operator <<(std::ostream & os, std::pair<T1, T2> const & v);

/// @brief Output function for std::map object.
template <typename T1, typename T2>
std::ostream & operator <<(std::ostream & os, std::map<T1, T2> const & m );

/// @brief Output function for std::list object.
template <typename T>
std::ostream & operator <<(std::ostream & os, std::list<T> const & l);

/// @brief Output function for std::set object.
template <typename T>
std::ostream & operator <<(std::ostream & os, std::set<T> const & s);

/// @brief Output function for std::vector object.
template <class T>
std::ostream & operator <<( std::ostream & os, std::vector<T> const & v);

}

namespace utility {

/// @brief Output function for utility::vector1 object.
template <platform::SSize L, class T>
std::ostream & operator <<(std::ostream & os, utility::vectorL<L, T> const & v) {
	os << "[";
	if ( v.size() ) {
		for ( size_t i=v.l(); i<=v.u(); ++i ) {
			os << v[i];
			if ( i < v.u() ) os << ", ";
		}
	}
	os << "]";
	return os;
}

} // namespace utility


namespace std { // inserting operator for ::std types in to std namespace

/// @brief Output function for std::pair object.
template <typename T1, typename T2>
std::ostream & operator <<(std::ostream & os, std::pair<T1, T2> const & v) {
	os << "(";
	os << v.first;
	os << ", ";
	os << v.second;
	os << ")";
	return os;
}


/// @brief Output function for std::map object.
template <typename T1, typename T2>
std::ostream & operator <<(std::ostream & os, std::map<T1, T2> const & m) {
	typedef typename std::map<T1, T2>::const_iterator ConstIterator;
	ConstIterator p;

	os << "{";

	for ( p=m.begin(); p!=m.end(); ++p ) {
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

	for ( p=l.begin(); p!=l.end(); ++p ) {
		os << *p << ", ";
	}

	os << "]]";
	return os;
}

/// @brief Output function for std::set object.
template <typename T>
std::ostream & operator <<(std::ostream & os, std::set<T> const & s) {
	typedef typename std::set<T>::const_iterator ConstIterator;
	ConstIterator p;

	os << "[";

	for ( p=s.begin(); p!=s.end(); ++p ) {
		if ( p != s.begin() ) os << ", ";
		os << *p;
	}

	os << "]";
	return os;
}

/// @brief Output function for std::vector object.
template <class T>
std::ostream & operator <<( std::ostream & os, std::vector<T> const & v)
{
	os << "[";
	for ( size_t i = 0; i < v.size(); ++i ) {
		os << v[i];
		if ( i < v.size()-1 ) os << ", ";
	}
	os << "]";
	return os;
}

} // namespace std


#endif // CXX11


#endif // INCLUDED_utility_stream_util_hh
