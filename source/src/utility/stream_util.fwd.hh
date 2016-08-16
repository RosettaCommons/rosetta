// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/stream_util.fwd.hh
/// @brief  Implemention of ostream operator << for various common types, forward declarations
/// @author Sergey Lyskov

#ifndef INCLUDED_utility_stream_util_fwd_hh
#define INCLUDED_utility_stream_util_fwd_hh

#include <utility/type_traits.fwd.hh>

#include <utility/vectorL.fwd.hh>
#include <utility/vector0.fwd.hh>
#include <utility/vector1.fwd.hh>

#include <vector>
#include <list>
#include <map>
#include <set>
#include <utility>

#ifdef CXX11

#include <type_traits>


namespace std { // inserting operator for ::std types in to std namespace

/// @brief Output function for std::map object.
template <typename T1, typename T2, typename std::enable_if< utility::has_insertion_operator_s<T1>::value  and  utility::has_insertion_operator_s<T2>::value >::type * = nullptr>
std::ostream & operator <<(std::ostream & os, std::map<T1, T2> const & m);


/// @brief Output function for std::list object.
template <typename T, typename std::enable_if< utility::has_insertion_operator_s<T>::value >::type * = nullptr>
std::ostream & operator <<(std::ostream & os, std::list<T> const & l);

/// @brief Output function for std::set object.
template <typename T, typename std::enable_if< utility::has_insertion_operator_s<T>::value >::type * = nullptr>
std::ostream & operator <<(std::ostream & os, std::set<T> const & s);

// forward declaration to allow output of composite types, like: vector1< std::pair >.
template <typename T1, typename T2, typename std::enable_if< utility::has_insertion_operator_s<T1>::value  and  utility::has_insertion_operator_s<T2>::value >::type * = nullptr>
std::ostream & operator <<(std::ostream & os, std::pair<T1, T2> const & v);

/// @brief Output function for std::vector object.
template <class T, typename std::enable_if< utility::has_insertion_operator_s<T>::value >::type * = nullptr>
std::ostream & operator <<( std::ostream & os, std::vector<T> const & v);

} // namespace std

namespace utility {
/// @brief Output function for utility::vectorL object.
template <platform::SSize L, class T, typename std::enable_if< utility::has_insertion_operator_s<T>::value >::type * = nullptr>
std::ostream & operator <<(std::ostream & os, utility::vectorL<L, T> const & v);

} // namespace utility

#endif



#endif // INCLUDED_utility_stream_util_fwd_hh
