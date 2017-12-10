// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/type_traits/type_traits.hh
/// @brief  Implemention of type traits related to Rosetta types
/// @author Sergey Lyskov

#ifndef INCLUDED_utility_type_traits_has_insertion_operator_hh
#define INCLUDED_utility_type_traits_has_insertion_operator_hh

#include <utility/type_traits/has_insertion_operator.fwd.hh>

#include <utility/stream_util.fwd.hh>

#include <ostream>
#include <utility>

namespace utility {
namespace type_traits {

namespace has_insertion_operator_implementation {
enum class False {};
struct any_type {
	template<typename T> any_type(T const&);
};
False operator<<(std::ostream const&, any_type const&);
} // namespace has_insertion_operator_implementation

template<typename T>
constexpr bool has_insertion_operator() {
	using namespace has_insertion_operator_implementation;
	return std::is_same< decltype(std::declval<std::ostream&>() << std::declval<T>()), std::ostream & >::value;
}
// Workaround for MSVC 2015
template<typename T>
struct has_insertion_operator_s {
	static const bool value = has_insertion_operator<T>();
};

} // namespace type_traits
} // namespace utility

#endif // INCLUDED_utility_type_traits_has_insertion_operator_hh
