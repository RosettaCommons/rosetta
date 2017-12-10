// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/type_traits.cxxtest.hh
/// @brief  Rosetta type_traits tests
/// @author Sergey Lyskov

#include <utility/type_traits/is_string_constructible.hh>

#include <cxxtest/TestSuite.h>

class type_traits_Tests : public CxxTest::TestSuite
{
public:
	/// @brief test for is_string_constructible type trait
	void test_is_string_constructible() {
		TS_ASSERT_EQUALS( utility::type_traits::is_string_constructible<char *>::value, true);
		TS_ASSERT_EQUALS( utility::type_traits::is_string_constructible<char const *>::value, true);

		TS_ASSERT_EQUALS( utility::type_traits::is_string_constructible<char (&)[]>::value, true);
		TS_ASSERT_EQUALS( utility::type_traits::is_string_constructible<char (&)[42]>::value, true);
		TS_ASSERT_EQUALS( utility::type_traits::is_string_constructible<char const (&)[42]>::value, true);

		TS_ASSERT_EQUALS( utility::type_traits::is_string_constructible<std::string>::value, true);
		TS_ASSERT_EQUALS( utility::type_traits::is_string_constructible<std::string const &>::value, true);

		struct A {};
		TS_ASSERT_EQUALS( utility::type_traits::is_string_constructible<A>::value, false);
		TS_ASSERT_EQUALS( utility::type_traits::is_string_constructible<int>::value, false);
		TS_ASSERT_EQUALS( utility::type_traits::is_string_constructible<char>::value, false);
	}
};
