// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/pointer/make_shared.cxxtest.hh
/// @brief  Rosetta make_shared tests
/// @author Sergey Lyskov

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/memory.hh>

#include <type_traits>

#include <cxxtest/TestSuite.h>

class make_shared_Tests : public CxxTest::TestSuite
{
public:

	void test_make_shared() {
		struct A {
			A(int, double const *, std::string &);
		};

		using A_SP = utility::pointer::shared_ptr<A>;

		TS_ASSERT( (std::is_same<A_SP, decltype( utility::pointer::make_shared<A>(1, nullptr, "1") ) >::value) );
	}

};
