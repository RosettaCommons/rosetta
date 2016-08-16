// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/BitSet.cxxtest.hh
/// @brief  BitSet.cxxtest: test suite for utility::BitSet
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Ron Jacak (ron.jacak@gmail.com)


// Package headers
#include <utility/BitSet.hh>
#include <cxxtest/TestSuite.h>


// Place the following definitions inside a namespace and then "use" that namespace in each of the
// individual test functions.  Without this namespace, the linker throws a multiple definition
// error for the operator "|" because there's one definition here and one in BitVector.cxxtest.hh
// with the same argument list but different return type. For some reason, the operator | is only
// defined here (perhaps so that the test cases are constructed very fast?) - not in the actual
// src/utility/BitSet.hh file.  One solution would be to replace the name of the function with
// something that would be unique to this file (and would only be defined once at link time) such as
// "setUnion" and do the same in BitVector.cxxtest.hh.  The alternative, which is what's done here,
// is to "separate" the two using a custom namespace. -ronj
//
namespace set {

enum Bit { zero, one, two, three, four, nine=9, ten };

inline
utility::BitSet< Bit >
operator |( Bit const & a, Bit const & b )
{
	return utility::BitSet< Bit >( a, b );
}

namespace b {

class Bit {
public:

	inline
	explicit
	Bit( int ii ) :
		i( ii )
	{}

	friend
	inline
	bool
	operator <( Bit const & a, Bit const & b )
	{
		return ( a.i < b.i );
	}

	friend
	inline
	utility::BitSet< Bit >
	operator |( Bit const & a, Bit const & b )
	{
		return utility::BitSet< Bit >( a, b );
	}

	int i;
};

} // namespace b
} // namespace vector


class BitSetTests : public CxxTest::TestSuite {

	public:

		/// @brief Constructor Test
		void test_BitSet_constructor() {

			using namespace set;

			{ // enum Bit: Copy construction
				utility::BitSet< Bit > s( one | nine |= ten );
				TS_ASSERT( s[ one ] );
				TS_ASSERT( ! s[ two ] );
				TS_ASSERT( s[ nine ] );
				TS_ASSERT( s[ ten ] );
				TS_ASSERT( s[ ten ] );
				TS_ASSERT( s.size() == 3 );
				TS_ASSERT( ! s.empty() );
			}

			{ // class Bit Bit: Copy construction
				b::Bit one(1), two(2), nine(9), ten(10);
				utility::BitSet< b::Bit > s( one | nine |= ten );
				TS_ASSERT( s[ one ] );
				TS_ASSERT( ! s[ two ] );
				TS_ASSERT( s[ nine ] );
				TS_ASSERT( s[ ten ] );
			}

			{ // enum Bit: Multiple bit construction
				utility::BitSet< Bit > s( one, nine, ten );
				TS_ASSERT( s[ one ] );
				TS_ASSERT( ! s[ two ] );
				TS_ASSERT( s[ nine ] );
				TS_ASSERT( s[ ten ] );
			}

		}


		/// @brief Assignment Test
		void test_BitSet_assignment() {

			using namespace set;

			{ // |= BitSet
				utility::BitSet< Bit > s( one | nine |= ten );
				utility::BitSet< Bit > const t( two );
				TS_ASSERT( s[ one ] );
				TS_ASSERT( ! s[ two ] );
				TS_ASSERT( s[ nine ] );
				TS_ASSERT( s[ ten ] );
				TS_ASSERT( t[ two ] );
				s |= t;
				TS_ASSERT( s[ two ] );
			}

			{ // += BitSet
				utility::BitSet< Bit > s( one | nine |= ten );
				utility::BitSet< Bit > const t( two );
				TS_ASSERT( s[ one ] );
				TS_ASSERT( ! s[ two ] );
				TS_ASSERT( s[ nine ] );
				TS_ASSERT( s[ ten ] );
				TS_ASSERT( t[ two ] );
				s += t;
				TS_ASSERT( s[ two ] );
			}

			{ // -= BitSet
				utility::BitSet< Bit > s( one | nine |= ten );
				utility::BitSet< Bit > const t( nine );
				TS_ASSERT( s[ one ] );
				TS_ASSERT( ! s[ two ] );
				TS_ASSERT( s[ nine ] );
				TS_ASSERT( s[ ten ] );
				TS_ASSERT( t[ nine ] );
				s -= t;
				TS_ASSERT( ! s[ nine ] );
			}

			{ // |=, -=, += Bit
				utility::BitSet< Bit > s( one | nine |= ten );
				TS_ASSERT( s[ one ] );
				TS_ASSERT( ! s[ two ] );
				TS_ASSERT( s[ nine ] );
				TS_ASSERT( s[ ten ] );
				s |= two;
				TS_ASSERT( s[ two ] );
				s -= nine;
				TS_ASSERT( ! s[ nine ] );
				s += nine;
				TS_ASSERT( s[ nine ] );
			}
		}


		/// @brief Method Test
		void test_BitSet_method() {

			using namespace set;

			{ // BitSet | BitSet
				utility::BitSet< Bit > const s( one | nine );
				utility::BitSet< Bit > const t( two | ten );
				utility::BitSet< Bit > const u( s | t );
				TS_ASSERT( u[ one ] );
				TS_ASSERT( u[ two ] );
				TS_ASSERT( u[ nine ] );
				TS_ASSERT( u[ ten ] );
				TS_ASSERT( ! u[ three ] );
			}

			{ // BitSet + BitSet
				utility::BitSet< Bit > const s( one | nine );
				utility::BitSet< Bit > const t( two | ten );
				utility::BitSet< Bit > const u( s + t );
				TS_ASSERT( ! s[ ten ] );
				TS_ASSERT( u[ one ] );
				TS_ASSERT( u[ two ] );
				TS_ASSERT( u[ nine ] );
				TS_ASSERT( u[ ten ] );
				TS_ASSERT( ! u[ three ] );
			}

			{ // BitSet - BitSet
				utility::BitSet< Bit > const s( one | nine );
				utility::BitSet< Bit > const t( two | nine );
				utility::BitSet< Bit > const u( s - t );
				TS_ASSERT( u[ one ] );
				TS_ASSERT( ! u[ two ] );
				TS_ASSERT( ! u[ nine ] );
				TS_ASSERT( ! u[ ten ] );
			}

			{ // Swap
				utility::BitSet< Bit > s( one | nine );
				utility::BitSet< Bit > t( two | nine );
				TS_ASSERT( s != t );
				TS_ASSERT( s[ one ] );
				TS_ASSERT( ! s[ two ] );
				TS_ASSERT( s[ nine ] );
				TS_ASSERT( t[ two ] );
				TS_ASSERT( ! t[ three ] );
				TS_ASSERT( t[ nine ] );
				swap( s, t );
				TS_ASSERT( t[ one ] );
				TS_ASSERT( ! t[ two ] );
				TS_ASSERT( t[ nine ] );
				TS_ASSERT( s[ two ] );
				TS_ASSERT( ! s[ three ] );
				TS_ASSERT( s[ nine ] );
			}

			{ // BitSet ==, != BitSet
				utility::BitSet< Bit > s( one | nine |= ten );
				utility::BitSet< Bit > t;
				TS_ASSERT( t.empty() );
				t = s;
				TS_ASSERT( s == t );
				t |= two;
				TS_ASSERT( s != t );
			}

		}

};


