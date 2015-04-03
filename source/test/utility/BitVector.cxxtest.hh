// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/BitVector.cxxtest.hh
/// @brief  BitVector.cxxtest: test suite for utility::BitVector
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Ron Jacak (ron.jacak@gmail.com)


// Package headers
#include <utility/BitVector.hh>
#include <cxxtest/TestSuite.h>


// See the comment in BitSet.cxxtest.hh for why this namespace is here. -ronj
namespace vector {

enum Bit { zero, one, two, three, four, nine=9, ten };

inline
utility::BitVector< Bit >
operator |( Bit const & a, Bit const & b )
{
	return utility::BitVector< Bit >( a, b );
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
	utility::BitVector< Bit >
	operator |( Bit const & a, Bit const & b )
	{
		return utility::BitVector< Bit >( a, b );
	}

	inline
	operator int() const
	{
		return i;
	}

	int i;
};

} // namespace b
} // namespace vector


class BitVectorTests : public CxxTest::TestSuite {

	public:

		/// @brief Constructor Test
		void test_BitVector_constructor() {

			using namespace vector;

			{ // enum Bit: Copy construction
				utility::BitVector< Bit > s( one | nine |= ten );
				TS_ASSERT( s[ one ] );
				TS_ASSERT( ! s[ two ] );
				TS_ASSERT( s[ nine ] );
				TS_ASSERT( s[ ten ] );
				TS_ASSERT( s.size() == 11 );
				TS_ASSERT( ! s.empty() );
			}

			{ // class Bit Bit: Copy construction
				b::Bit one(1), two(2), nine(9), ten(10);
				utility::BitVector< b::Bit > s( one | nine |= ten );
				TS_ASSERT( s[ one ] );
				TS_ASSERT( ! s[ two ] );
				TS_ASSERT( s[ nine ] );
				TS_ASSERT( s[ ten ] );
			}

			{ // enum Bit: Multiple bit construction
				utility::BitVector< Bit > s( one, nine, ten );
				TS_ASSERT( s[ one ] );
				TS_ASSERT( ! s[ two ] );
				TS_ASSERT( s[ nine ] );
				TS_ASSERT( s[ ten ] );
			}

		}

		/// @brief Assignment Test
		void test_BitVector_assignment() {

			using namespace vector;

			{ // |= BitVector
				utility::BitVector< Bit > s( one | nine |= ten );
				utility::BitVector< Bit > const t( two );
				TS_ASSERT( s[ one ] );
				TS_ASSERT( ! s[ two ] );
				TS_ASSERT( s[ nine ] );
				TS_ASSERT( s[ ten ] );
				TS_ASSERT( t[ two ] );
				s |= t;
				TS_ASSERT( s[ two ] );
			}

			{ // += BitVector
				utility::BitVector< Bit > s( one | nine |= ten );
				utility::BitVector< Bit > const t( two );
				TS_ASSERT( s[ one ] );
				TS_ASSERT( ! s[ two ] );
				TS_ASSERT( s[ nine ] );
				TS_ASSERT( s[ ten ] );
				TS_ASSERT( t[ two ] );
				s += t;
				TS_ASSERT( s[ two ] );
			}

			{ // -= BitVector
				utility::BitVector< Bit > s( one | nine |= ten );
				utility::BitVector< Bit > const t( nine );
				TS_ASSERT( s[ one ] );
				TS_ASSERT( ! s[ two ] );
				TS_ASSERT( s[ nine ] );
				TS_ASSERT( s[ ten ] );
				TS_ASSERT( t[ nine ] );
				s -= t;
				TS_ASSERT( ! s[ nine ] );
			}

			{ // |=, -=, += Bit
				utility::BitVector< Bit > s( one | nine |= ten );
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
		void test_BitVector_method() {

			using namespace vector;

			{ // BitVector | BitVector
				utility::BitVector< Bit > const s( one | nine );
				utility::BitVector< Bit > const t( two | ten );
				utility::BitVector< Bit > const u( s | t );
				TS_ASSERT( u[ one ] );
				TS_ASSERT( u[ two ] );
				TS_ASSERT( u[ nine ] );
				TS_ASSERT( u[ ten ] );
				TS_ASSERT( ! u[ three ] );
			}

			{ // BitVector + BitVector
				utility::BitVector< Bit > const s( one | nine );
				utility::BitVector< Bit > const t( two | ten );
				utility::BitVector< Bit > const u( s + t );
				TS_ASSERT( ! s[ ten ] );
				TS_ASSERT( u[ one ] );
				TS_ASSERT( u[ two ] );
				TS_ASSERT( u[ nine ] );
				TS_ASSERT( u[ ten ] );
				TS_ASSERT( ! u[ three ] );
			}

			{ // BitVector - BitVector
				utility::BitVector< Bit > const s( one | nine );
				utility::BitVector< Bit > const t( two | nine );
				utility::BitVector< Bit > const u( s - t );
				TS_ASSERT( u[ one ] );
				TS_ASSERT( ! u[ two ] );
				TS_ASSERT( ! u[ nine ] );
				TS_ASSERT( ! u[ ten ] );
			}

			{ // Swap
				utility::BitVector< Bit > s( one | nine );
				utility::BitVector< Bit > t( two | nine );
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

			{ // BitVector ==, != BitVector
				utility::BitVector< Bit > s( one | nine |= ten );
				utility::BitVector< Bit > t;
				TS_ASSERT( t.empty() );
				t = s;
				TS_ASSERT( s == t );
				t |= two;
				TS_ASSERT( s != t );
			}

		}

};  // class BitVectorTests

