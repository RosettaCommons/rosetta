// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/match/Hasher.cxxtest.hh
/// @brief  test suite for protocols::match::SixDHasher
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

// Utility headers
#include <utility/fixedsizearray1.hh>
#include <utility/exit.hh>

/// Project headers
#include <core/types.hh>
#include <protocols/match/MatchSet.hh>

// C++ headers
#include <string>
#include <iostream>

#include <numeric/HomogeneousTransform.hh>
#include <numeric/geometry/hashing/SixDHasher.hh>

// Boost headers
#include <boost/unordered_map.hpp>
#include <boost/algorithm/string/predicate.hpp>

//Auto Headers
#include <protocols/match/Hit.hh>
#include <utility/vector1.hh>
#include <boost/unordered/unordered_map_fwd.hpp>


// The code in namespace hash_examples are copy-righted under the boost liscense, as is
// struct word_info and the testing code contained in the unit test test_boost_unordered_map.
//
// Copyright 2006-2008 Daniel James.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
struct word_info {
	int tag;
	explicit word_info(int t = 0) : tag(t) {}
};

namespace hash_examples
{
struct iequal_to
	: std::binary_function<std::string, std::string, bool>
{
	iequal_to() {}
	explicit iequal_to(std::locale const& l) : locale_(l) {}

	template <typename String1, typename String2>
	bool operator()(String1 const& x1, String2 const& x2) const
	{
		return boost::algorithm::iequals(x1, x2, locale_);
	}
private:
	std::locale locale_;
};

struct ihash
	: std::unary_function<std::string, std::size_t>
{
	ihash() {}
	explicit ihash(std::locale const& l) : locale_(l) {}

	template <typename String>
	std::size_t operator()(String const& x) const
	{
		std::size_t seed = 0;

		for ( typename String::const_iterator it = x.begin();
				it != x.end(); ++it ) {
			boost::hash_combine(seed, std::toupper(*it, locale_));
		}

		return seed;
	}
private:
	std::locale locale_;
};
}


// --------------- Test Class --------------- //

class HasherTests : public CxxTest::TestSuite {

public:

	typedef core::Vector Vector;
	typedef core::Size   Size;
	typedef core::Real   Real;
	typedef numeric::geometry::BoundingBox< Vector > BoundingBox;
	typedef numeric::HomogeneousTransform< Real >    HTReal;


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	//FIX THIS ASAP

	// --------------- Test Cases --------------- //
	void test_boost_unordered_map() {

		boost::unordered_map<std::string, word_info,
			hash_examples::ihash, hash_examples::iequal_to> idictionary;

		TS_ASSERT(idictionary.empty());

		idictionary["one"] = word_info(1);
		TS_ASSERT(idictionary.size() == 1);
		TS_ASSERT(idictionary.find("ONE") != idictionary.end() &&
			idictionary.find("ONE") == idictionary.find("one"));

		idictionary.insert(std::make_pair("ONE", word_info(2)));
		TS_ASSERT(idictionary.size() == 1);
		TS_ASSERT(idictionary.find("ONE") != idictionary.end() &&
			idictionary.find("ONE")->first == "one" &&
			idictionary.find("ONE")->second.tag == 1);

		idictionary["One"] = word_info(3);
		TS_ASSERT(idictionary.size() == 1);
		TS_ASSERT(idictionary.find("ONE") != idictionary.end() &&
			idictionary.find("ONE")->first == "one" &&
			idictionary.find("ONE")->second.tag == 3);

		idictionary["two"] = word_info(4);
		TS_ASSERT(idictionary.size() == 2);
		TS_ASSERT(idictionary.find("two") != idictionary.end() &&
			idictionary.find("TWO")->first == "two" &&
			idictionary.find("Two")->second.tag == 4);

	}

	void test_SixDCoordinateBinner_offset_theta() {
		using namespace protocols::match;

		Vector lower( 12.5, 16.25, 4.25 );
		Vector upper( 15.5, 20, 8.5 );
		BoundingBox bb( lower, upper );
		Real6 binwidths;
		binwidths[ 1 ] = binwidths[ 2 ] = binwidths[ 3 ] = 0.25;
		binwidths[ 4 ] = binwidths[ 5 ] = binwidths[ 6 ] = 10;

		Size3 euler_offsets( 0 );

		numeric::geometry::hashing::SixDCoordinateBinner binner( bb, euler_offsets, binwidths );

		Real6 p;
		p[ 1 ] = 13.6; p[ 2 ] = 19.4; p[ 3 ] = 5.3;
		p[ 4 ] = 50;   p[ 5 ] = 123;  p[ 6 ] = 76;

		numeric::geometry::hashing::Bin6D bin = binner.bin6( p );
		//std::cout << "bin: ";
		//for ( Size ii = 1; ii <= 6; ++ii ) std::cout << bin[ ii ] << " ";
		//std::cout << std::endl;
		//bin:  4 12 4 5 12 7
		TS_ASSERT( bin[ 1 ] == 4 );
		TS_ASSERT( bin[ 2 ] == 12 );
		TS_ASSERT( bin[ 3 ] == 4 );
		TS_ASSERT( bin[ 4 ] == 5 );
		TS_ASSERT( bin[ 5 ] == 12 );
		TS_ASSERT( bin[ 6 ] == 7 );

		HTReal zrot1, xrot2A, xrot2B, xrot2C, zrot3;

		zrot1.set_zaxis_rotation_deg( 155 );
		xrot2A.set_xaxis_rotation_deg( 88 );
		xrot2B.set_xaxis_rotation_deg( 91 );
		xrot2C.set_xaxis_rotation_deg( 96 );
		zrot3.set_zaxis_rotation_deg( 44 );

		HTReal tA = zrot1 * xrot2A * zrot3;
		HTReal tB = zrot1 * xrot2B * zrot3;
		HTReal tC = zrot1 * xrot2C * zrot3;

		Vector eulerA = tA.euler_angles_deg(), eulerB = tB.euler_angles_deg(), eulerC = tC.euler_angles_deg();

		Real6 pA = p;
		pA[ 4 ] = eulerA(1);   pA[ 5 ] = eulerA(2);  pA[ 6 ] = eulerA(3);
		numeric::geometry::hashing::Bin6D binA = binner.bin6( pA );
		//std::cout << "binA: ";
		//for ( Size ii = 1; ii <= 6; ++ii ) std::cout << bin[ ii ] << " ";
		//std::cout << std::endl;
		boost::uint64_t indexA = binner.bin_index( pA );
		//std::cout << "indexA: " << indexA << std::endl;

		Real6 pB = p;
		pB[ 4 ] = eulerB(1);   pB[ 5 ] = eulerB(2);  pB[ 6 ] = eulerB(3);
		numeric::geometry::hashing::Bin6D binB = binner.bin6( pB );
		//std::cout << "binB: ";
		//for ( Size ii = 1; ii <= 6; ++ii ) std::cout << bin[ ii ] << " ";
		//std::cout << std::endl;
		boost::uint64_t indexB = binner.bin_index( pB );
		//std::cout << "indexB: " << indexB << std::endl;

		Real6 pC = p;
		pC[ 4 ] = eulerC(1);   pC[ 5 ] = eulerC(2);  pC[ 6 ] = eulerC(3);
		numeric::geometry::hashing::Bin6D binC = binner.bin6( pC );
		//std::cout << "binC: ";
		//for ( Size ii = 1; ii <= 6; ++ii ) std::cout << bin[ ii ] << " ";
		//std::cout << std::endl;
		boost::uint64_t indexC = binner.bin_index( pC );
		//std::cout << "indexC: " << indexC << std::endl;

		TS_ASSERT( indexA != indexB );
		TS_ASSERT( indexA != indexC );
		TS_ASSERT( indexB == indexC );

		numeric::geometry::hashing::Bin6D recoverA = binner.bin_from_index( indexA );
		numeric::geometry::hashing::Bin6D recoverB = binner.bin_from_index( indexB );
		numeric::geometry::hashing::Bin6D recoverC = binner.bin_from_index( indexC );

		for ( Size ii = 1; ii <= 6; ++ii ) { TS_ASSERT( recoverA[ ii ] == binA[ ii ] ); }
		for ( Size ii = 1; ii <= 6; ++ii ) { TS_ASSERT( recoverB[ ii ] == binB[ ii ] ); }
		for ( Size ii = 1; ii <= 6; ++ii ) { TS_ASSERT( recoverC[ ii ] == binC[ ii ] ); }


		Size3 euler_offsets2( 0 );
		euler_offsets2[ 3 ] = 1;

		numeric::geometry::hashing::SixDCoordinateBinner binner2( bb, euler_offsets2, binwidths );

		boost::uint64_t indexA2 = binner2.bin_index( pA );
		boost::uint64_t indexB2 = binner2.bin_index( pB );
		boost::uint64_t indexC2 = binner2.bin_index( pC );

		TS_ASSERT( indexA2 == indexB2 );
		TS_ASSERT( indexA2 != indexC2 );
		TS_ASSERT( indexB2 != indexC2 );

		//std::cout << "recover A ";
		//for ( Size ii = 1; ii <= 6; ++ii ) { std::cout << " " << recoverA[ ii ]; }
		//std::cout << std::endl << "recover B ";
		//for ( Size ii = 1; ii <= 6; ++ii ) { std::cout << " " << recoverB[ ii ]; }
		//std::cout << std::endl << "recover C ";
		//for ( Size ii = 1; ii <= 6; ++ii ) { std::cout << " " << recoverC[ ii ]; }
		//std::cout << std::endl;

		//std::cout << "tA: " << std::endl;
		//std::cout << "  " << tA.xx() << " " << tA.yx() << " " << tA.zx() << " " << tA.px() << std::endl;
		//std::cout << "  " << tA.xy() << " " << tA.yy() << " " << tA.zy() << " " << tA.py() << std::endl;
		//std::cout << "  " << tA.xz() << " " << tA.yz() << " " << tA.zz() << " " << tA.pz() << std::endl;

		//std::cout << "tB: " << std::endl;
		//std::cout << "  " << tB.xx() << " " << tB.yx() << " " << tB.zx() << " " << tB.px() << std::endl;
		//std::cout << "  " << tB.xy() << " " << tB.yy() << " " << tB.zy() << " " << tB.py() << std::endl;
		//std::cout << "  " << tB.xz() << " " << tB.yz() << " " << tB.zz() << " " << tB.pz() << std::endl;

		//std::cout << "tC: " << std::endl;
		//std::cout << "  " << tC.xx() << " " << tC.yx() << " " << tC.zx() << " " << tC.px() << std::endl;
		//std::cout << "  " << tC.xy() << " " << tC.yy() << " " << tC.zy() << " " << tC.py() << std::endl;
		//std::cout << "  " << tC.xz() << " " << tC.yz() << " " << tC.zz() << " " << tC.pz() << std::endl;


		//HTReal tAprime, tBprime, tCprime;
		//tAprime.from_euler_angles_deg( eulerA );
		//tBprime.from_euler_angles_deg( eulerB );
		//tCprime.from_euler_angles_deg( eulerC );

		//std::cout << "tAprime: " << std::endl;
		//std::cout << "  " << tAprime.xx() << " " << tAprime.yx() << " " << tAprime.zx() << " " << tAprime.px() << std::endl;
		//std::cout << "  " << tAprime.xy() << " " << tAprime.yy() << " " << tAprime.zy() << " " << tAprime.py() << std::endl;
		//std::cout << "  " << tAprime.xz() << " " << tAprime.yz() << " " << tAprime.zz() << " " << tAprime.pz() << std::endl;

		//std::cout << "tBprime: " << std::endl;
		//std::cout << "  " << tBprime.xx() << " " << tBprime.yx() << " " << tBprime.zx() << " " << tBprime.px() << std::endl;
		//std::cout << "  " << tBprime.xy() << " " << tBprime.yy() << " " << tBprime.zy() << " " << tBprime.py() << std::endl;
		//std::cout << "  " << tBprime.xz() << " " << tBprime.yz() << " " << tBprime.zz() << " " << tBprime.pz() << std::endl;

		//std::cout << "tCprime: " << std::endl;
		//std::cout << "  " << tCprime.xx() << " " << tCprime.yx() << " " << tCprime.zx() << " " << tCprime.px() << std::endl;
		//std::cout << "  " << tCprime.xy() << " " << tCprime.yy() << " " << tCprime.zy() << " " << tCprime.py() << std::endl;
		//std::cout << "  " << tCprime.xz() << " " << tCprime.yz() << " " << tCprime.zz() << " " << tCprime.pz() << std::endl;


		//std::cout << "EulerA: ";
		//for ( Size ii = 1; ii <= 3; ++ii ) std::cout << eulerA( ii ) << " ";
		//std::cout << "EulerB: ";
		//for ( Size ii = 1; ii <= 3; ++ii ) std::cout << eulerB( ii ) << " ";
		//std::cout << std::endl;
		//std::cout << "EulerC: ";
		//for ( Size ii = 1; ii <= 3; ++ii ) std::cout << eulerC( ii ) << " ";
		//std::cout << std::endl;

		/* FIX THIS ASAP
		using namespace core;

		xyzbin_hash hash_12_15_10( 12, 15, 10 );
		xyzbin_equals equals;

		boost::unordered_map< xyzbin, std::string, xyzbin_hash, xyzbin_equals >
		volume_hash( 100, hash_12_15_10 );

		xyzbin bin;
		bin[ 1 ] = 3; bin[ 2 ] = 4; bin[ 3 ] = 8;
		Size hashed_val = hash_12_15_10( bin );
		TS_ASSERT( hashed_val == 3 + 12*4 + 12*15*8 );

		volume_hash[ bin ] = "w00t";
		TS_ASSERT( "w00t" == volume_hash.find( bin )->second );

		typedef boost::unordered_map< xyzbin, std::string, xyzbin_hash, xyzbin_equals >::const_iterator volume_iter;

		volume_iter iter = volume_hash.begin();
		TS_ASSERT( equals( iter->first, bin ) );
		++iter;
		TS_ASSERT( iter == volume_hash.end() );
		*/
	}


	void test_SixDCoordinateBinner_offset_psi() {
		using namespace protocols::match;

		Vector lower( 12.5, 16.25, 4.25 );
		Vector upper( 15.5, 20, 8.5 );
		BoundingBox bb( lower, upper );
		Real6 binwidths;
		binwidths[ 1 ] = binwidths[ 2 ] = binwidths[ 3 ] = 0.25;
		binwidths[ 4 ] = binwidths[ 5 ] = binwidths[ 6 ] = 10;

		Size3 euler_offsets( 0 );

		numeric::geometry::hashing::SixDCoordinateBinner binner( bb, euler_offsets, binwidths );

		Real6 pA, pB, pC;
		pA[ 1 ] = 13.6; pA[ 2 ] = 19.4; pA[ 3 ] = 5.3;
		pA[ 4 ] = 50;   pA[ 5 ] = 127;  pA[ 6 ] = 76;

		pB[ 1 ] = 13.6; pB[ 2 ] = 19.4; pB[ 3 ] = 5.3;
		pB[ 4 ] = 50;   pB[ 5 ] = 130;  pB[ 6 ] = 76;

		pC[ 1 ] = 13.6; pC[ 2 ] = 19.4; pC[ 3 ] = 5.3;
		pC[ 4 ] = 50;   pC[ 5 ] = 136;  pC[ 6 ] = 76;

		boost::uint64_t indexA = binner.bin_index( pA );
		boost::uint64_t indexB = binner.bin_index( pB );
		boost::uint64_t indexC = binner.bin_index( pC );

		TS_ASSERT( indexA != indexB );
		TS_ASSERT( indexA != indexC );
		TS_ASSERT( indexB == indexC );


		Size3 euler_offsets2( 0 );
		euler_offsets2[ 2 ] = 1;

		numeric::geometry::hashing::SixDCoordinateBinner binner2( bb, euler_offsets2, binwidths );

		boost::uint64_t indexA2 = binner2.bin_index( pA );
		boost::uint64_t indexB2 = binner2.bin_index( pB );
		boost::uint64_t indexC2 = binner2.bin_index( pC );

		TS_ASSERT( indexA2 == indexB2 );
		TS_ASSERT( indexA2 != indexC2 );
		TS_ASSERT( indexB2 != indexC2 );

	}

	void test_SixDCoordinateBinner_offset_phi() {

		using namespace protocols::match;

		Vector lower( 12.5, 16.25, 4.25 );
		Vector upper( 15.5, 20, 8.5 );
		BoundingBox bb( lower, upper );
		Real6 binwidths;
		binwidths[ 1 ] = binwidths[ 2 ] = binwidths[ 3 ] = 0.25;
		binwidths[ 4 ] = binwidths[ 5 ] = binwidths[ 6 ] = 10;

		Size3 euler_offsets( 0 );

		numeric::geometry::hashing::SixDCoordinateBinner binner( bb, euler_offsets, binwidths );

		Real6 pA, pB, pC;
		pA[ 1 ] = 13.6; pA[ 2 ] = 19.4; pA[ 3 ] = 5.3;
		pA[ 4 ] = 49;   pA[ 5 ] = 127;  pA[ 6 ] = 176;

		pB[ 1 ] = 13.6; pB[ 2 ] = 19.4; pB[ 3 ] = 5.3;
		pB[ 4 ] = 50;   pB[ 5 ] = 127;  pB[ 6 ] = 176;

		pC[ 1 ] = 13.6; pC[ 2 ] = 19.4; pC[ 3 ] = 5.3;
		pC[ 4 ] = 57;   pC[ 5 ] = 127;  pC[ 6 ] = 176;

		boost::uint64_t indexA = binner.bin_index( pA );
		boost::uint64_t indexB = binner.bin_index( pB );
		boost::uint64_t indexC = binner.bin_index( pC );

		TS_ASSERT( indexA != indexB );
		TS_ASSERT( indexA != indexC );
		TS_ASSERT( indexB == indexC );


		Size3 euler_offsets2( 0 );
		euler_offsets2[ 1 ] = 1;

		numeric::geometry::hashing::SixDCoordinateBinner binner2( bb, euler_offsets2, binwidths );

		boost::uint64_t indexA2 = binner2.bin_index( pA );
		boost::uint64_t indexB2 = binner2.bin_index( pB );
		boost::uint64_t indexC2 = binner2.bin_index( pC );

		TS_ASSERT( indexA2 == indexB2 );
		TS_ASSERT( indexA2 != indexC2 );
		TS_ASSERT( indexB2 != indexC2 );

	}

	void test_SixDCoordinateBinner_wrap_theta() {
		using namespace protocols::match;

		Vector lower( 12.5, 16.25, 4.25 );
		Vector upper( 15.5, 20, 8.5 );
		BoundingBox bb( lower, upper );
		Real6 binwidths;
		binwidths[ 1 ] = binwidths[ 2 ] = binwidths[ 3 ] = 0.25;
		binwidths[ 4 ] = binwidths[ 5 ] = binwidths[ 6 ] = 10;

		Size3 euler_offsets( 0 );

		numeric::geometry::hashing::SixDCoordinateBinner binner( bb, euler_offsets, binwidths );

		Real6 pA, pB, pC;
		pA[ 1 ] = 13.6; pA[ 2 ] = 19.4; pA[ 3 ] = 5.3;
		pB[ 1 ] = 13.6; pB[ 2 ] = 19.4; pB[ 3 ] = 5.3;
		pC[ 1 ] = 13.6; pC[ 2 ] = 19.4; pC[ 3 ] = 5.3;

		HTReal zrot1, xrot2A, xrot2B, xrot2C, zrot3;

		zrot1.set_zaxis_rotation_deg( 155 );
		xrot2A.set_xaxis_rotation_deg( 173 );
		xrot2B.set_xaxis_rotation_deg( 177 );
		xrot2C.set_xaxis_rotation_deg( 182 );
		zrot3.set_zaxis_rotation_deg( 44 );

		HTReal tA = zrot1 * xrot2A * zrot3;
		HTReal tB = zrot1 * xrot2B * zrot3;
		HTReal tC = zrot1 * xrot2C * zrot3;

		Vector eulerA = tA.euler_angles_deg(), eulerB = tB.euler_angles_deg(), eulerC = tC.euler_angles_deg();
		for ( Size ii = 1; ii <= 2; ++ii ) if ( eulerA( ii ) < 0 ) eulerA( ii ) += 360;
		for ( Size ii = 1; ii <= 2; ++ii ) if ( eulerB( ii ) < 0 ) eulerB( ii ) += 360;
		for ( Size ii = 1; ii <= 2; ++ii ) if ( eulerC( ii ) < 0 ) eulerC( ii ) += 360;

		//for ( Size ii = 1; ii <= 3; ++ii ) std::cout << " " << eulerA( ii ); std::cout << std::endl;
		//for ( Size ii = 1; ii <= 3; ++ii ) std::cout << " " << eulerB( ii ); std::cout << std::endl;
		//for ( Size ii = 1; ii <= 3; ++ii ) std::cout << " " << eulerC( ii ); std::cout << std::endl;

		pA[ 4 ] = eulerA(1);   pA[ 5 ] = eulerA(2);  pA[ 6 ] = eulerA(3);
		numeric::geometry::hashing::Bin6D binA = binner.bin6( pA );
		//std::cout << "binA: ";
		//for ( Size ii = 1; ii <= 6; ++ii ) std::cout << binA[ ii ] << " "; std::cout << std::endl;
		boost::uint64_t indexA = binner.bin_index( pA );
		//std::cout << "indexA: " << indexA << std::endl;

		pB[ 4 ] = eulerB(1);   pB[ 5 ] = eulerB(2);  pB[ 6 ] = eulerB(3);
		numeric::geometry::hashing::Bin6D binB = binner.bin6( pB );
		//std::cout << "binB: ";
		//for ( Size ii = 1; ii <= 6; ++ii ) std::cout << binB[ ii ] << " "; std::cout << std::endl;
		boost::uint64_t indexB = binner.bin_index( pB );
		//std::cout << "indexB: " << indexB << std::endl;

		pC[ 4 ] = eulerC(1);   pC[ 5 ] = eulerC(2);  pC[ 6 ] = eulerC(3);
		numeric::geometry::hashing::Bin6D binC = binner.bin6( pC );
		//std::cout << "binC: ";
		//for ( Size ii = 1; ii <= 6; ++ii ) std::cout << binC[ ii ] << " "; std::cout << std::endl;
		boost::uint64_t indexC = binner.bin_index( pC );
		//std::cout << "indexC: " << indexC << std::endl;

		//binA: 4 12 4 15 4 17
		//indexA: 28656593
		//binB: 4 12 4 15 4 17
		//indexB: 28656593
		//binC: 4 12 4 33 22 17
		//indexC: 28668581
		TS_ASSERT( binA[ 1 ] == 4 );  TS_ASSERT( binA[ 2 ] == 12 ); TS_ASSERT( binA[ 3 ] == 4 );
		TS_ASSERT( binA[ 4 ] == 15 ); TS_ASSERT( binA[ 5 ] == 4 );  TS_ASSERT( binA[ 6 ] == 17 );
		TS_ASSERT( indexA == 28656593 );

		TS_ASSERT( binB[ 1 ] == 4 );  TS_ASSERT( binB[ 2 ] == 12 ); TS_ASSERT( binB[ 3 ] == 4 );
		TS_ASSERT( binB[ 4 ] == 15 ); TS_ASSERT( binB[ 5 ] == 4 );  TS_ASSERT( binB[ 6 ] == 17 );
		TS_ASSERT( indexB == 28656593 );

		TS_ASSERT( binC[ 1 ] == 4 );  TS_ASSERT( binC[ 2 ] == 12 ); TS_ASSERT( binC[ 3 ] == 4 );
		TS_ASSERT( binC[ 4 ] == 33 ); TS_ASSERT( binC[ 5 ] == 22 ); TS_ASSERT( binC[ 6 ] == 17 );
		TS_ASSERT( indexC == 28668581 );

		//boost::uint64_t indexA = binner.bin_index( pA );
		//boost::uint64_t indexB = binner.bin_index( pB );
		//boost::uint64_t indexC = binner.bin_index( pC );

		TS_ASSERT( indexA == indexB );
		TS_ASSERT( indexA != indexC );
		TS_ASSERT( indexB != indexC );


		Size3 euler_offsets2( 0 );
		euler_offsets2[ 3 ] = 1;

		numeric::geometry::hashing::SixDCoordinateBinner binner2( bb, euler_offsets2, binwidths );
		//std::cout << "binner2 dimsizes:";
		//for ( Size ii = 1; ii <= 6; ++ii ) std::cout << " " << binner2.dimsizes()[ ii ]; std::cout << std::endl;

		boost::uint64_t indexA2 = binner2.bin_index( pA );
		boost::uint64_t indexB2 = binner2.bin_index( pB );
		boost::uint64_t indexC2 = binner2.bin_index( pC );

		TS_ASSERT( indexA2 != indexB2 );
		TS_ASSERT( indexA2 != indexC2 );
		TS_ASSERT( indexB2 == indexC2 );

		numeric::geometry::hashing::Bin6D binA2 = binner2.bin6( pA );
		//std::cout << "binA2: ";
		//for ( Size ii = 1; ii <= 6; ++ii ) std::cout << binA2[ ii ] << " "; std::cout << std::endl;
		//std::cout << "indexA2: " << indexA2 << std::endl;

		numeric::geometry::hashing::Bin6D binB2 = binner2.bin6( pB );
		//std::cout << "binB2: ";
		//for ( Size ii = 1; ii <= 6; ++ii ) std::cout << binB2[ ii ] << " "; std::cout << std::endl;
		//std::cout << "indexB2: " << indexB2 << std::endl;

		numeric::geometry::hashing::Bin6D binC2 = binner2.bin6( pC );
		//std::cout << "binC2: ";
		//for ( Size ii = 1; ii <= 6; ++ii ) std::cout << binC2[ ii ] << " "; std::cout << std::endl;
		//std::cout << "indexC2: " << indexC2 << std::endl;


		TS_ASSERT( binA2[ 1 ] == 4 );  TS_ASSERT( binA2[ 2 ] == 12 ); TS_ASSERT( binA2[ 3 ] == 4 );
		TS_ASSERT( binA2[ 4 ] == 15 ); TS_ASSERT( binA2[ 5 ] == 4 );  TS_ASSERT( binA2[ 6 ] == 17 );
		TS_ASSERT( indexA2 == 30248625 );

		TS_ASSERT( binB2[ 1 ] == 4 );  TS_ASSERT( binB2[ 2 ] == 12 ); TS_ASSERT( binB2[ 3 ] == 4 );
		TS_ASSERT( binB2[ 4 ] == 15 ); TS_ASSERT( binB2[ 5 ] == 4 );  TS_ASSERT( binB2[ 6 ] == 18 );
		TS_ASSERT( indexB2 == 30248626 );

		TS_ASSERT( binC2[ 1 ] == 4 );  TS_ASSERT( binC2[ 2 ] == 12 ); TS_ASSERT( binC2[ 3 ] == 4 );
		TS_ASSERT( binC2[ 4 ] == 15 ); TS_ASSERT( binC2[ 5 ] == 4 );  TS_ASSERT( binC2[ 6 ] == 18 );
		TS_ASSERT( indexC2 == 30248626 );

		//binA2: 4 12 4 15 4 17
		//indexA2: 28656593
		//binB2: 4 12 4 15 4 18
		//indexB2: 30248626
		//binC2: 4 12 4 15 4 18
		//indexC2: 30248626

	}

	void test_SixDCoordinateBinner_wrap_theta2() {
		using namespace protocols::match;

		Vector lower( 12.5, 16.25, 4.25 );
		Vector upper( 15.5, 20, 8.5 );
		BoundingBox bb( lower, upper );
		Real6 binwidths;
		binwidths[ 1 ] = binwidths[ 2 ] = binwidths[ 3 ] = 0.25;
		binwidths[ 4 ] = binwidths[ 5 ] = binwidths[ 6 ] = 10;

		Size3 euler_offsets( 0 );

		numeric::geometry::hashing::SixDCoordinateBinner binner( bb, euler_offsets, binwidths );

		Real6 pA, pB, pC;
		pA[ 1 ] = 13.6; pA[ 2 ] = 19.4; pA[ 3 ] = 5.3;
		pB[ 1 ] = 13.6; pB[ 2 ] = 19.4; pB[ 3 ] = 5.3;
		pC[ 1 ] = 13.6; pC[ 2 ] = 19.4; pC[ 3 ] = 5.3;

		HTReal zrot1, xrot2A, xrot2B, xrot2C, zrot3;

		zrot1.set_zaxis_rotation_deg( 155 );
		xrot2A.set_xaxis_rotation_deg( 6 );
		xrot2B.set_xaxis_rotation_deg( 3 );
		xrot2C.set_xaxis_rotation_deg( -2 );
		zrot3.set_zaxis_rotation_deg( 44 );

		HTReal tA = zrot1 * xrot2A * zrot3;
		HTReal tB = zrot1 * xrot2B * zrot3;
		HTReal tC = zrot1 * xrot2C * zrot3;

		Vector eulerA = tA.euler_angles_deg(), eulerB = tB.euler_angles_deg(), eulerC = tC.euler_angles_deg();
		for ( Size ii = 1; ii <= 2; ++ii ) if ( eulerA( ii ) < 0 ) eulerA( ii ) += 360;
		for ( Size ii = 1; ii <= 2; ++ii ) if ( eulerB( ii ) < 0 ) eulerB( ii ) += 360;
		for ( Size ii = 1; ii <= 2; ++ii ) if ( eulerC( ii ) < 0 ) eulerC( ii ) += 360;

		//for ( Size ii = 1; ii <= 3; ++ii ) std::cout << " " << eulerA( ii ); std::cout << std::endl;
		//for ( Size ii = 1; ii <= 3; ++ii ) std::cout << " " << eulerB( ii ); std::cout << std::endl;
		//for ( Size ii = 1; ii <= 3; ++ii ) std::cout << " " << eulerC( ii ); std::cout << std::endl;

		pA[ 4 ] = eulerA(1);   pA[ 5 ] = eulerA(2);  pA[ 6 ] = eulerA(3);
		numeric::geometry::hashing::Bin6D binA = binner.bin6( pA );
		//std::cout << "binA: ";
		//for ( Size ii = 1; ii <= 6; ++ii ) std::cout << binA[ ii ] << " "; std::cout << std::endl;
		boost::uint64_t indexA = binner.bin_index( pA );
		//std::cout << "indexA: " << indexA << std::endl;

		pB[ 4 ] = eulerB(1);   pB[ 5 ] = eulerB(2);  pB[ 6 ] = eulerB(3);
		numeric::geometry::hashing::Bin6D binB = binner.bin6( pB );
		//std::cout << "binB: ";
		//for ( Size ii = 1; ii <= 6; ++ii ) std::cout << binB[ ii ] << " "; std::cout << std::endl;
		boost::uint64_t indexB = binner.bin_index( pB );
		//std::cout << "indexB: " << indexB << std::endl;

		pC[ 4 ] = eulerC(1);   pC[ 5 ] = eulerC(2);  pC[ 6 ] = eulerC(3);
		numeric::geometry::hashing::Bin6D binC = binner.bin6( pC );
		//std::cout << "binC: ";
		//for ( Size ii = 1; ii <= 6; ++ii ) std::cout << binC[ ii ] << " "; std::cout << std::endl;
		boost::uint64_t indexC = binner.bin_index( pC );
		//std::cout << "indexC: " << indexC << std::endl;

		//binA: 4 12 4 15 4 0
		//indexA: 28656576
		//binB: 4 12 4 15 4 0
		//indexB: 28656576
		//binC: 4 12 4 33 22 0
		//indexC: 28668564

		TS_ASSERT( binA[ 1 ] == 4 );  TS_ASSERT( binA[ 2 ] == 12 ); TS_ASSERT( binA[ 3 ] == 4 );
		TS_ASSERT( binA[ 4 ] == 15 ); TS_ASSERT( binA[ 5 ] == 4 );  TS_ASSERT( binA[ 6 ] == 0 );
		TS_ASSERT( indexA == 28656576 );

		TS_ASSERT( binB[ 1 ] == 4 );  TS_ASSERT( binB[ 2 ] == 12 ); TS_ASSERT( binB[ 3 ] == 4 );
		TS_ASSERT( binB[ 4 ] == 15 ); TS_ASSERT( binB[ 5 ] == 4 );  TS_ASSERT( binB[ 6 ] == 0 );
		TS_ASSERT( indexB == 28656576 );

		TS_ASSERT( binC[ 1 ] == 4 );  TS_ASSERT( binC[ 2 ] == 12 ); TS_ASSERT( binC[ 3 ] == 4 );
		TS_ASSERT( binC[ 4 ] == 33 ); TS_ASSERT( binC[ 5 ] == 22 ); TS_ASSERT( binC[ 6 ] == 0 );
		TS_ASSERT( indexC == 28668564 );


		//boost::uint64_t indexA = binner.bin_index( pA );
		//boost::uint64_t indexB = binner.bin_index( pB );
		//boost::uint64_t indexC = binner.bin_index( pC );

		TS_ASSERT( indexA == indexB );
		TS_ASSERT( indexA != indexC );
		TS_ASSERT( indexB != indexC );


		Size3 euler_offsets2( 0 );
		euler_offsets2[ 3 ] = 1;

		numeric::geometry::hashing::SixDCoordinateBinner binner2( bb, euler_offsets2, binwidths );
		//std::cout << "binner2 dimsizes:";
		//for ( Size ii = 1; ii <= 6; ++ii ) std::cout << " " << binner2.dimsizes()[ ii ]; std::cout << std::endl;

		boost::uint64_t indexA2 = binner2.bin_index( pA );
		boost::uint64_t indexB2 = binner2.bin_index( pB );
		boost::uint64_t indexC2 = binner2.bin_index( pC );

		TS_ASSERT( indexA2 != indexB2 );
		TS_ASSERT( indexA2 != indexC2 );
		TS_ASSERT( indexB2 == indexC2 );

		numeric::geometry::hashing::Bin6D binA2 = binner2.bin6( pA );
		//std::cout << "binA2: ";
		//for ( Size ii = 1; ii <= 6; ++ii ) std::cout << binA2[ ii ] << " "; std::cout << std::endl;
		//std::cout << "indexA2: " << indexA2 << std::endl;

		numeric::geometry::hashing::Bin6D binB2 = binner2.bin6( pB );
		//std::cout << "binB2: ";
		//for ( Size ii = 1; ii <= 6; ++ii ) std::cout << binB2[ ii ] << " "; std::cout << std::endl;
		//std::cout << "indexB2: " << indexB2 << std::endl;

		numeric::geometry::hashing::Bin6D binC2 = binner2.bin6( pC );
		//std::cout << "binC2: ";
		//for ( Size ii = 1; ii <= 6; ++ii ) std::cout << binC2[ ii ] << " "; std::cout << std::endl;
		//std::cout << "indexC2: " << indexC2 << std::endl;

		//binA2: 4 12 4 15 4 1
		//indexA2: 30248609
		//binB2: 4 12 4 15 4 0
		//indexB2: 30248608
		//binC2: 4 12 4 15 4 0
		//indexC2: 30248608

		TS_ASSERT( binA2[ 1 ] == 4 );  TS_ASSERT( binA2[ 2 ] == 12 ); TS_ASSERT( binA2[ 3 ] == 4 );
		TS_ASSERT( binA2[ 4 ] == 15 ); TS_ASSERT( binA2[ 5 ] == 4 );  TS_ASSERT( binA2[ 6 ] == 1 );
		TS_ASSERT( indexA2 == 30248609 );

		TS_ASSERT( binB2[ 1 ] == 4 );  TS_ASSERT( binB2[ 2 ] == 12 ); TS_ASSERT( binB2[ 3 ] == 4 );
		TS_ASSERT( binB2[ 4 ] == 15 ); TS_ASSERT( binB2[ 5 ] == 4 );  TS_ASSERT( binB2[ 6 ] == 0 );
		TS_ASSERT( indexB2 == 30248608 );

		TS_ASSERT( binC2[ 1 ] == 4 );  TS_ASSERT( binC2[ 2 ] == 12 ); TS_ASSERT( binC2[ 3 ] == 4 );
		TS_ASSERT( binC2[ 4 ] == 15 ); TS_ASSERT( binC2[ 5 ] == 4 );  TS_ASSERT( binC2[ 6 ] == 0 );
		TS_ASSERT( indexC2 == 30248608 );

	}


	void test_SixDCoordinateBinner_wrap_theta_and_phi() {
		using namespace protocols::match;

		Vector lower( 12.5, 16.25, 4.25 );
		Vector upper( 15.5, 20, 8.5 );
		BoundingBox bb( lower, upper );
		Real6 binwidths;
		binwidths[ 1 ] = binwidths[ 2 ] = binwidths[ 3 ] = 0.25;
		binwidths[ 4 ] = binwidths[ 5 ] = binwidths[ 6 ] = 10;

		Size3 euler_offsets( 0 );

		numeric::geometry::hashing::SixDCoordinateBinner binner( bb, euler_offsets, binwidths );

		Real6 pAa, pBa, pCa;
		Real6 pAb, pBb, pCb;
		pAa[ 1 ] = 13.6; pAa[ 2 ] = 19.4; pAa[ 3 ] = 5.3;
		pBa = pCa = pAb = pBb = pCb = pAa;


		HTReal zrot1a, zrot1b, xrot2A, xrot2B, xrot2C, zrot3;

		zrot1a.set_zaxis_rotation_deg( 155 );
		zrot1b.set_zaxis_rotation_deg( 162 );
		xrot2A.set_xaxis_rotation_deg( 173 );
		xrot2B.set_xaxis_rotation_deg( 177 );
		xrot2C.set_xaxis_rotation_deg( 182 );
		zrot3.set_zaxis_rotation_deg( 44 );

		HTReal tAa = zrot1a * xrot2A * zrot3;
		HTReal tBa = zrot1a * xrot2B * zrot3;
		HTReal tCa = zrot1a * xrot2C * zrot3;
		HTReal tAb = zrot1b * xrot2A * zrot3;
		HTReal tBb = zrot1b * xrot2B * zrot3;
		HTReal tCb = zrot1b * xrot2C * zrot3;

		Vector eulerAa = tAa.euler_angles_deg(), eulerBa = tBa.euler_angles_deg(), eulerCa = tCa.euler_angles_deg();
		Vector eulerAb = tAb.euler_angles_deg(), eulerBb = tBb.euler_angles_deg(), eulerCb = tCb.euler_angles_deg();

		for ( Size ii = 1; ii <= 2; ++ii ) if ( eulerAa( ii ) < 0 ) eulerAa( ii ) += 360;
		for ( Size ii = 1; ii <= 2; ++ii ) if ( eulerBa( ii ) < 0 ) eulerBa( ii ) += 360;
		for ( Size ii = 1; ii <= 2; ++ii ) if ( eulerCa( ii ) < 0 ) eulerCa( ii ) += 360;

		for ( Size ii = 1; ii <= 2; ++ii ) if ( eulerAb( ii ) < 0 ) eulerAb( ii ) += 360;
		for ( Size ii = 1; ii <= 2; ++ii ) if ( eulerBb( ii ) < 0 ) eulerBb( ii ) += 360;
		for ( Size ii = 1; ii <= 2; ++ii ) if ( eulerCb( ii ) < 0 ) eulerCb( ii ) += 360;

		//for ( Size ii = 1; ii <= 3; ++ii ) std::cout << " " << eulerAa( ii ); std::cout << std::endl;
		//for ( Size ii = 1; ii <= 3; ++ii ) std::cout << " " << eulerBa( ii ); std::cout << std::endl;
		//for ( Size ii = 1; ii <= 3; ++ii ) std::cout << " " << eulerCa( ii ); std::cout << std::endl;

		//for ( Size ii = 1; ii <= 3; ++ii ) std::cout << " " << eulerAb( ii ); std::cout << std::endl;
		//for ( Size ii = 1; ii <= 3; ++ii ) std::cout << " " << eulerBb( ii ); std::cout << std::endl;
		//for ( Size ii = 1; ii <= 3; ++ii ) std::cout << " " << eulerCb( ii ); std::cout << std::endl;

		pAa[ 4 ] = eulerAa(1);   pAa[ 5 ] = eulerAa(2);  pAa[ 6 ] = eulerAa(3);
		//numeric::geometry::hashing::Bin6D binAa = binner.bin6( pAa );
		//std::cout << "binAa: ";
		//for ( Size ii = 1; ii <= 6; ++ii ) std::cout << binAa[ ii ] << " "; std::cout << std::endl;
		boost::uint64_t indexAa = binner.bin_index( pAa );
		//std::cout << "indexAa: " << indexAa << std::endl;

		pBa[ 4 ] = eulerBa(1);   pBa[ 5 ] = eulerBa(2);  pBa[ 6 ] = eulerBa(3);
		//numeric::geometry::hashing::Bin6D binBa = binner.bin6( pBa );
		//std::cout << "binBa: ";
		//for ( Size ii = 1; ii <= 6; ++ii ) std::cout << binBa[ ii ] << " "; std::cout << std::endl;
		boost::uint64_t indexBa = binner.bin_index( pBa );
		//std::cout << "indexBa: " << indexBa << std::endl;

		pCa[ 4 ] = eulerCa(1);   pCa[ 5 ] = eulerCa(2);  pCa[ 6 ] = eulerCa(3);
		//numeric::geometry::hashing::Bin6D binCa = binner.bin6( pCa );
		//std::cout << "binCa: ";
		//for ( Size ii = 1; ii <= 6; ++ii ) std::cout << binCa[ ii ] << " ";
		//std::cout << std::endl;
		boost::uint64_t indexCa = binner.bin_index( pCa );
		//std::cout << "indexCa: " << indexCa << std::endl;

		pAb[ 4 ] = eulerAb(1);   pAb[ 5 ] = eulerAb(2);  pAb[ 6 ] = eulerAb(3);
		//numeric::geometry::hashing::Bin6D binAb = binner.bin6( pAb );
		//std::cout << "binAb: ";
		//for ( Size ii = 1; ii <= 6; ++ii ) std::cout << binAb[ ii ] << " ";
		//std::cout << std::endl;
		boost::uint64_t indexAb = binner.bin_index( pAb );
		//std::cout << "indexAb: " << indexAb << std::endl;

		pBb[ 4 ] = eulerBb(1);   pBb[ 5 ] = eulerBb(2);  pBb[ 6 ] = eulerBb(3);
		//numeric::geometry::hashing::Bin6D binBb = binner.bin6( pBb );
		//std::cout << "binBb: ";
		//for ( Size ii = 1; ii <= 6; ++ii ) std::cout << binBb[ ii ] << " ";
		//std::cout << std::endl;
		boost::uint64_t indexBb = binner.bin_index( pBb );
		//std::cout << "indexBb: " << indexBb << std::endl;

		pCb[ 4 ] = eulerCb(1);   pCb[ 5 ] = eulerCb(2);  pCb[ 6 ] = eulerCb(3);
		//numeric::geometry::hashing::Bin6D binCb = binner.bin6( pCb );
		//std::cout << "binCb: ";
		//for ( Size ii = 1; ii <= 6; ++ii ) std::cout << binCb[ ii ] << " ";
		//std::cout << std::endl;
		boost::uint64_t indexCb = binner.bin_index( pCb );
		//std::cout << "indexCb: " << indexCb << std::endl;


		//boost::uint64_t indexA = binner.bin_index( pA );
		//boost::uint64_t indexB = binner.bin_index( pB );
		//boost::uint64_t indexC = binner.bin_index( pC );

		TS_ASSERT( indexAa == indexBa );
		TS_ASSERT( indexAa != indexCa );
		TS_ASSERT( indexAa != indexAb );
		TS_ASSERT( indexAa != indexBb );
		TS_ASSERT( indexAa != indexCb );
		TS_ASSERT( indexBa != indexCa );
		TS_ASSERT( indexBa != indexAb );
		TS_ASSERT( indexBa != indexBb );
		TS_ASSERT( indexBa != indexCb );
		TS_ASSERT( indexCa != indexAb );
		TS_ASSERT( indexCa != indexBb );
		TS_ASSERT( indexCa != indexCb );
		TS_ASSERT( indexAb == indexBb );
		TS_ASSERT( indexAb != indexCb );
		TS_ASSERT( indexBb != indexCb );

		// wrap theta only
		Size3 euler_offsets2( 0 );
		euler_offsets2[ 3 ] = 1;

		numeric::geometry::hashing::SixDCoordinateBinner binner2( bb, euler_offsets2, binwidths );
		//std::cout << "binner2 dimsizes:";
		//for ( Size ii = 1; ii <= 6; ++ii ) std::cout << " " << binner2.dimsizes()[ ii ]; std::cout << std::endl;

		boost::uint64_t indexAa2 = binner2.bin_index( pAa );
		boost::uint64_t indexBa2 = binner2.bin_index( pBa );
		boost::uint64_t indexCa2 = binner2.bin_index( pCa );
		boost::uint64_t indexAb2 = binner2.bin_index( pAb );
		boost::uint64_t indexBb2 = binner2.bin_index( pBb );
		boost::uint64_t indexCb2 = binner2.bin_index( pCb );

		TS_ASSERT( indexAa2 != indexBa2 );
		TS_ASSERT( indexAa2 != indexCa2 );
		TS_ASSERT( indexAa2 != indexAb2 );
		TS_ASSERT( indexAa2 != indexBb2 );
		TS_ASSERT( indexAa2 != indexCb2 );
		TS_ASSERT( indexBa2 == indexCa2 );
		TS_ASSERT( indexBa2 != indexAb2 );
		TS_ASSERT( indexBa2 != indexBb2 );
		TS_ASSERT( indexBa2 != indexCb2 );
		TS_ASSERT( indexCa2 != indexAb2 );
		TS_ASSERT( indexCa2 != indexBb2 );
		TS_ASSERT( indexCa2 != indexCb2 );
		TS_ASSERT( indexAb2 != indexBb2 );
		TS_ASSERT( indexAb2 != indexCb2 );
		TS_ASSERT( indexBb2 == indexCb2 );

		// wrap phi only
		Size3 euler_offsets3( 0 );
		euler_offsets3[ 1 ] = 1;

		numeric::geometry::hashing::SixDCoordinateBinner binner3( bb, euler_offsets3, binwidths );
		//std::cout << "binner3 dimsizes:";
		//for ( Size ii = 1; ii <= 6; ++ii ) std::cout << " " << binner3.dimsizes()[ ii ]; std::cout << std::endl;

		boost::uint64_t indexAa3 = binner3.bin_index( pAa );
		boost::uint64_t indexBa3 = binner3.bin_index( pBa );
		boost::uint64_t indexCa3 = binner3.bin_index( pCa );
		boost::uint64_t indexAb3 = binner3.bin_index( pAb );
		boost::uint64_t indexBb3 = binner3.bin_index( pBb );
		boost::uint64_t indexCb3 = binner3.bin_index( pCb );

		TS_ASSERT( indexAa3 == indexBa3 );
		TS_ASSERT( indexAa3 != indexCa3 );
		TS_ASSERT( indexAa3 == indexAb3 );
		TS_ASSERT( indexAa3 == indexBb3 );
		TS_ASSERT( indexAa3 != indexCb3 );
		TS_ASSERT( indexBa3 != indexCa3 );
		TS_ASSERT( indexBa3 == indexAb3 );
		TS_ASSERT( indexBa3 == indexBb3 );
		TS_ASSERT( indexBa3 != indexCb3 );
		TS_ASSERT( indexCa3 != indexAb3 );
		TS_ASSERT( indexCa3 != indexBb3 );
		TS_ASSERT( indexCa3 == indexCb3 );
		TS_ASSERT( indexAb3 == indexBb3 );
		TS_ASSERT( indexAb3 != indexCb3 );
		TS_ASSERT( indexBb3 != indexCb3 );

		// wrap both theta and phi
		Size3 euler_offsets4( 0 );
		euler_offsets4[ 1 ] = 1;
		euler_offsets4[ 3 ] = 1;

		numeric::geometry::hashing::SixDCoordinateBinner binner4( bb, euler_offsets4, binwidths );
		//std::cout << "binner4 dimsizes:";
		//for ( Size ii = 1; ii <= 6; ++ii ) std::cout << " " << binner4.dimsizes()[ ii ]; std::cout << std::endl;

		boost::uint64_t indexAa4 = binner4.bin_index( pAa );
		boost::uint64_t indexBa4 = binner4.bin_index( pBa );
		boost::uint64_t indexCa4 = binner4.bin_index( pCa );
		boost::uint64_t indexAb4 = binner4.bin_index( pAb );
		boost::uint64_t indexBb4 = binner4.bin_index( pBb );
		boost::uint64_t indexCb4 = binner4.bin_index( pCb );

		TS_ASSERT( indexAa4 != indexBa4 );
		TS_ASSERT( indexAa4 != indexCa4 );
		TS_ASSERT( indexAa4 == indexAb4 );
		TS_ASSERT( indexAa4 != indexBb4 );
		TS_ASSERT( indexAa4 != indexCb4 );
		TS_ASSERT( indexBa4 == indexCa4 );
		TS_ASSERT( indexBa4 != indexAb4 );
		TS_ASSERT( indexBa4 == indexBb4 );
		TS_ASSERT( indexBa4 == indexCb4 );
		TS_ASSERT( indexCa4 != indexAb4 );
		TS_ASSERT( indexCa4 == indexBb4 );
		TS_ASSERT( indexCa4 == indexCb4 );
		TS_ASSERT( indexAb4 != indexBb4 );
		TS_ASSERT( indexAb4 != indexCb4 );
		TS_ASSERT( indexBb4 == indexCb4 );
	}


	void test_SixDCoordinateBinner_halfbin_index_no_offset() {
		using namespace protocols::match;

		Vector lower( 12.5, 16.25, 4.25 );
		Vector upper( 15.5, 20, 8.5 );
		BoundingBox bb( lower, upper );
		Real6 binwidths;
		binwidths[ 1 ] = binwidths[ 2 ] = binwidths[ 3 ] = 0.25;
		binwidths[ 4 ] = binwidths[ 5 ] = binwidths[ 6 ] = 10;

		Size3 euler_offsets( 0 );

		numeric::geometry::hashing::SixDCoordinateBinner binner( bb, euler_offsets, binwidths );

		Real6 pA, pB, pC, pD, pE;
		pA[ 1 ] = 13.6; pA[ 2 ] = 19.4; pA[ 3 ] = 5.3;
		pA[ 4 ] = 49;   pA[ 5 ] = 127;  pA[ 6 ] = 176;

		pB[ 1 ] = 13.7; pB[ 2 ] = 19.3; pB[ 3 ] = 5.2;
		pB[ 4 ] = 50;   pB[ 5 ] = 122;  pB[ 6 ] = 176;

		pC[ 1 ] = 13.6; pC[ 2 ] = 19.4; pC[ 3 ] = 5.3;
		pC[ 4 ] = 57;   pC[ 5 ] = 127;  pC[ 6 ] = 171;

		pD[ 1 ] = 13.6; pD[ 2 ] = 19.4; pD[ 3 ] = 5.3;
		pD[ 4 ] = 249;  pD[ 5 ] = 127;  pD[ 6 ] = 176;

		pE[ 1 ] = 13.7; pE[ 2 ] = 19.3; pE[ 3 ] = 5.2;
		pE[ 4 ] = 256;  pE[ 5 ] = 122;  pE[ 6 ] = 176;


		numeric::geometry::hashing::Bin6D hbA = binner.halfbin6( pA );
		numeric::geometry::hashing::Bin6D hbB = binner.halfbin6( pB );
		numeric::geometry::hashing::Bin6D hbC = binner.halfbin6( pC );
		numeric::geometry::hashing::Bin6D hbD = binner.halfbin6( pD );
		numeric::geometry::hashing::Bin6D hbE = binner.halfbin6( pE );

		TS_ASSERT( hbA[ 1 ] == 0 ); TS_ASSERT( hbA[ 2 ] == 1 ); TS_ASSERT( hbA[ 3 ] == 0 );
		TS_ASSERT( hbA[ 4 ] == 1 ); TS_ASSERT( hbA[ 5 ] == 1 ); TS_ASSERT( hbA[ 6 ] == 1 );

		TS_ASSERT( hbB[ 1 ] == 1 ); TS_ASSERT( hbB[ 2 ] == 0 ); TS_ASSERT( hbB[ 3 ] == 1 );
		TS_ASSERT( hbB[ 4 ] == 0 ); TS_ASSERT( hbB[ 5 ] == 0 ); TS_ASSERT( hbB[ 6 ] == 1 );

		TS_ASSERT( hbC[ 1 ] == 0 ); TS_ASSERT( hbC[ 2 ] == 1 ); TS_ASSERT( hbC[ 3 ] == 0 );
		TS_ASSERT( hbC[ 4 ] == 1 ); TS_ASSERT( hbC[ 5 ] == 1 ); TS_ASSERT( hbC[ 6 ] == 0 );

		TS_ASSERT( hbD[ 1 ] == 0 ); TS_ASSERT( hbD[ 2 ] == 1 ); TS_ASSERT( hbD[ 3 ] == 0 );
		TS_ASSERT( hbD[ 4 ] == 1 ); TS_ASSERT( hbD[ 5 ] == 1 ); TS_ASSERT( hbD[ 6 ] == 1 );

		TS_ASSERT( hbE[ 1 ] == 1 ); TS_ASSERT( hbE[ 2 ] == 0 ); TS_ASSERT( hbE[ 3 ] == 1 );
		TS_ASSERT( hbE[ 4 ] == 1 ); TS_ASSERT( hbE[ 5 ] == 0 ); TS_ASSERT( hbE[ 6 ] == 1 );

	}

	void test_SixDCoordinateBinner_halfbin_index_euclid_offset() {
		using namespace protocols::match;

		Vector lower( 12.5, 16.25, 4.25 );
		lower += Vector( 0.125 );
		Vector upper( 15.5, 20, 8.5 );
		BoundingBox bb( lower, upper );
		Real6 binwidths;
		binwidths[ 1 ] = binwidths[ 2 ] = binwidths[ 3 ] = 0.25;
		binwidths[ 4 ] = binwidths[ 5 ] = binwidths[ 6 ] = 10;

		Size3 euler_offsets( 0 );

		numeric::geometry::hashing::SixDCoordinateBinner binner( bb, euler_offsets, binwidths );

		Real6 pA, pB, pC;
		pA[ 1 ] = 13.6; pA[ 2 ] = 19.4; pA[ 3 ] = 5.3;
		pA[ 4 ] = 49;   pA[ 5 ] = 127;  pA[ 6 ] = 176;

		pB[ 1 ] = 13.7; pB[ 2 ] = 19.3; pB[ 3 ] = 5.2;
		pB[ 4 ] = 50;   pB[ 5 ] = 122;  pB[ 6 ] = 176;

		pC[ 1 ] = 13.6; pC[ 2 ] = 19.4; pC[ 3 ] = 5.3;
		pC[ 4 ] = 57;   pC[ 5 ] = 127;  pC[ 6 ] = 171;

		numeric::geometry::hashing::Bin6D hbA = binner.halfbin6( pA );
		numeric::geometry::hashing::Bin6D hbB = binner.halfbin6( pB );
		numeric::geometry::hashing::Bin6D hbC = binner.halfbin6( pC );

		TS_ASSERT( hbA[ 1 ] == 1 ); TS_ASSERT( hbA[ 2 ] == 0 ); TS_ASSERT( hbA[ 3 ] == 1 );
		TS_ASSERT( hbA[ 4 ] == 1 ); TS_ASSERT( hbA[ 5 ] == 1 ); TS_ASSERT( hbA[ 6 ] == 1 );

		TS_ASSERT( hbB[ 1 ] == 0 ); TS_ASSERT( hbB[ 2 ] == 1 ); TS_ASSERT( hbB[ 3 ] == 0 );
		TS_ASSERT( hbB[ 4 ] == 0 ); TS_ASSERT( hbB[ 5 ] == 0 ); TS_ASSERT( hbB[ 6 ] == 1 );

		TS_ASSERT( hbC[ 1 ] == 1 ); TS_ASSERT( hbC[ 2 ] == 0 ); TS_ASSERT( hbC[ 3 ] == 1 );
		TS_ASSERT( hbC[ 4 ] == 1 ); TS_ASSERT( hbC[ 5 ] == 1 ); TS_ASSERT( hbC[ 6 ] == 0 );

	}

	void test_SixDCoordinateBinner_halfbin_index_euler_phi_psi_offset() {
		using namespace protocols::match;

		Vector lower( 12.5, 16.25, 4.25 );
		Vector upper( 15.5, 20, 8.5 );
		BoundingBox bb( lower, upper );
		Real6 binwidths;
		binwidths[ 1 ] = binwidths[ 2 ] = binwidths[ 3 ] = 0.25;
		binwidths[ 4 ] = binwidths[ 5 ] = binwidths[ 6 ] = 10;

		Size3 euler_offsets( 1 );
		euler_offsets[ 3 ] = 0;

		numeric::geometry::hashing::SixDCoordinateBinner binner( bb, euler_offsets, binwidths );

		Real6 pA, pB, pC, pD, pE;
		pA[ 1 ] = 13.6; pA[ 2 ] = 19.4; pA[ 3 ] = 5.3;
		pA[ 4 ] = 49;   pA[ 5 ] = 127;  pA[ 6 ] = 176;

		pB[ 1 ] = 13.7; pB[ 2 ] = 19.3; pB[ 3 ] = 5.2;
		pB[ 4 ] = 50.1;   pB[ 5 ] = 122;  pB[ 6 ] = 176;

		pC[ 1 ] = 13.6; pC[ 2 ] = 19.4; pC[ 3 ] = 5.3;
		pC[ 4 ] = 57;   pC[ 5 ] = 127;  pC[ 6 ] = 171;

		pD[ 1 ] = 13.6; pD[ 2 ] = 19.4; pD[ 3 ] = 5.3;
		pD[ 4 ] = 249;  pD[ 5 ] = 127;  pD[ 6 ] = 176;

		pE[ 1 ] = 13.7; pE[ 2 ] = 19.3; pE[ 3 ] = 5.2;
		pE[ 4 ] = 256;  pE[ 5 ] = 122;  pE[ 6 ] = 176;

		numeric::geometry::hashing::Bin6D hbA = binner.halfbin6( pA );
		numeric::geometry::hashing::Bin6D hbB = binner.halfbin6( pB );
		numeric::geometry::hashing::Bin6D hbC = binner.halfbin6( pC );
		numeric::geometry::hashing::Bin6D hbD = binner.halfbin6( pD );
		numeric::geometry::hashing::Bin6D hbE = binner.halfbin6( pE );

		TS_ASSERT( hbA[ 1 ] == 0 ); TS_ASSERT( hbA[ 2 ] == 1 ); TS_ASSERT( hbA[ 3 ] == 0 );
		TS_ASSERT( hbA[ 4 ] == 0 ); TS_ASSERT( hbA[ 5 ] == 0 ); TS_ASSERT( hbA[ 6 ] == 1 );

		TS_ASSERT( hbB[ 1 ] == 1 ); TS_ASSERT( hbB[ 2 ] == 0 ); TS_ASSERT( hbB[ 3 ] == 1 );
		TS_ASSERT( hbB[ 4 ] == 1 ); TS_ASSERT( hbB[ 5 ] == 1 ); TS_ASSERT( hbB[ 6 ] == 1 );

		TS_ASSERT( hbC[ 1 ] == 0 ); TS_ASSERT( hbC[ 2 ] == 1 ); TS_ASSERT( hbC[ 3 ] == 0 );
		TS_ASSERT( hbC[ 4 ] == 0 ); TS_ASSERT( hbC[ 5 ] == 0 ); TS_ASSERT( hbC[ 6 ] == 0 );

		TS_ASSERT( hbD[ 1 ] == 0 ); TS_ASSERT( hbD[ 2 ] == 1 ); TS_ASSERT( hbD[ 3 ] == 0 );
		TS_ASSERT( hbD[ 4 ] == 0 ); TS_ASSERT( hbD[ 5 ] == 0 ); TS_ASSERT( hbD[ 6 ] == 1 );

		TS_ASSERT( hbE[ 1 ] == 1 ); TS_ASSERT( hbE[ 2 ] == 0 ); TS_ASSERT( hbE[ 3 ] == 1 );
		TS_ASSERT( hbE[ 4 ] == 0 ); TS_ASSERT( hbE[ 5 ] == 1 ); TS_ASSERT( hbE[ 6 ] == 1 );
	}

	void test_SixDCoordinateBinner_halfbin_index_theta_offset_euclid_offset() {
		using namespace protocols::match;

		Vector lower( 12.5, 16.25, 4.25 );
		lower += Vector( 0.125 );
		Vector upper( 15.5, 20, 8.5 );
		BoundingBox bb( lower, upper );
		Real6 binwidths;
		binwidths[ 1 ] = binwidths[ 2 ] = binwidths[ 3 ] = 0.25;
		binwidths[ 4 ] = binwidths[ 5 ] = binwidths[ 6 ] = 10;

		Size3 euler_offsets( 0 ); euler_offsets[ 3 ] = 1;

		numeric::geometry::hashing::SixDCoordinateBinner binner( bb, euler_offsets, binwidths );

		Real6 pA, pB, pC, pD, pE;
		pA[ 1 ] = 13.6; pA[ 2 ] = 19.4; pA[ 3 ] = 5.3;
		pA[ 4 ] = 49;   pA[ 5 ] = 127;  pA[ 6 ] = 176;

		pB[ 1 ] = 13.7; pB[ 2 ] = 19.3; pB[ 3 ] = 5.2;
		pB[ 4 ] = 50;   pB[ 5 ] = 122;  pB[ 6 ] = 176;

		pC[ 1 ] = 13.6; pC[ 2 ] = 19.4; pC[ 3 ] = 5.3;
		pC[ 4 ] = 57;   pC[ 5 ] = 127;  pC[ 6 ] = 171;

		pD[ 1 ] = 13.6; pD[ 2 ] = 19.4; pD[ 3 ] = 5.3;
		pD[ 4 ] = 249;  pD[ 5 ] = 127;  pD[ 6 ] = 176;

		pE[ 1 ] = 13.7; pE[ 2 ] = 19.3; pE[ 3 ] = 5.2;
		pE[ 4 ] = 256;  pE[ 5 ] = 122;  pE[ 6 ] = 176;

		numeric::geometry::hashing::Bin6D hbA = binner.halfbin6( pA );
		numeric::geometry::hashing::Bin6D hbB = binner.halfbin6( pB );
		numeric::geometry::hashing::Bin6D hbC = binner.halfbin6( pC );
		numeric::geometry::hashing::Bin6D hbD = binner.halfbin6( pD );
		numeric::geometry::hashing::Bin6D hbE = binner.halfbin6( pE );

		TS_ASSERT( hbA[ 1 ] == 1 ); TS_ASSERT( hbA[ 2 ] == 0 ); TS_ASSERT( hbA[ 3 ] == 1 );
		TS_ASSERT( hbA[ 4 ] == 1 ); TS_ASSERT( hbA[ 5 ] == 1 ); TS_ASSERT( hbA[ 6 ] == 0 );

		TS_ASSERT( hbB[ 1 ] == 0 ); TS_ASSERT( hbB[ 2 ] == 1 ); TS_ASSERT( hbB[ 3 ] == 0 );
		TS_ASSERT( hbB[ 4 ] == 0 ); TS_ASSERT( hbB[ 5 ] == 0 ); TS_ASSERT( hbB[ 6 ] == 0 );

		TS_ASSERT( hbC[ 1 ] == 1 ); TS_ASSERT( hbC[ 2 ] == 0 ); TS_ASSERT( hbC[ 3 ] == 1 );
		TS_ASSERT( hbC[ 4 ] == 1 ); TS_ASSERT( hbC[ 5 ] == 1 ); TS_ASSERT( hbC[ 6 ] == 1 );

		TS_ASSERT( hbD[ 1 ] == 1 ); TS_ASSERT( hbD[ 2 ] == 0 ); TS_ASSERT( hbD[ 3 ] == 1 );
		TS_ASSERT( hbD[ 4 ] == 1 ); TS_ASSERT( hbD[ 5 ] == 1 ); TS_ASSERT( hbD[ 6 ] == 0 );

		TS_ASSERT( hbE[ 1 ] == 0 ); TS_ASSERT( hbE[ 2 ] == 1 ); TS_ASSERT( hbE[ 3 ] == 0 );
		TS_ASSERT( hbE[ 4 ] == 1 ); TS_ASSERT( hbE[ 5 ] == 0 ); TS_ASSERT( hbE[ 6 ] == 0 );

	}

	void test_SixDCoordinateBinner_halfbin_index_euler_offset() {
		using namespace protocols::match;

		Vector lower( 12.5, 16.25, 4.25 );
		Vector upper( 15.5, 20, 8.5 );
		BoundingBox bb( lower, upper );
		Real6 binwidths;
		binwidths[ 1 ] = binwidths[ 2 ] = binwidths[ 3 ] = 0.25;
		binwidths[ 4 ] = binwidths[ 5 ] = binwidths[ 6 ] = 10;

		Size3 euler_offsets( 1 );

		numeric::geometry::hashing::SixDCoordinateBinner binner( bb, euler_offsets, binwidths );

		Real6 pA, pB, pC, pD, pE;
		pA[ 1 ] = 13.6; pA[ 2 ] = 19.4; pA[ 3 ] = 5.3;
		pA[ 4 ] = 49;   pA[ 5 ] = 127;  pA[ 6 ] = 176;

		pB[ 1 ] = 13.7; pB[ 2 ] = 19.3; pB[ 3 ] = 5.2;
		pB[ 4 ] = 50.1;   pB[ 5 ] = 122;  pB[ 6 ] = 176;

		pC[ 1 ] = 13.6; pC[ 2 ] = 19.4; pC[ 3 ] = 5.3;
		pC[ 4 ] = 57;   pC[ 5 ] = 127;  pC[ 6 ] = 171;

		pD[ 1 ] = 13.6; pD[ 2 ] = 19.4; pD[ 3 ] = 5.3;
		pD[ 4 ] = 249;  pD[ 5 ] = 127;  pD[ 6 ] = 176;

		pE[ 1 ] = 13.7; pE[ 2 ] = 19.3; pE[ 3 ] = 5.2;
		pE[ 4 ] = 256;  pE[ 5 ] = 122;  pE[ 6 ] = 176;

		numeric::geometry::hashing::Bin6D hbA = binner.halfbin6( pA );
		numeric::geometry::hashing::Bin6D hbB = binner.halfbin6( pB );
		numeric::geometry::hashing::Bin6D hbC = binner.halfbin6( pC );
		numeric::geometry::hashing::Bin6D hbD = binner.halfbin6( pD );
		numeric::geometry::hashing::Bin6D hbE = binner.halfbin6( pE );

		TS_ASSERT( hbA[ 1 ] == 0 ); TS_ASSERT( hbA[ 2 ] == 1 ); TS_ASSERT( hbA[ 3 ] == 0 );
		TS_ASSERT( hbA[ 4 ] == 0 ); TS_ASSERT( hbA[ 5 ] == 0 ); TS_ASSERT( hbA[ 6 ] == 0 );

		TS_ASSERT( hbB[ 1 ] == 1 ); TS_ASSERT( hbB[ 2 ] == 0 ); TS_ASSERT( hbB[ 3 ] == 1 );
		TS_ASSERT( hbB[ 4 ] == 1 ); TS_ASSERT( hbB[ 5 ] == 1 ); TS_ASSERT( hbB[ 6 ] == 0 );

		TS_ASSERT( hbC[ 1 ] == 0 ); TS_ASSERT( hbC[ 2 ] == 1 ); TS_ASSERT( hbC[ 3 ] == 0 );
		TS_ASSERT( hbC[ 4 ] == 0 ); TS_ASSERT( hbC[ 5 ] == 0 ); TS_ASSERT( hbC[ 6 ] == 1 );

		TS_ASSERT( hbD[ 1 ] == 0 ); TS_ASSERT( hbD[ 2 ] == 1 ); TS_ASSERT( hbD[ 3 ] == 0 );
		TS_ASSERT( hbD[ 4 ] == 0 ); TS_ASSERT( hbD[ 5 ] == 0 ); TS_ASSERT( hbD[ 6 ] == 0 );

		TS_ASSERT( hbE[ 1 ] == 1 ); TS_ASSERT( hbE[ 2 ] == 0 ); TS_ASSERT( hbE[ 3 ] == 1 );
		TS_ASSERT( hbE[ 4 ] == 0 ); TS_ASSERT( hbE[ 5 ] == 1 ); TS_ASSERT( hbE[ 6 ] == 0 );

	}

	void test_HitNeighborFinder_find_connected_components_1()
	{
		using namespace protocols::match;
		HitNeighborFinder finder;

		Vector lower( 12.5, 16.25, 4.25 );
		Vector upper( 15.5, 20, 8.5 );

		BoundingBox bb( lower, upper );
		finder.set_bounding_box( bb );

		Real6 binwidths;
		binwidths[ 1 ] = binwidths[ 2 ] = binwidths[ 3 ] = 0.25; finder.set_uniform_xyz_bin_width( 1 );
		binwidths[ 4 ] = binwidths[ 5 ] = binwidths[ 6 ] = 10;   finder.set_uniform_euler_angle_bin_width( 10 );

		finder.initialize();

		Real6 pA, pB;//, pC, pD, pE;
		pA[ 1 ] = 13.6; pA[ 2 ] = 19.4; pA[ 3 ] = 5.3;
		pA[ 4 ] = 49;   pA[ 5 ] = 127;  pA[ 6 ] = 176;

		pB[ 1 ] = 13.7; pB[ 2 ] = 19.3; pB[ 3 ] = 5.2;
		pB[ 4 ] = 50.1; pB[ 5 ] = 122;  pB[ 6 ] = 176;

		Hit hA, hB; hA.second() = pA; hB.second() = pB;
		std::list< Hit > hitlist;
		hitlist.push_back( hA ); hitlist.push_back( hB );
		finder.add_hits( hitlist );

		utility::vector1< std::list< Hit const * > > hit_ccs = finder.connected_components();
		TS_ASSERT( hit_ccs.size() == 1 ); // both hits should be in the same connected component.
	}

	void test_HitNeighborFinder_find_connected_components_2()
	{
		using namespace protocols::match;
		HitNeighborFinder finder;

		Vector lower( 12.5, 16.25, 4.25 );
		Vector upper( 15.5, 20, 8.5 );

		BoundingBox bb( lower, upper );
		finder.set_bounding_box( bb );

		finder.set_uniform_xyz_bin_width( 1 );
		finder.set_uniform_euler_angle_bin_width( 10 );

		finder.initialize();

		Real6 pA, pB, pC; // pD, pE;
		pA[ 1 ] = 13.6; pA[ 2 ] = 19.4; pA[ 3 ] = 5.3;
		pA[ 4 ] = 49;   pA[ 5 ] = 127;  pA[ 6 ] = 176;

		pB[ 1 ] = 13.7; pB[ 2 ] = 19.3; pB[ 3 ] = 5.2;
		pB[ 4 ] = 50.1;   pB[ 5 ] = 122;  pB[ 6 ] = 176;

		pC[ 1 ] = 13.6; pC[ 2 ] = 19.4; pC[ 3 ] = 5.3;
		pC[ 4 ] = 62;   pC[ 5 ] = 127;  pC[ 6 ] = 171;

		Hit hA, hB, hC; hA.second() = pA; hB.second() = pB; hC.second() = pC;
		std::list< Hit > hitlist;
		hitlist.push_back( hA ); hitlist.push_back( hB ); hitlist.push_back( hC );
		finder.add_hits( hitlist );

		utility::vector1< std::list< Hit const * > > hit_ccs = finder.connected_components();
		TS_ASSERT( hit_ccs.size() == 2 ); // This should find two connected components
	}

	void test_HitNeighborFinder_find_connected_components_3()
	{
		using namespace protocols::match;
		HitNeighborFinder finder;

		Vector lower( 12.5, 16.25, 4.25 );
		Vector upper( 15.5, 20, 8.5 );

		BoundingBox bb( lower, upper );
		finder.set_bounding_box( bb );

		finder.set_uniform_xyz_bin_width( 1 );
		finder.set_uniform_euler_angle_bin_width( 10 );

		finder.initialize();

		Real6 pA, pB, pC, pD; //pE;
		pA[ 1 ] = 13.6; pA[ 2 ] = 19.4; pA[ 3 ] = 5.3;
		pA[ 4 ] = 49;   pA[ 5 ] = 127;  pA[ 6 ] = 176;

		pB[ 1 ] = 13.7; pB[ 2 ] = 19.3; pB[ 3 ] = 5.2;
		pB[ 4 ] = 50.1;   pB[ 5 ] = 122;  pB[ 6 ] = 176;

		pC[ 1 ] = 13.6; pC[ 2 ] = 19.4; pC[ 3 ] = 5.3;
		pC[ 4 ] = 62;   pC[ 5 ] = 127;  pC[ 6 ] = 171;

		pD[ 1 ] = 13.6; pD[ 2 ] = 19.4; pD[ 3 ] = 5.3; // pD will span the pB-pC gap
		pD[ 4 ] = 55;   pD[ 5 ] = 127;  pD[ 6 ] = 171;

		Hit hA, hB, hC, hD; hA.second() = pA; hB.second() = pB; hC.second() = pC; hD.second() = pD;
		std::list< Hit > hitlist;
		hitlist.push_back( hA ); hitlist.push_back( hB ); hitlist.push_back( hC ); hitlist.push_back( hD );
		finder.add_hits( hitlist );

		utility::vector1< std::list< Hit const * > > hit_ccs = finder.connected_components();
		TS_ASSERT( hit_ccs.size() == 1 ); // This should find only one connected component, now
	}

	/// Test that wrapping theta works
	void test_HitNeighborFinder_find_connected_components_4()
	{
		using namespace protocols::match;
		HitNeighborFinder finder;

		Vector lower( 12.5, 16.25, 4.25 );
		Vector upper( 15.5, 20, 8.5 );

		BoundingBox bb( lower, upper );
		finder.set_bounding_box( bb );

		finder.set_uniform_xyz_bin_width( 1 );
		finder.set_uniform_euler_angle_bin_width( 10 );

		finder.initialize();

		Real6 pA, pB, pC, pD; //pE;
		pA[ 1 ] = 13.6; pA[ 2 ] = 19.4; pA[ 3 ] = 5.3;
		pA[ 4 ] = 49;   pA[ 5 ] = 127;  pA[ 6 ] = 176;

		pB[ 1 ] = 13.7; pB[ 2 ] = 19.3; pB[ 3 ] = 5.2;
		pB[ 4 ] = 230.1;pB[ 5 ] = 302;  pB[ 6 ] = 176; // this should wrap to be a neighbor of A


		Hit hA, hB; hA.second() = pA; hB.second() = pB;
		std::list< Hit > hitlist;
		hitlist.push_back( hA ); hitlist.push_back( hB );
		finder.add_hits( hitlist );

		utility::vector1< std::list< Hit const * > > hit_ccs = finder.connected_components();
		TS_ASSERT( hit_ccs.size() == 1 ); // This should find only one connected component, now
	}


	void test_MatchCounter_two_geomcsts_expect_one_match()
	{
		using namespace protocols::match;
		MatchCounter counter;

		Vector lower( 12.5, 16.25, 4.25 );
		Vector upper( 15.5, 20, 8.5 );

		BoundingBox bb( lower, upper );
		counter.set_n_geometric_constraints( 2 );
		counter.set_bounding_box( bb );

		counter.set_uniform_xyz_bin_width( 1 );
		counter.set_uniform_euler_angle_bin_width( 10 );

		counter.initialize();

		Real6 pA, pB; //pC, pD; pE;
		pA[ 1 ] = 13.6; pA[ 2 ] = 19.4; pA[ 3 ] = 5.3;
		pA[ 4 ] = 49;   pA[ 5 ] = 127;  pA[ 6 ] = 176;

		pB[ 1 ] = 13.7; pB[ 2 ] = 19.3; pB[ 3 ] = 5.2;
		pB[ 4 ] = 50.1;   pB[ 5 ] = 122;  pB[ 6 ] = 176;

		//pC[ 1 ] = 13.6; pC[ 2 ] = 19.4; pC[ 3 ] = 5.3;
		//pC[ 4 ] = 62;   pC[ 5 ] = 127;  pC[ 6 ] = 171;

		//pD[ 1 ] = 13.6; pD[ 2 ] = 19.4; pD[ 3 ] = 5.3; // pD will span the pB-pC gap
		//pD[ 4 ] = 55;   pD[ 5 ] = 127;  pD[ 6 ] = 171;

		Hit hA, hB/*, hC, hD*/; hA.second() = pA; hB.second() = pB;// hC.second() = pC; hD.second() = pD;
		std::list< Hit > hitlist_a, hitlist_b;
		hitlist_a.push_back( hA ); hitlist_b.push_back( hB );// hitlist.push_back( hC ); hitlist.push_back( hD );
		counter.add_hits( 1, hitlist_a );
		counter.add_hits( 2, hitlist_b );

		Size n_matches = counter.count_n_matches();
		TS_ASSERT( n_matches == 1 );
	}

	void test_MatchCounter_three_geomcsts_expect_one_match()
	{
		using namespace protocols::match;
		MatchCounter counter;

		Vector lower( 12, 16, 4 );
		Vector upper( 15, 20, 8 );

		BoundingBox bb( lower, upper );
		counter.set_n_geometric_constraints( 3 );
		counter.set_bounding_box( bb );

		counter.set_uniform_xyz_bin_width( 1 );
		counter.set_uniform_euler_angle_bin_width( 10 );

		counter.initialize();

		Real6 pA, pB, pC;// pD; pE;
		pA[ 1 ] = 13.6; pA[ 2 ] = 19.4; pA[ 3 ] = 5.3;
		pA[ 4 ] = 49;   pA[ 5 ] = 127;  pA[ 6 ] = 176;

		pB[ 1 ] = 13.7; pB[ 2 ] = 19.3; pB[ 3 ] = 5.2;
		pB[ 4 ] = 50.1;   pB[ 5 ] = 122;  pB[ 6 ] = 176;

		pC[ 1 ] = 13.6; pC[ 2 ] = 19.4; pC[ 3 ] = 5.3;
		pC[ 4 ] = 48;   pC[ 5 ] = 126;  pC[ 6 ] = 171;

		//pD[ 1 ] = 13.6; pD[ 2 ] = 19.4; pD[ 3 ] = 5.3; // pD will span the pB-pC gap
		//pD[ 4 ] = 55;   pD[ 5 ] = 127;  pD[ 6 ] = 171;

		Hit hA, hB, hC/*, hD*/; hA.second() = pA; hB.second() = pB; hC.second() = pC;// hD.second() = pD;
		std::list< Hit > hitlist_a, hitlist_b, hitlist_c;
		hitlist_a.push_back( hA ); hitlist_b.push_back( hB ); hitlist_c.push_back( hC );// hitlist.push_back( hD );
		counter.add_hits( 1, hitlist_a );
		counter.add_hits( 2, hitlist_b );
		counter.add_hits( 3, hitlist_c );

		Size n_matches = counter.count_n_matches();
		TS_ASSERT( n_matches == 1 );
	}

	void test_MatchCounter_three_geomcsts_expect_zero_matches()
	{
		using namespace protocols::match;
		MatchCounter counter;

		Vector lower( 12, 16, 4 );
		Vector upper( 15, 20, 8 );

		BoundingBox bb( lower, upper );
		counter.set_n_geometric_constraints( 3 );
		counter.set_bounding_box( bb );

		counter.set_uniform_xyz_bin_width( 1 );
		counter.set_uniform_euler_angle_bin_width( 10 );

		counter.initialize();

		Real6 pA, pB, pC; // pD; pE;
		pA[ 1 ] = 13.6; pA[ 2 ] = 19.4; pA[ 3 ] = 5.3;
		pA[ 4 ] = 49;   pA[ 5 ] = 127;  pA[ 6 ] = 176;

		pB[ 1 ] = 14.1; pB[ 2 ] = 19.3; pB[ 3 ] = 5.2;
		pB[ 4 ] = 52;   pB[ 5 ] = 122;  pB[ 6 ] = 176;

		pC[ 1 ] = 13.1; pC[ 2 ] = 19.4; pC[ 3 ] = 5.3; // pC is near pA but too far from pB.
		pC[ 4 ] = 48;   pC[ 5 ] = 126;  pC[ 6 ] = 171;

		//pD[ 1 ] = 13.6; pD[ 2 ] = 19.4; pD[ 3 ] = 5.3; // pD will span the pB-pC gap
		//pD[ 4 ] = 55;   pD[ 5 ] = 127;  pD[ 6 ] = 171;

		Hit hA, hB, hC/*, hD*/; hA.second() = pA; hB.second() = pB; hC.second() = pC;// hD.second() = pD;
		std::list< Hit > hitlist_a, hitlist_b, hitlist_c;
		hitlist_a.push_back( hA ); hitlist_b.push_back( hB ); hitlist_c.push_back( hC );// hitlist.push_back( hD );
		counter.add_hits( 1, hitlist_a );
		counter.add_hits( 2, hitlist_b );
		counter.add_hits( 3, hitlist_c );

		Size n_matches = counter.count_n_matches();
		TS_ASSERT( n_matches == 0 );
	}


};
