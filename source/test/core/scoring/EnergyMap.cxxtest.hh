// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ScoreTest.cxxtest.hh
/// @brief  unified scoring test.
/// @author Sergey Lyskov

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>


#include <core/types.hh>

// Unit headers
#include <core/scoring/EnergyMap.hh>

#include <test/UTracer.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <numeric/xyzVector.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>


using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.EnergyMap.cxxtest");

// using declarations
using namespace std;
using namespace core;
using namespace scoring;

///////////////////////////////////////////////////////////////////////////
/// @name EnergyMapTest
/// @brief: Test the functionality of the EnergyMap class
///////////////////////////////////////////////////////////////////////////
class EnergyMapTest : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	/// @brief fill an energy map with Reals st emap[ii] = sin[ii]
	void initialize_energymap_sin( EnergyMap & emap ) {
		for ( Size ii = 1; ii <= n_score_types; ++ii ) {
			emap[ (ScoreType) ii ] = sin( ii );
		}
	}

	/// @brief fill an energy map with Reals st emap[ii] = cos[ii]
	void initialize_energymap_cos( EnergyMap & emap ) {
		for ( Size ii = 1; ii <= n_score_types; ++ii ) {
			emap[ (ScoreType) ii ] = cos( ii );
		}
	}


	/// @brief test that the energy map initializes all of its energies to zero
	void test_EnergyMap_Constructor() {
		EnergyMap emap;
		for ( Size ii = 1; ii <= n_score_types; ++ii ) {
			TS_ASSERT_EQUALS( emap[ (ScoreType) ii ], 0.0 );
		}
	}

	/// @brief test that the energy map copies all of its data in its assignment operator
	void test_EnergyMap_AssignmentOperator() {
		EnergyMap emap;
		initialize_energymap_sin( emap );
		EnergyMap copy_emap;
		copy_emap = emap;

		for ( Size ii = 1; ii <= n_score_types; ++ii ) {
			TS_ASSERT_EQUALS( emap[ (ScoreType) ii ], copy_emap[ (ScoreType) ii ] );
		}
	}

	/// @brief test that the non-constant begin iterator returned is the first element
	/// in the energy map, and that the dereferenced iterator may be assigned
	void test_EnergyMap_begin() {
		EnergyMap emap;
		initialize_energymap_sin( emap );
		EnergyMap::iterator emap_iter = emap.begin();
		TS_ASSERT_EQUALS( *emap_iter, sin( 1 ) );
		*emap_iter = sin( 2 );
		TS_ASSERT_DIFFERS( emap[ (ScoreType) 1 ], sin(1) );
		TS_ASSERT_EQUALS(  emap[ (ScoreType) 1 ], sin(2) );
	}

	/// @brief test that the constant begin iterator returned is the first element
	/// in the energy map
	void test_EnergyMap_constbegin() {
		EnergyMap emap;
		initialize_energymap_sin( emap );
		EnergyMap::const_iterator emap_citer = (const_cast< EnergyMap const & > (emap)).begin();
		TS_ASSERT_EQUALS( *emap_citer, sin( 1 ) );
	}

	/// @brief test that the end iterator returned is one element
	/// after the end of the list from 1 to n_score_types
	void test_EnergyMap_end() {
		EnergyMap emap;
		EnergyMap::const_iterator emap_iter = emap.begin();
		EnergyMap::const_iterator emap_end = emap.end();
		for ( Size ii = 1; ii < n_score_types; ++ii ) {
			++emap_iter;
		}
		TS_ASSERT_DIFFERS( emap_iter, emap_end );
		++emap_iter;
		TS_ASSERT_EQUALS( emap_iter, emap_end );
	}

	/// @brief test that the constant end iterator returned is one element
	/// after the end of the list from 1 to n_score_types
	void test_EnergyMap_constend() {
		EnergyMap emap;
		EnergyMap::const_iterator emap_citer = (const_cast< EnergyMap const & > (emap)).begin();
		EnergyMap::const_iterator emap_cend = (const_cast< EnergyMap const & > (emap)).end();
		for ( Size ii = 1; ii < n_score_types; ++ii ) {
			++emap_citer;
		}
		TS_ASSERT_DIFFERS( emap_citer, emap_cend );
		++emap_citer;
		TS_ASSERT_EQUALS( emap_citer, emap_cend );
	}

	/// @brief test that getting a value returns what's expected
	void test_EnergyMap_get() {
		EnergyMap emap;
		initialize_energymap_sin( emap );
		for ( Size ii = 1; ii <= n_score_types; ++ii ) {
			TS_ASSERT_EQUALS( sin( ii ), emap.get( (ScoreType) ii ) );
		}
	}

	/// @brief Test that setting a value happens correctly
	void test_EnergyMap_set() {
		EnergyMap emap;
		EnergyMap emap_compare;
		initialize_energymap_sin( emap );
		initialize_energymap_sin( emap_compare );
		for ( Size ii = 1; ii <= n_score_types; ++ii ) {
			emap_compare[ (ScoreType) ii ] += 2; //move emap_compare outside of range of emap
			TS_ASSERT_DIFFERS( emap_compare[ (ScoreType) ii ], emap[ (ScoreType) ii ] );
			emap.set( (ScoreType) ii, sin(ii) + 2 );
			TS_ASSERT_DELTA( emap_compare[ (ScoreType) ii ], emap[ (ScoreType) ii ], 1E10-6 );
			/// There is more precision in a double than 10^-6...
			/// however this unit test does not presume that
			/// Real is typedefed to double.
		}
	}

	/// @brief test that the [] operator works correctly.  If this test fails,
	/// most other tests will fail.  Relies on get() functionality.  get() unit test
	/// also relies on [] functionality -- circular dependency unavoidable.
	void test_EnergyMap_operator_square_brackets_nonconst() {
		EnergyMap emap;
		for ( Size ii = 1; ii <= n_score_types; ++ii ) {
			emap[ (ScoreType) ii ] = sin( ii );
		}

		for ( Size ii = 1; ii <= n_score_types; ++ii ) {
			TS_ASSERT_EQUALS( emap.get( (ScoreType) ii ), sin( ii ));
			TS_ASSERT_EQUALS( emap[ (ScoreType) ii ], sin( ii ));
		}

	}

	/// @brief test that the const [] operator works correctly.  If this test fails,
	/// most other tests will fail. Relies on non-const [] operator.
	void test_EnergyMap_operator_square_brackets_const() {
		EnergyMap emap;
		for ( Size ii = 1; ii <= n_score_types; ++ii ) {
			emap[ (ScoreType) ii ] = sin( ii );
		}

		for ( Size ii = 1; ii <= n_score_types; ++ii ) {
			TS_ASSERT_EQUALS( (const_cast< EnergyMap const & > (emap) )[ (ScoreType) ii ], sin( ii ));
		}

	}

	/// @brief test that the zero( ScoreTypes ) method zeros out exactly the right valutes
	/// relies on the correct behavior of the ScoreTypes class.
	void test_EnergyMap_zero_with_params() {
		EnergyMap emap;
		initialize_energymap_sin( emap );
		ScoreTypes scoretypes_to_zero;
		for ( Size ii = 1; ii <= n_score_types; ++ii ) {
			if ( ii % 2 == 1 ) {
				scoretypes_to_zero.push_back( (ScoreType) ii );
			}
		}

		emap.zero( scoretypes_to_zero );
		for ( Size ii = 1; ii <= n_score_types; ++ii ) {
			if ( ii % 2 == 1 ) {
				TS_ASSERT_EQUALS( emap[ (ScoreType) ii ], 0.0 );
			} else {
				TS_ASSERT_EQUALS( emap[ (ScoreType) ii ], sin( ii ) );
			}
		}

	}

	/// @brief test that the zero() method assigns all values in the EnergyMap to 0.0
	void test_EnergyMap_zero() {
		EnergyMap emap;
		initialize_energymap_sin( emap );
		emap.zero();
		for ( Size ii = 1; ii <= n_score_types; ++ii ) {
			TS_ASSERT_EQUALS( emap[ (ScoreType) ii ], 0.0 );
		}
	}

	/// @brief test that the clear() method assigns all values in the EnergyMap to 0.0
	void test_EnergyMap_clear() {
		EnergyMap emap;
		initialize_energymap_sin( emap );
		emap.clear();
		for ( Size ii = 1; ii <= n_score_types; ++ii ) {
			TS_ASSERT_EQUALS( emap[ (ScoreType) ii ], 0.0 );
		}
	}

	void test_EnergyMap_dot() {
		EnergyMap emap1;
		EnergyMap emap2;
		initialize_energymap_sin( emap1 );
		initialize_energymap_cos( emap2 );
		Real dot_product_1_dot_2 = emap1.dot( emap2 );
		Real dot_product_2_dot_1 = emap2.dot( emap1 );
		TS_ASSERT_DELTA( dot_product_1_dot_2, dot_product_2_dot_1, 1E-6 );

		Real expected_dot( 0.0 );
		for ( Size ii = 1; ii <= n_score_types; ++ii ) {
			expected_dot += cos( ii ) * sin( ii );
		}
		TS_ASSERT_DELTA( dot_product_1_dot_2, expected_dot, 1E-6 );
		TS_ASSERT_DELTA( dot_product_2_dot_1, expected_dot, 1E-6 );
	}

	/// @brief test that operator += accumulates the contents of one energy map into another,
	/// leaving the rhs untouched
	void test_EnergyMap_plus_equals_operator() {
		EnergyMap emap1;
		EnergyMap emap2;
		initialize_energymap_sin( emap1 );
		initialize_energymap_cos( emap2 );

		emap1 += emap2;
		for ( Size ii = 1; ii <= n_score_types; ++ii ) {
			TS_ASSERT_DELTA( emap1[ (ScoreType) ii ], cos(ii) + sin(ii), 1E-6 );
			TS_ASSERT_DELTA( emap2[ (ScoreType) ii ], cos(ii), 1E-6 );
		}
	}

	/// @brief test that operator -= decrements the contents of one energy map from another,
	/// leaving the rhs untouched
	void test_EnergyMap_minus_equals_operator() {
		EnergyMap emap1;
		EnergyMap emap2;
		initialize_energymap_sin( emap1 );
		initialize_energymap_cos( emap2 );

		emap1 -= emap2;
		for ( Size ii = 1; ii <= n_score_types; ++ii ) {
			TS_ASSERT_DELTA( emap1[ (ScoreType) ii ], sin(ii) - cos(ii), 1E-6 );
			TS_ASSERT_DELTA( emap2[ (ScoreType) ii ], cos(ii), 1E-6 );
		}
	}

	/// @brief test that operator *= multiplies the contents of one energy map into another,
	/// leaving the rhs untouched
	void test_EnergyMap_times_equals_operator() {
		EnergyMap emap1;
		EnergyMap emap2;
		initialize_energymap_sin( emap1 );
		initialize_energymap_cos( emap2 );

		emap1 *= emap2;
		for ( Size ii = 1; ii <= n_score_types; ++ii ) {
			TS_ASSERT_DELTA( emap1[ (ScoreType) ii ], sin(ii) * cos(ii), 1E-6 );
			TS_ASSERT_DELTA( emap2[ (ScoreType) ii ], cos(ii), 1E-6 );
		}
	}

	/// @brief test that accumulation from one energy map into another over a subset of the
	/// score types accumulates exactly those score types, and leaves the input emap untouched
	void test_EnergyMap_accumulate() {
		EnergyMap emap1, emap2;
		initialize_energymap_sin( emap1 );
		initialize_energymap_cos( emap2 );

		ScoreTypes scoretypes_to_accumulate;
		for ( Size ii = 1; ii <= n_score_types; ++ii ) {
			if ( ii % 2 == 1 ) {
				scoretypes_to_accumulate.push_back( (ScoreType) ii );
			}
		}

		emap1.accumulate( emap2, scoretypes_to_accumulate );
		for ( Size ii = 1; ii <= n_score_types; ++ii ) {
			if ( ii % 2 == 1 ) {
				TS_ASSERT_DELTA( emap1[ (ScoreType) ii ], sin( ii ) + cos( ii ), 1E-6 );
			} else {
				TS_ASSERT_EQUALS( emap1[ (ScoreType) ii ], sin( ii ) );
			}
			TS_ASSERT_EQUALS( emap2[ (ScoreType) ii ], cos(ii) );
		}
	}

	/// @brief test that accumulation with a weight factor from one energy map into
	/// another over a subset of the score types accumulates and scales exactly
	/// those score types, and leaves the input emap untouched
	void test_EnergyMap_accumulate_with_weight_factor() {
		EnergyMap emap1, emap2;
		initialize_energymap_sin( emap1 );
		initialize_energymap_cos( emap2 );

		ScoreTypes scoretypes_to_accumulate;
		for ( Size ii = 1; ii <= n_score_types; ++ii ) {
			if ( ii % 2 == 1 ) {
				scoretypes_to_accumulate.push_back( (ScoreType) ii );
			}
		}

		emap1.accumulate( emap2, scoretypes_to_accumulate, 0.5 );
		for ( Size ii = 1; ii <= n_score_types; ++ii ) {
			if ( ii % 2 == 1 ) {
				TS_ASSERT_DELTA( emap1[ (ScoreType) ii ], sin( ii ) + 0.5 * cos( ii ), 1E-6 );
			} else {
				TS_ASSERT_EQUALS( emap1[ (ScoreType) ii ], sin( ii ) );
			}
			TS_ASSERT_EQUALS( emap2[ (ScoreType) ii ], cos(ii) );
		}
	}

	/// @brief test that norm_squared computes the square magnitude of an emap vector.
	void test_EnergyMap_norm_squared() {
		EnergyMap emap;
		initialize_energymap_sin( emap );
		Real expected_normsquared( 0.0 );

		ScoreTypes scoretypes_to_square;
		for ( Size ii = 1; ii <= n_score_types; ++ii ) {
			if ( ii % 2 == 1 ) {
				scoretypes_to_square.push_back( (ScoreType) ii );
			}
		}


		for ( Size ii = 1; ii <= n_score_types; ++ii ) {
			if ( ii % 2 == 1 ) {
				expected_normsquared += sin(ii)*sin(ii);
			}
		}

		TS_ASSERT_DELTA( emap.norm_squared(scoretypes_to_square), expected_normsquared, 1E-6 );
	}
};
