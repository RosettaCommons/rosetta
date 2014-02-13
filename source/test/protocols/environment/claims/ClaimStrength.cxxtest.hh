// Test headers
#include <cxxtest/TestSuite.h>

// Package headers
#include <protocols/environment/claims/ClaimStrength.hh>

//C++ headers
#include <algorithm>
#include <vector>

#include <iostream>

// --------------- Test Class --------------- //

class ClaimStrengthTest : public CxxTest::TestSuite {

public:

  // Shared data elements go here.


  // --------------- Fixtures --------------- //

  // Shared initialization goes here.
  void setUp() {
  }

  // Shared finalization goes here.
  void tearDown() {
  }

  // --------------- Test Cases --------------- //
  /// @brief angle tests
  void test_ClaimStrength_sort (){

    using namespace protocols::environment::claims;

    std::vector< ClaimStrength > cs;

    ClaimStrength c1( ClaimStrength::CAN_CONTROL, -1 );
    ClaimStrength c2( ClaimStrength::CAN_CONTROL, 0 );
    ClaimStrength c3( ClaimStrength::CAN_CONTROL, 1 );
    ClaimStrength c4( ClaimStrength::MUST_CONTROL, 0 );
    ClaimStrength c5( ClaimStrength::EXCLUSIVE, 0 );

    std::cout << "WTF" << std::endl;

    cs[3] = c1;
    cs[4] = c2;
    cs[1] = c3;
    cs[2] = c4;
    cs[0] = c5;

    std::sort( cs.begin(), cs.end() );

    TS_ASSERT( cs[0] == c1 );
    TS_ASSERT( cs[1] == c2 );
    TS_ASSERT( cs[2] == c3 );
    TS_ASSERT( cs[3] == c4 );
    TS_ASSERT( cs[4] == c5 );
    TS_ASSERT( false );
    TS_FAIL( "Oh no!" );
  }

};
