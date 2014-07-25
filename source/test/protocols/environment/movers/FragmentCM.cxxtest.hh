// Test headers
#include <cxxtest/TestSuite.h>

// Package headers


//Other headers
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/types.hh>

#include <test/core/init_util.hh>

//C++ headers
#include <iostream>

using namespace core; 

// ---------------- Toy Movers --------------- //

namespace protocols {
namespace environment {

// --------------- Test Class --------------- //

class FragmentCM : public CxxTest::TestSuite {
public:

  // Shared data elements go here.
  core::pose::Pose pose;
  utility::vector1< core::Real > init_phis;

  core::fragment::FragSetOP frags_small;
  core::fragment::FragSetOP frags_large;

  // --------------- Fixtures --------------- //

  // Shared initialization goes here.
  void setUp() {
    core_init();

    using namespace protocols::environment;
    using namespace core::environment;

    core::pose::make_pose_from_sequence(pose, "FRIENDLYFRIENDS", "fa_standard");

    //store phi values for comparison in the test
    for( core::Size i = 1; i <= pose.total_residue(); ++i ){
      init_phis.push_back( pose.phi( i ) );
    }

    //frags_small = frag_io.read_data( "" );
    //frags_large = frag_io.read_data( "" );
  }

  // Shared finalization goes here.
  void tearDown() {
  }

  // --------------- Test Cases --------------- //
  void test_large() {

  }


};
