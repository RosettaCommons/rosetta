// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:

// Test headers
#include <cxxtest/TestSuite.h>

// Package headers
#include <core/environment/DofPassport.hh>

#include <protocols/environment/claims/EnvClaim.hh>
#include <protocols/environment/claims/XYZClaim.hh>

#include <protocols/environment/ProtectedConformation.hh>
#include <protocols/environment/EnvExcn.hh>
#include <protocols/environment/DofUnlock.hh>
#include <protocols/environment/ClaimingMover.hh>
#include <protocols/environment/Environment.hh>

//Other headers
#include <core/conformation/Conformation.hh>

#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/types.hh>

#include <test/core/init_util.hh>

//C++ headers
#include <iostream>
#include <boost/bind/bind.hpp>

using namespace core; 
using namespace core::conformation; 

class Tester : public protocols::environment::ClaimingMover {
public:

  Tester() : claim_( NULL ) {}
  Tester( protocols::environment::claims::EnvClaimOP claim ): claim_( claim ) {
    claim_->set_owner( this );
  }

  virtual void apply( core::pose::Pose& ){}

  // This use of bool apply allows us to apply an aribtrary functi
  virtual void apply( core::pose::Pose& pose, boost::function< void() > f ) {
    protocols::environment::DofUnlock activation( pose.conformation(), passport() );
    f();
  }

  virtual protocols::environment::claims::EnvClaims yield_claims( core::pose::Pose const&,
                                                                  basic::datacache::WriteableCacheableMapOP ) {
    protocols::environment::claims::EnvClaims claims;
    if( claim_ ) claims.push_back( claim_ );
    return claims;
  }

  virtual std::string get_name() const { return "TESTER"; }

private:
  protocols::environment::claims::EnvClaimOP claim_;
};

typedef utility::pointer::owning_ptr< Tester > TesterOP;
typedef utility::pointer::owning_ptr< Tester const > TesterCOP;

// --------------- Test Class --------------- //

class ProtectedConf : public CxxTest::TestSuite {
public:

  // Shared data elements go here.
  core::pose::Pose pose;

  // --------------- Fixtures --------------- //

  // Shared initialization goes here.
  void setUp() {
    core_init();

    using namespace protocols::environment;
    using namespace core::environment;

    core::pose::make_pose_from_sequence(pose, "FRIENDLYFRIENDS", "fa_standard");
  }

  // Shared finalization goes here.
  void tearDown() {
  }

  void test_replace_residue() {

    typedef utility::vector1<std::pair<std::string, std::string> > AtomPairVector;
    using namespace protocols::environment::claims;

    core::Size const SEQPOS = 5;

    protocols::environment::Environment env( "env" );
    TesterOP noclaim_test = new Tester();
    XYZClaimOP claim = new XYZClaim( NULL, core::environment::LocalPosition( "BASE", SEQPOS ) );
    claim->strength( CAN_CONTROL, DOES_NOT_CONTROL );
    TesterOP xyzclaim_test = new Tester( claim );


    env.register_mover( noclaim_test );
    env.register_mover( xyzclaim_test );

    core::pose::Pose prot_pose;
    TS_ASSERT_THROWS_NOTHING( prot_pose = env.start( pose ) );

    core::conformation::Residue const& new_same_rsd( pose.residue( SEQPOS ) );
    AtomPairVector ap_vect;
    ap_vect.push_back( std::make_pair( "CA", "CA" ) );
    ap_vect.push_back( std::make_pair( "N", "N" ) );
    ap_vect.push_back( std::make_pair( "C", "C" ) );

    // Make aliases for the two versions of replace_residue (make a complicated test a *bit* more readable)
    void (core::pose::Pose::*repres_bool) ( Size const, Residue const&, bool const) = &core::pose::Pose::replace_residue;
    void (core::pose::Pose::*repres_apvect) ( int const, Residue const&, AtomPairVector const& ) = &core::pose::Pose::replace_residue;

    //Verify: if we don't chage anything with this call, we pass
    TS_ASSERT_THROWS_NOTHING( noclaim_test->apply( prot_pose, boost::bind( repres_apvect, &prot_pose, SEQPOS, new_same_rsd, ap_vect ) ) );
    TS_ASSERT_THROWS_NOTHING( noclaim_test->apply( prot_pose, boost::bind( repres_bool, &prot_pose, SEQPOS, new_same_rsd, true ) ) );

    //Verify: if we change something without claiming, we get an exception
    core::conformation::Residue new_diff_rsd( pose.residue( SEQPOS ) );
    new_diff_rsd.set_xyz( 1, numeric::xyzVector< Real >( 1.0, 1.0, 1.0 ) );
    new_diff_rsd.set_xyz( 2, numeric::xyzVector< Real >( 2.0, 2.0, 2.0 ) );

    TS_ASSERT_THROWS( noclaim_test->apply( prot_pose, boost::bind( repres_apvect, &prot_pose, SEQPOS, new_diff_rsd, ap_vect ) ),
                     protocols::environment::EXCN_Env_Security_Exception );

    for( Size i = 1; i <= prot_pose.residue( SEQPOS ).natoms(); ++i ){
      TS_ASSERT_LESS_THAN( ( prot_pose.residue( SEQPOS ).xyz( i ) - pose.residue( SEQPOS ).xyz( i ) ).length(), 1e-6 );
      TS_ASSERT_EQUALS( prot_pose.fold_tree(), pose.fold_tree() );
    }

    TS_ASSERT_THROWS( noclaim_test->apply( prot_pose, boost::bind( repres_bool, &prot_pose, SEQPOS, new_diff_rsd, true ) ),
                     protocols::environment::EXCN_Env_Security_Exception );
    for( Size i = 1; i <= prot_pose.residue( SEQPOS ).natoms(); ++i ){
      TS_ASSERT_LESS_THAN( ( new_same_rsd.xyz( i ) - pose.residue( SEQPOS ).xyz( i ) ).length(), 1e-6 );
      TS_ASSERT_EQUALS( prot_pose.fold_tree(), pose.fold_tree() );
    }

    //Verify: if we change something WITH claiming, we pass
    TS_ASSERT_THROWS_NOTHING( xyzclaim_test->apply( prot_pose, boost::bind( repres_bool, &prot_pose, SEQPOS, new_diff_rsd, true ) ) );
  }
};
