// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file topo_driver.cc
/// @brief  Sandbox executable for the Environment and associated classes.
/// @author Justin Porter

#include <basic/Tracer.hh>


#include <basic/Tracer.hh>

#include <numeric/xyzVector.hh>

#include <core/id/types.hh>

#include <core/pack/task/residue_selector/ChainSelector.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/kinematics/Jump.hh>
#include <core/conformation/Residue.hh>

#include <protocols/environment/Environment.hh>
#include <protocols/environment/CoMTrackerCM.hh>
#include <protocols/environment/EnvExcn.hh>
#include <protocols/environment/DofUnlock.hh>

#include <protocols/environment/claims/EnvClaim.hh>
#include <protocols/environment/claims/XYZClaim.hh>

#include <protocols/abinitio/abscript/StructPerturberCM.hh>

#include <devel/init.hh>

#include <boost/bind/bind.hpp>

#define TS_ASSERT( x ) assert( x )
#define TS_ASSERT_THROWS_NOTHING( x ) x
#define TS_ASSERT_EQUALS( x , y ) TS_ASSERT( x == y )
#define TS_ASSERT_DIFFERS( x , y ) TS_ASSERT( x != y )
#define TS_ASSERT_DELTA( x , y, d ) TS_ASSERT( std::abs( x - y ) < d )
#define TS_ASSERT_THROWS( x , y ) try{ x; } catch( y ){}
#define TS_ASSERT_LESS_THAN( x, y ) TS_ASSERT( x < y );
#define TS_TRACE( x ) std::cout << x << std::endl;

static thread_local basic::Tracer tr( "main" );

class Tester : public protocols::environment::ClaimingMover {
public:

  Tester() : claim_( NULL ) {}
  Tester( protocols::environment::claims::EnvClaimOP claim ): claim_( claim ) {
    claim_->set_owner( this );
  }

  virtual void apply( core::pose::Pose& pose ){}

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

void test_replace_residue( core::pose::Pose const& pose ) {

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

  TS_ASSERT_THROWS( noclaim_test->apply( prot_pose, boost::bind( repres_apvect, &prot_pose, SEQPOS, new_diff_rsd, ap_vect ) ),
                    protocols::environment::EXCN_Env_Security_Exception );

  for( Size i = 1; i <= prot_pose.residue( SEQPOS ).natoms(); ++i ){
    TS_ASSERT_LESS_THAN( ( new_same_rsd.xyz( i ) - pose.residue( SEQPOS ).xyz( i ) ).length(), 1e-6 );
  }


  TS_ASSERT_THROWS( noclaim_test->apply( prot_pose, boost::bind( repres_bool, &prot_pose, SEQPOS, new_diff_rsd, true ) ),
                    protocols::environment::EXCN_Env_Security_Exception );
  for( Size i = 1; i <= prot_pose.residue( SEQPOS ).natoms(); ++i ){
    TS_ASSERT_LESS_THAN( ( new_same_rsd.xyz( i ) - pose.residue( SEQPOS ).xyz( i ) ).length(), 1e-6 );
  }

  //Verify: if we change something WITH claiming, we pass
  TS_ASSERT_THROWS_NOTHING( xyzclaim_test->apply( prot_pose, boost::bind( repres_bool, &prot_pose, SEQPOS, new_diff_rsd, true ) ) );
}

int main( int argc, char** argv ){
  devel::init( argc, argv );

  core::pose::Pose pose;
  core::pose::make_pose_from_sequence(pose, "FRIENDLYFRIENDS", "fa_standard");

//  for( Size i = 1; i <= pose.total_residue(); ++i ){
//    pose.set_phi( i, -65 );
//    pose.set_psi( i, -41 );
//  }
//
//  pose.append_pose_by_jump( *pose.clone(), 5 );
//  core::kinematics::Jump j = pose.jump( 1 );
//  j.gaussian_move(1, 50, 0);
//  pose.set_jump(1, j);

  try {
    test_replace_residue( pose );
  } catch ( utility::excn::EXCN_Msg_Exception excn ){
    std::cout << excn << std::endl;
    std::exit( 1 );
  }

  std::cout << "Done!" << std::endl;

}
