// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <core/environment/DofPassport.hh>

#include <protocols/environment/claims/EnvClaim.fwd.hh>
#include <protocols/environment/claims/TorsionClaim.hh>
#include <protocols/environment/claims/JumpClaim.hh>

#include <protocols/environment/claims/BrokerElements.hh>

#include <protocols/environment/ProtectedConformation.hh>
#include <protocols/environment/EnvExcn.hh>
#include <protocols/environment/DofUnlock.hh>
#include <protocols/environment/ClientMover.hh>
#include <protocols/environment/Environment.hh>
#include <protocols/environment/EnvClaimBroker.hh>

//Other headers
#include <core/conformation/Conformation.hh>

#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/types.hh>

#include <test/core/init_util.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/WriteableCacheableMap.hh>
#include <core/pose/datacache/CacheableDataType.hh>


using namespace protocols::environment;

// ---------------- Toy Movers --------------- //

namespace protocols {
namespace environment {

const core::Size CLAIMED_RESID = 10;
const core::Size UNCLAIMED_RESID = 9;
const core::Real NEW_PHI = 23.0;

const core::Size JUMP_START = 3;
const core::Size JUMP_END = 7;
class ToyMover : public protocols::environment::ClientMover {
protected:

  ToyMover( bool claim, bool move ):
    claim_( claim ),
    move_( move )
  {}

  //Yes I know this is against the law
  bool claim_;
  bool move_;
};

class TorsionMover : public ToyMover {
public:
  TorsionMover( bool claim, bool move,
    claims::ControlStrength c_str = claims::MUST_CONTROL,
    claims::ControlStrength i_str = claims::DOES_NOT_CONTROL ):
    ToyMover( claim, move ),
    resnum_( CLAIMED_RESID ),
    control_str_( c_str ),
    init_str_( i_str )
  {}

	using ToyMover::apply;
  void apply( Pose& pose, Size resid ){
    DofUnlock activation( pose.conformation(), passport() );
    if( move_ ) {
      pose.set_phi( resid, NEW_PHI );
    }
  }

  virtual void apply( Pose& pose ){
    apply( pose, resnum_ );
  }

  virtual void missing_unlock_apply( Pose& pose ){
    if( move_ ){
      pose.set_phi( resnum_, NEW_PHI );
    }
  }

  virtual std::string get_name() const{
    return "TorsionMover";
  }

  virtual claims::EnvClaims yield_claims( core::pose::Pose const&,
            basic::datacache::WriteableCacheableMapOP ){
    using core::environment::LocalPosition;
    claims::EnvClaims claims;

    if( claim_ ){
      claims::TorsionClaimOP new_claim( new claims::TorsionClaim(
		utility::pointer::dynamic_pointer_cast< protocols::environment::ClientMover >(get_self_ptr()),
		LocalPosition( "BASE", resnum_ ) ) );
      new_claim->strength( control_str_, init_str_ );
      claims.push_back( new_claim );
    }

    return claims;
  }

  virtual void initialize( core::pose::Pose& pose ){
    DofUnlock activation( pose.conformation(), passport() );
    if( move_ ){
      pose.set_phi( resnum_, NEW_PHI );
    }
  }

  Size resnum() {
    return resnum_;
  }

private:
  Size resnum_;
  claims::ControlStrength control_str_;
  claims::ControlStrength init_str_;
};

typedef utility::pointer::shared_ptr< TorsionMover > TorsionMoverOP;

} //environment
} //protocols
