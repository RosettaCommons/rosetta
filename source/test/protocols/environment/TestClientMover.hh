// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <protocols/environment/ClientMover.hh>
#include <protocols/environment/DofUnlock.hh>

#include <protocols/environment/claims/EnvClaim.hh>

#include <core/environment/DofPassport.hh>

#include <core/pose/Pose.hh>

#include <boost/bind/bind.hpp>

namespace protocols {
namespace environment {

class Tester : public protocols::environment::ClientMover {
public:

  Tester() : claim_( /* NULL */ ) {}

  void init( protocols::environment::claims::EnvClaimOP claim ) {
	claim_ = claim;
    claim_->set_owner( utility::pointer::dynamic_pointer_cast< protocols::environment::ClientMover >(get_self_ptr()) );
  }

  virtual void apply( core::pose::Pose& ){}

  // This use of bool apply allows us to apply an aribtrary functi
	using ClientMover::apply;
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

  void claim( protocols::environment::claims::EnvClaimOP claim ){
    claim_ = claim;
  }

	protocols::environment::claims::EnvClaimOP claim() {
		return claim_;
	}

  virtual std::string get_name() const { return "TESTER"; }

private:
  protocols::environment::claims::EnvClaimOP claim_;
};

typedef utility::pointer::shared_ptr< Tester > TesterOP;
typedef utility::pointer::shared_ptr< Tester const > TesterCOP;

}
}
