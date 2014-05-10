// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/abinitio/abscript/StructPerturberCMCM.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/abinitio/abscript/StructPerturberCM.hh>
#include <protocols/abinitio/abscript/StructPerturberCMCreator.hh>

// Package headers
#include <core/environment/DofPassport.hh>

#include <protocols/environment/DofUnlock.hh>
#include <protocols/environment/claims/EnvClaim.hh>
#include <protocols/environment/claims/TorsionClaim.hh>

#include <core/kinematics/MoveMap.hh>
// Project headers
#include <core/id/TorsionID.hh>
#include <core/id/DOF_ID.hh>


//Utility Headers
#include <utility/tag/Tag.hh>
#include <numeric/random/random.hh>

// tracer
#include <basic/Tracer.hh>

// C++ Headers

// ObjexxFCL Headers

static numeric::random::RandomGenerator RG(249845846);
static basic::Tracer tr("protocols.abinitio.abscript.StructPerturberCM", basic::t_info);

namespace protocols {
namespace abinitio {
namespace abscript {

using namespace core::environment;
using namespace protocols::environment;

// creator
std::string
StructPerturberCMCreator::keyname() const {
  return StructPerturberCMCreator::mover_name();
}

protocols::moves::MoverOP
StructPerturberCMCreator::create_mover() const {
  return new StructPerturberCM;
}

std::string
StructPerturberCMCreator::mover_name() {
  return "StructPerturberCM";
}

StructPerturberCM::StructPerturberCM():
  Parent()
{}

claims::EnvClaims StructPerturberCM::yield_claims( core::pose::Pose const& in_pose,
                                                   basic::datacache::WriteableCacheableMapOP ){
  claims::EnvClaims claims;

  claims::TorsionClaimOP claim = new claims::TorsionClaim( this, label(), std::make_pair( 1, in_pose.total_residue() ) );
  claim->ctrl_strength( claims::CAN_CONTROL );
  claim->init_strength( claims::DOES_NOT_INITIALIZE );

  claims.push_back( claim );

  return claims;
}

void StructPerturberCM::parse_my_tag( utility::tag::TagCOP tag,
                                      basic::datacache::DataMap&,
                                      protocols::filters::Filters_map const&,
                                      protocols::moves::Movers_map const&,
                                      core::pose::Pose const& ){
  magnitude( tag->getOption< core::Real >( "magnitude", 2.0 ) );
  label( tag->getOption< std::string >( "label", "BASE" ) );
}

void StructPerturberCM::apply( core::pose::Pose& pose ){
  if( passport() ){
    DofUnlock unlock( pose.conformation(), passport() );
    core::kinematics::MoveMapOP mm = passport()->render_movemap( pose.conformation() );

    for( Size i = 1; i <= pose.total_residue(); ++i ){
      core::id::DOF_ID phi = pose.conformation().dof_id_from_torsion_id( core::id::TorsionID( i, core::id::BB, core::id::phi_torsion ) );
      core::id::DOF_ID psi = pose.conformation().dof_id_from_torsion_id( core::id::TorsionID( i, core::id::BB, core::id::psi_torsion ) );

      tr.Trace << std::endl << "Purturbing resid " << i << ": ";
      if( passport()->dof_access( phi ) ){
        tr.Trace << "phi ( " << pose.phi( i ) << " -> ";
        pose.set_phi( i, pose.phi( i ) + RG.gaussian() * magnitude_ );
        tr.Trace << pose.phi( i ) << ", ";
      } if( passport()->dof_access( psi ) ) {
        tr.Trace << "psi ( " << pose.psi( i ) << " -> ";
        pose.set_psi( i, pose.psi( i ) + RG.gaussian() * magnitude_ );
        tr.Trace << pose.psi( i );
      }
      tr.Trace << std::endl;
    }
  } else {
    for( Size i = 1; i <= pose.total_residue(); ++i ){
        pose.set_phi( i, pose.phi( i ) + RG.gaussian() * magnitude_ );
        pose.set_psi( i, pose.psi( i ) + RG.gaussian() * magnitude_ );
    }
  }
}

std::string StructPerturberCM::get_name() const {
  return "StructPerturberCM";
}

moves::MoverOP StructPerturberCM::clone() const {
  return new StructPerturberCM( *this );
}

} // abscript
} // abinitio
} // protocols
