// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/loop_mover/LoopCM.cc
/// @brief  loop mover base class
/// @author Justin R. Porter

// Unit headers
#include <protocols/loops/loop_mover/LoopCM.hh>
#include <protocols/loops/loop_mover/LoopCMCreator.hh>

// Package Headers
#include <protocols/loops/loop_mover/perturb/LoopMover_CCD.hh>
#include <protocols/loops/loop_mover/perturb/LoopMover_KIC.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_CCD.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_KIC.hh>

#include <protocols/loops/Loops.hh>
#include <protocols/loops/Loop.hh>

#include <protocols/environment/claims/TorsionClaim.hh>
#include <protocols/environment/claims/JumpClaim.hh>
#include <protocols/environment/DofUnlock.hh>
#include <protocols/environment/EnvExcn.hh>

#include <core/environment/LocalPosition.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/pack/task/residue_selector/ResidueSelector.hh>

//Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

#include <basic/datacache/DataMap.hh>

/// ObjexxFCL headers

// C++ Headers

namespace protocols {
namespace loops {
namespace loop_mover {

static thread_local basic::Tracer tr( "protocols.loops.loop_mover.LoopCM", basic::t_info );

// creator
std::string
LoopCMCreator::keyname() const {
  return LoopCMCreator::mover_name();
}

protocols::moves::MoverOP
LoopCMCreator::create_mover() const {
  return environment::ClaimingMoverOP( new LoopCM );
}

std::string
LoopCMCreator::mover_name() {
  return "LoopCM";
}

std::string const CCD( "CCD" );
std::string const KIC( "KIC" );
std::string const PERTURB( "perturb" );
std::string const REFINE( "refine" );
std::string const UNSET( "[NOT SET]" );

LoopCM::LoopCM() :
  Parent(),
  algorithm_( UNSET ),
  style_( UNSET )
{}

void LoopCM::parse_my_tag( TagCOP const tag,
                           basic::datacache::DataMap& datamap,
                           Filters_map const&,
                           moves::Movers_map const&,
                           Pose const& ){
  using core::pack::task::residue_selector::ResidueSelector;

  std::string const algorithm = tag->getOption< std::string >( "algorithm" );
  if( algorithm == KIC ){
    algorithm_ = KIC;
  } else if ( algorithm == CCD ){
    algorithm_ = CCD;
  } else {
    std::ostringstream ss;
    ss << "In '" << this->get_name() << "', the value '" << algorithm << "' of the option 'algorithm' is not valid. "
       << "Valid options are '" << KIC << "' and '" << CCD << "'." << std::endl;
    throw utility::excn::EXCN_RosettaScriptsOption( ss.str() );
  }

  std::string const style = tag->getOption< std::string >( "style" );
  if( style == PERTURB ){
    style_ = PERTURB;
  } else if( style == REFINE ){
    style_ = REFINE;
  } else {
    std::ostringstream ss;
    ss << "In '" << this->get_name() << "', the value '" << style << "' of the option 'style' is not valid. "
       << "Valid options are '" << PERTURB << "' and '" << REFINE << "'." << std::endl;
    throw utility::excn::EXCN_RosettaScriptsOption( ss.str() );
  }

  std::string const selector = tag->getOption< std::string >( "selector" );
  selector_ = datamap.get< ResidueSelector const* >( "ResidueSelector", selector );
}

void LoopCM::build_mover( LoopsOP loops ) {

  if( algorithm_ == UNSET ) {
    throw utility::excn::EXCN_BadInput( "Value of algorithm was unset in '" + get_name() + "'." );
  } else if( style_ == UNSET ) {
    throw utility::excn::EXCN_BadInput( "Value of style was unset in '" + get_name() + "'." );
  }

  if( algorithm_ == KIC ){
    if( style_ == PERTURB ){
      mover_ = new perturb::LoopMover_Perturb_KIC( loops );
    } else if( style_ == REFINE ){
      mover_ = new refine::LoopMover_Refine_KIC( loops );
    }
  } else if ( algorithm_ == CCD ){
    if( style_ == PERTURB ){
      mover_ = new perturb::LoopMover_Perturb_CCD( loops );
    } else if( style_ == REFINE ){
      mover_ = new refine::LoopMover_Refine_CCD( loops );
    }
  }

  assert( mover_ );
}

bool ang_delta( core::Real const& a, core::Real const& b ){
  core::Real const TOLERANCE = 1e-6;
  return std::abs( std::cos( a ) - std::cos( b ) ) > TOLERANCE;
}

void LoopCM::apply( core::pose::Pose& in_pose ){

  core::pose::Pose tmp_pose( in_pose );
  tmp_pose.set_new_conformation( core::conformation::ConformationOP( new core::conformation::Conformation( tmp_pose.conformation() ) ) );
  mover_->apply( tmp_pose );

  // Copy result into protected conformation in in_pose
  environment::DofUnlock unlock( in_pose.conformation(), passport() );
  core::Real const TOLERANCE = 1e-6;

  for ( Size i = 1; i <= in_pose.total_residue(); ++i ) {
    try {
      if( in_pose.residue( i ).is_protein() ){
        if( ang_delta( tmp_pose.omega( i ), in_pose.omega( i ) ) ){
          in_pose.set_omega( i, tmp_pose.omega( i ) );
        }
        if( ang_delta( tmp_pose.phi( i ), in_pose.phi( i ) ) ){
          in_pose.set_phi( i, tmp_pose.phi( i ) );
        }
        if( ang_delta( tmp_pose.psi( i ), in_pose.psi( i ) ) ) {
          in_pose.set_psi( i, tmp_pose.psi( i ) );
        }
        for( Size j = 1; j <= in_pose.conformation().residue( i ).nchi(); ++j ){
          if( ang_delta( in_pose.chi( (int) j, i ), tmp_pose.chi( (int) j, i ) ) > TOLERANCE){
            //TODO: figure out why this fails, and adjust mover parameters accordingly
//            in_pose.set_chi( (int) j, i, tmp_pose.chi( (int) j, i ) );
          }
        }
      }
    } catch ( environment::EXCN_Env_Security_Exception const& e ){
      tr.Error << "[ERROR] Unauthorized changes occurred during loop closure by mover'" << this->get_name()
               << "': (attempt to write to resid " << i << ")." << std::endl;
      throw e;
    }
  }

}

environment::claims::EnvClaims LoopCM::yield_claims( core::pose::Pose const& pose,
                                                     basic::datacache::WriteableCacheableMapOP ) {
  using namespace environment::claims;
  using namespace core::environment;
  EnvClaims claims;

  utility::vector1< bool > selection( pose.total_residue(), false );
  selector_->apply( pose, selection );

  LoopsOP loops = new Loops( selection );
  if( loops->empty() ){
    std::ostringstream ss;
    ss << "The mover " << get_name() << " couldn't build a loops object from the selection it was given." << std::endl;
    throw utility::excn::EXCN_BadInput( ss.str() );
  }
  build_mover( loops );

  core::environment::LocalPositions pos_list;

  Size const LOOP_EXTEND = 1;
  for( Loops::const_iterator loop = loops->loops().begin(); loop != loops->loops().end(); ++loop ){
    tr.Debug << " Claiming torsions in residues " << std::max( (Size) 1, loop->start() - LOOP_EXTEND )
             << "-" << std::min( pose.total_residue(), loop->stop() + LOOP_EXTEND ) << "." << std::endl;
    for( Size i = std::max( (Size) 1, loop->start()-LOOP_EXTEND );
         i <= std::min( pose.total_residue(), loop->stop()+LOOP_EXTEND ); ++i ){
      pos_list.push_back( new core::environment::LocalPosition( "BASE", i ) );
    }
  }

  environment::ClaimingMoverOP this_ptr = utility::pointer::static_pointer_cast< ClaimingMover >( get_self_ptr() );
  TorsionClaimOP torsion = new TorsionClaim( this_ptr, pos_list );
  torsion->claim_backbone( true );
  torsion->claim_sidechain( true );
  torsion->strength( MUST_CONTROL, DOES_NOT_CONTROL );

  claims.push_back( torsion );

  for( Size i = 1; i <= loops->num_loop(); ++i ){
    Loop const& loop = loops->loops()[i];
    JumpClaimOP jump = new JumpClaim( this_ptr,
                                      get_name()+"_"+utility::to_string( i ),
                                      LocalPosition( "BASE", loop.start() ),
                                      LocalPosition( "BASE", loop.stop() ),
                                      LocalPosition( "BASE", loop.cut() ) );
    jump->strength( EXCLUSIVE, EXCLUSIVE );
    jump->physical( false );
    claims.push_back( jump );
    tr.Debug << "    " << get_name() << " claiming Loop: " << loop.start() << "-" << loop.stop()
             << " w/ cut @ " << loop.cut() << std::endl;
  }

  return claims;
}

void LoopCM::passport_updated() {
  if( has_passport() ){
    core::kinematics::MoveMapOP mm = passport()->render_movemap();
    mover_->false_movemap( mm );
  } else if( mover_ ) {
    core::kinematics::MoveMapOP mm = new core::kinematics::MoveMap();
    mm->set_bb( true );
    mm->set_chi( true );
    mover_->false_movemap( mm );
  } else {
    tr.Trace << get_name() << " ignoring passport_updated call as it has no mover (yet, hopefully)." << std::endl;
  }
}

std::string LoopCM::get_name() const {
  return "LoopCM(" + algorithm_ + +"_" + style_ + ")";
}


} // namespace loop_mover
} // namespace loops
} // namespace protocols
