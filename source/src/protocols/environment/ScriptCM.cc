// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/ScriptCM.cc
/// @author Justin R. Porter
/// @author Brian D. Weitzner
/// @author Oliver F. Lange

// Unit Headers
#include <protocols/environment/ScriptCM.hh>
#include <protocols/environment/ScriptCMCreator.hh>

// Package headers

// Project headers
#include <core/environment/DofPassport.hh>
#include <protocols/environment/DofUnlock.hh>

#include <core/kinematics/MoveMap.hh>

#include <protocols/moves/MoverFactory.hh>


//Utility Headers
#include <utility/tag/Tag.hh>
#include <numeric/random/random.hh>

// tracer
#include <basic/Tracer.hh>

#ifdef WIN32
  #include <basic/datacache/WriteableCacheableMap.hh>
#endif

// C++ Headers

// ObjexxFCL Headers

namespace protocols {
namespace environment {

std::string const NOT_SET = "[NOT_SET]";

static thread_local basic::Tracer tr( "protocols.environment.ScriptCM", basic::t_info );

using namespace core::environment;
using namespace protocols::environment;

// creator
std::string
ScriptCMCreator::keyname() const {
  return ScriptCMCreator::mover_name();
}

protocols::moves::MoverOP
ScriptCMCreator::create_mover() const {
  return new ScriptCM;
}

std::string
ScriptCMCreator::mover_name() {
  return "ScriptCM";
}

ScriptCM::ScriptCM():
  ClaimingMover(),
  name_( NOT_SET ),
  client_( NULL )
{}

void ScriptCM::passport_updated(){
  core::kinematics::MoveMapOP mm = new core::kinematics::MoveMap();

  if( has_passport() ){
    passport()->render_movemap( mm );
  } else {
    mm->set_bb( true );
    mm->set_chi( true );
    mm->set_jump( true );
  }

  client_->set_movemap( mm );

}

void ScriptCM::initialize( core::pose::Pose& pose ){
  DofUnlock activation( pose.conformation(), passport() );
  client_->initialize( pose );
}

void ScriptCM::apply( core::pose::Pose& pose ){
  DofUnlock activation( pose.conformation(), passport() );
  client_->apply( pose );
}

void ScriptCM::parse_my_tag( utility::tag::TagCOP tag,
                             basic::datacache::DataMap& datamap,
                             protocols::filters::Filters_map const& filters,
                             protocols::moves::Movers_map const& mover_map,
                             core::pose::Pose const& pose ) {
  name_ = tag->getOption< std::string >( "name" );

  for( utility::vector0< utility::tag::TagCOP>::const_iterator tag_it = tag->getTags().begin();
      tag_it != tag->getTags().end(); ++tag_it ){
    TagCOP subtag = *tag_it;
    if( (*tag_it)->getName() == "Mover" ){
      std::string const& client_name = subtag->getOption< std::string >( "name" );
      tr.Debug << " Interpreting tag with name " << subtag->getName() << " as existing mover with name '"
               << client_name << "'" << std::endl;
      if( mover_map.find( client_name ) != mover_map.end() ){
        set_client( mover_map.find( client_name )->second );
      } else {
        throw utility::excn::EXCN_RosettaScriptsOption( "Undefined mover '"+client_name+"'." );
      }
    } else if( claims::EnvClaim::is_claim( subtag->getName() ) ) {
      tr.Debug << " Interpreting tag with name " << subtag->getName() << " as new claim." << std::endl;
      add_claim( claims::EnvClaim::make_claim( subtag->getName(), utility::pointer::static_pointer_cast< ClaimingMover >( get_self_ptr() ), subtag, datamap ) );
    } else {
      tr.Debug << " Interpreting tag with name " << subtag->getName() << " as a new mover." << std::endl;

      set_client( moves::MoverFactory::get_instance()->newMover( subtag, datamap, filters, mover_map, pose ) );

      if( subtag->hasOption( "name" ) ){
        tr.Warning << "[WARNING] Mover " << subtag->getOption< std::string >( "name" ) << " will not be availiable to"
                   << " reference by name. It will exist only within " << this->get_name() << std::endl;
      }
    }
  }
}

void ScriptCM::set_client( moves::MoverOP mover_in ) {
  if( client_ ){
    throw utility::excn::EXCN_RosettaScriptsOption( "The ScriptCM '" + this->get_name() + "' cannot contain >1 client mover." );
  }

  moves::MoveMapMoverOP mover_ptr = dynamic_cast< moves::MoveMapMover * >( mover_in.get() );

  if ( !mover_ptr ){
    throw utility::excn::EXCN_RosettaScriptsOption( "The "+mover_in->type()+" named '"+mover_in->get_name()+
                                                    "' doesn't implement MoveMapMover and can't be used by the ScriptCM." );
  }

  client_ = mover_ptr;
}

void ScriptCM::add_claim( claims::EnvClaimOP claim ) {
  tr.Debug << "  " << this->get_name() << " added EnvClaim " << *claim << "added." << std::endl;
  claim_list_.push_back( claim );
}


claims::EnvClaims ScriptCM::yield_claims( core::pose::Pose const&,
                                          basic::datacache::WriteableCacheableMapOP ){
  tr.Debug << this->get_name() << " yielding " << claim_list_.size() << " EnvClaims." << std::endl;
  return claim_list_;
}

std::string ScriptCM::get_name() const {
  return "ScriptCM("+name()+")";
}

moves::MoverOP ScriptCM::fresh_instance() const {
  return ClaimingMoverOP( new ScriptCM() );
}

moves::MoverOP ScriptCM::clone() const{
  return ClaimingMoverOP( new ScriptCM( *this ) );
}

} // environment
} // protocols
