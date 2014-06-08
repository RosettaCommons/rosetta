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

//Utility Headers
#include <utility/tag/Tag.hh>
#include <numeric/random/random.hh>

// tracer
#include <basic/Tracer.hh>

// C++ Headers

// ObjexxFCL Headers

namespace protocols {
namespace environment {

static basic::Tracer tr("protocols.environment.ScriptCM", basic::t_info);
static numeric::random::RandomGenerator RG(9143235);

using namespace core::environment;
using namespace protocols::environment;

// creator
std::string
ScriptCMCreator::keyname() const {
  return ScriptCMCreator::mover_name();
}

protocols::moves::MoverOP
ScriptCMCreator::create_mover() const {
  return ClaimingMoverOP( new ScriptCM );
}

std::string
ScriptCMCreator::mover_name() {
  return "ScriptCM";
}

ScriptCM::ScriptCM():
  ClaimingMover(),
  name_( "[NOT_SET]" ),
  client_( NULL )
{}

void ScriptCM::passport_updated(){
  core::kinematics::MoveMapOP mm;
  if( passport() ){
    mm = passport()->render_movemap();
    // client_->set_movemap( mm );
  } else {
    mm = new core::kinematics::MoveMap();
    mm->set_bb( true );
    mm->set_chi( true );
    mm->set_jump( true );
  }

  // client_->set_movemap( mm );

}

void ScriptCM::initialize( core::pose::Pose& pose ){
  DofUnlock activeation( pose.conformation(), passport() );
  // client_->initialize( pose );
}

void ScriptCM::apply( core::pose::Pose& pose ){
  DofUnlock activation( pose.conformation(), passport() );
  client_->apply( pose );
}

void ScriptCM::parse_my_tag( utility::tag::TagCOP tag,
                             basic::datacache::DataMap&,
                             protocols::filters::Filters_map const&,
                             protocols::moves::Movers_map const& mover_map,
                             core::pose::Pose const& ) {
  name_ = tag->getOption< std::string >( "name" );

  for( utility::vector0< utility::tag::TagCOP>::const_iterator tag_it = tag->getTags().begin();
      tag_it != tag->getTags().end(); ++tag_it ){
    TagCOP subtag = *tag_it;
    if( (*tag_it)->getName() == "Mover" ){
      std::string const& client_name = subtag->getOption< std::string >( "name" );
      if( mover_map.find( client_name ) != mover_map.end() ){
        client_ = mover_map.find( client_name )->second;
      } else {
        throw utility::excn::EXCN_RosettaScriptsOption( "The "+type()+" named "+name()+
                                                        " couldn't find the mover named "+client_name+"." );
      }
    } else {
      //TODO: parse subtags
    }
  }
}

claims::EnvClaims ScriptCM::yield_claims( core::pose::Pose const&,
                                          basic::datacache::WriteableCacheableMapOP ){
  claims::EnvClaims claim_list;

  return claim_list;
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
