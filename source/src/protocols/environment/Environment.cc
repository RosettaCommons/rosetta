// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/Environment.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/environment/Environment.hh>

// Package headers
#include <protocols/environment/ClaimingMover.hh>
#include <protocols/environment/ProtectedConformation.hh>

#include <protocols/environment/EnvExcn.hh>

// Project headers
#include <core/environment/DofPassport.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/WriteableCacheableMap.hh>

#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>

#include <core/chemical/VariantType.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>

// tracer
#include <basic/Tracer.hh>

// C++ Headers
#include <set>
#include <boost/foreach.hpp>

// ObjexxFCL Headers

static basic::Tracer tr("protocols.environment.Environment", basic::t_info);

namespace protocols {
namespace environment {

Environment::Environment( std::string name ):
  Parent( name ),
  broker_( NULL )
{}

Environment::~Environment() {
  BOOST_FOREACH( ProtectedConformationAP pconf, pconfs_ ){
    pconf->env_destruction();
  }
}

void Environment::register_mover( ClaimingMoverOP mover ){
  if( !is_registered( mover ) ){
    registered_movers_.insert( mover );
    if( mover->is_loop_closer() ){
      register_loop_closer( mover );
    }
    std::set< ClaimingMoverOP > submovers;
    mover->yield_submovers( submovers );
    register_movers( submovers.begin(), submovers.end() );
  }
}

core::pose::Pose Environment::start( core::pose::Pose const& in_pose ){

  // This cast of reference to OP is ok because we know that Conformations are stored as pointers anyway.
  ConformationCOP in_conf = &in_pose.conformation();

  //figure out if we've got a superenvironment
  ProtectedConformation const* conf_ptr = dynamic_cast<ProtectedConformation const*>( in_conf.get() );
  if( conf_ptr ) {
    set_superenv( conf_ptr->environment().get() );
  } else {
    set_superenv( NULL );
  }

  tr.Debug << "Start environment: '" << name() << "'" << std::endl;

  return this->broker( in_pose );
}

core::pose::Pose Environment::end( core::pose::Pose const& pose ){
  ProtectedConformationCOP conf = dynamic_cast< ProtectedConformation const* >( &pose.conformation() );
  if( !conf ){
    tr.Error << "[ERROR] Environment::end recieved a pose that contains an unprotcted Conformation."
             << std::endl;
    throw utility::excn::EXCN_BadInput( "Nonprotected pose came in to Environment::end" );
  }

  core::pose::Pose new_pose = pose;

  // Setting a new conformation clears Pose::datacache, and we want to protect the WRITABLE_DATA.
  if( new_pose.data().has( core::pose::datacache::CacheableDataType::WRITEABLE_DATA ) ){

    if( ! broker_->result().pose.data().has( core::pose::datacache::CacheableDataType::WRITEABLE_DATA ) ){
      tr.Warning << "[WARNING] The brokered pose contained WriteableCacheable data in the Pose DataCache, but this "
                 << "data is no longer present in the DataCache at environment closing. This is probably "
                 << "because a call was made to pose::set_new_conformation somewhere (anywhere), which clears the "
                 << "cache. Certain features of one or more of the ClaimingMovers registered to the Environment '"
                 << name() << "' are nonfunctional." << std::endl;
    }

    basic::datacache::WriteableCacheableMapOP data_map = dynamic_cast< basic::datacache::WriteableCacheableMap* >( new_pose.data().get_raw_ptr( core::pose::datacache::CacheableDataType::WRITEABLE_DATA ) );
    new_pose.set_new_conformation( end( conf ) );
    new_pose.data().set( core::pose::datacache::CacheableDataType::WRITEABLE_DATA, data_map );
  } else {
    new_pose.set_new_conformation( end(conf ) );
  }

  tr.Debug << "Finish environment " << name() << " with sequence: " << new_pose.annotated_sequence()
           << " and fold tree: " << new_pose.fold_tree() << std::endl;

  return new_pose;
}

core::conformation::ConformationOP Environment::end( ProtectedConformationCOP conf ){
  tr.Debug << "End environment: '" << name() << "'" << std::endl;

  core::pose::Pose pose;

  ProtectedConformationOP outer_conf = new ProtectedConformation( *conf );
  pose.set_new_conformation( outer_conf );

  remove_nonpermenant_features( pose );

  core::conformation::ConformationOP ret_conf = new Conformation( pose.conformation() );
  ret_conf->fold_tree( *input_ft_ );

  // Reprotect Conformation if there's a superenvironment.
  if( superenv() ){
    ret_conf = new ProtectedConformation( superenv(), *ret_conf );
  }

  cancel_passports();

  return ret_conf;
}

void Environment::assign_passport( ClaimingMoverOP mover, core::environment::DofPassportCOP passport ){
  mover->push_passport( this , passport );
}

void Environment::cancel_passports(){
  BOOST_FOREACH( ClaimingMoverOP mover, registered_movers_ ){
    tr.Error << "stripping passport from " << mover->get_name() << std::endl;
    mover->pop_passport( this );
  }
}

  /// ENV OPEN/CLOSING
void Environment::remove_nonpermenant_features( core::pose::Pose& pose ){
  if( loop_closer_ ){
    loop_closer_->apply( pose );
  }

  //reset chainbreak varaints
  tr.Trace << "Removing chainbreak variants from sequence "
           << pose.annotated_sequence() << std::endl;
  BOOST_FOREACH( std::set< core::Size >::key_type cut_res,
                 broker()->result().inserted_cut_variants ) {
    remove_chainbreak_variants( pose, cut_res, cut_res+1 );
    tr.Trace << "Chainbreak variants @ " << cut_res << "/" << cut_res+1 << " removed:"
             << pose.annotated_sequence() << std::endl;
  }

  //remove virtual residues

  //Reset fold tree
}

core::pose::Pose Environment::broker( core::pose::Pose const& in_pose ){
  tr.Debug << "Beginning Broking for environment " << this->name() << " with "
           << registered_movers_.size() << " registered movers." << std::endl;

  std::map< ClaimingMoverOP, core::environment::DofPassportOP > mover_passports;
  input_ft_ = new core::kinematics::FoldTree( in_pose.fold_tree() );

  BOOST_FOREACH( ClaimingMoverOP mover, registered_movers_ ){
    assert( mover_passports.find( mover ) == mover_passports.end() );
    mover_passports[ mover ] = Parent::issue_passport( mover->get_name() );
    assign_passport( mover, mover_passports[ mover ] );
  }

  core::pose::Pose out_pose;
  try{
    broker_ = new EnvClaimBroker( *this, mover_passports, in_pose );
    out_pose = broker_->result().pose;
  } catch ( ... ){
    cancel_passports();
    throw;
  }

  return out_pose;
}

void Environment::remove_chainbreak_variants( core::pose::Pose& pose, core::Size down_res, core::Size up_res ) const {
  using namespace core::chemical;
  using namespace core::conformation;

  Conformation& conf = pose.conformation();

  Residue const& rsd_lower( conf.residue( down_res ) );
  Residue const& rsd_upper( conf.residue( up_res ) );
  ResidueTypeSet const& rsd_set( rsd_lower.residue_type_set() );

  ResidueType const& new_type_lower( rsd_set.get_residue_type_with_variant_removed( rsd_lower.type(),
                                                                                    CUTPOINT_LOWER ) );
  ResidueType const& new_type_upper( rsd_set.get_residue_type_with_variant_removed( rsd_upper.type(),
                                                                                    CUTPOINT_UPPER ) );

  ResidueOP new_lower( ResidueFactory::create_residue( new_type_lower, rsd_lower, conf ) );
  ResidueOP new_upper( ResidueFactory::create_residue( new_type_upper, rsd_upper, conf ) );

  copy_residue_coordinates_and_rebuild_missing_atoms( rsd_lower, *new_lower, conf );
  copy_residue_coordinates_and_rebuild_missing_atoms( rsd_upper, *new_upper, conf );

  conf.replace_residue( rsd_lower.seqpos(), *new_lower, true );
  conf.replace_residue( rsd_upper.seqpos(), *new_upper, true );
}

bool Environment::is_registered( ClaimingMoverOP mover ) const {
  return ( std::find( registered_movers_.begin(), registered_movers_.end(), mover ) != registered_movers_.end() );
}

void Environment::register_loop_closer( ClaimingMoverOP closer ) {
  if( !loop_closer_ ){
    this->loop_closer_ = closer;
    register_mover( closer );
  } else {
    tr.Error << "[ERROR] Multiple ClaimingMovers designated as LoopCloser ("
    << closer->get_name() << " and " << loop_closer_->get_name()
    << ") in Environment " << name() << std::endl;
    throw utility::excn::EXCN_BadInput( "LoopCloser " + closer->get_name() +
                                       " cannot be added because Environment " +
                                       this->name() + " already has a LoopCloser " +
                                       loop_closer_->get_name() );
  }
}

EnvironmentCAP Environment::superenv() const{
  if( Parent::superenv() ){
    EnvironmentCAP env = dynamic_cast< Environment const* >( Parent::superenv().get() );
    assert( env != 0 );
    return env;
  } else {
    return 0;
  }
}

} // environment
} // protocols
