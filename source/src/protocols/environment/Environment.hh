// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/Enironment.hh
/// @brief An environment that automatically distributes rights to shared degrees of freedom (e.g. fold tree)
///
/// @author Justin Porter

#ifndef INCLUDED_protocols_environment_Environment_hh
#define INCLUDED_protocols_environment_Environment_hh

// Unit Headers
#include <protocols/environment/Environment.fwd.hh>
#include <core/environment/EnvCore.hh>

// Package headers
#include <core/environment/DofPassport.fwd.hh>

#include <protocols/environment/ClaimingMover.fwd.hh>
#include <protocols/environment/ProtectedConformation.fwd.hh>
#include <protocols/environment/EnvClaimBroker.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Conformation.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <string>
#include <list>
#include <set>

// ObjexxFCL Headers

namespace protocols {
namespace environment {

class Environment : public core::environment::EnvCore {
  typedef core::environment::EnvCore Parent;

  typedef core::conformation::Conformation Conformation;
  typedef core::conformation::ConformationOP ConformationOP;
  typedef core::conformation::ConformationCOP ConformationCOP;

public:
  Environment( std::string name );

  virtual ~Environment();

  //@brief register the given mover and, recursively, and submovers yielded by yield_submovers.
  void register_mover( ClaimingMoverOP );

  //@brief register multiple movers (repetitively calls register_mover)
  //@note must be implemented in header for linker
  template< class Iterator >
  void register_movers( Iterator start, Iterator end ){
    for( Iterator it = start; it != end; ++it ){
      register_mover( *it );
    }
  }

  void register_loop_closer( ClaimingMoverOP );

  bool is_registered( ClaimingMoverOP ) const;

  core::pose::Pose start( core::pose::Pose const& );

  core::pose::Pose end( core::pose::Pose const& );

  EnvironmentCAP superenv() const;

  EnvClaimBrokerCOP broker() const { return broker_; }

  void pconf_destruction( ProtectedConformationAP ptr ) const {
    pconfs_.erase( ptr );
  }

  void pconf_creation( ProtectedConformationAP ptr ) const {
    pconfs_.insert( ptr );
  }

private:

  core::conformation::ConformationOP end( ProtectedConformationCOP );

  core::pose::Pose broker( core::pose::Pose const& );

  void remove_nonpermenant_features( core::pose::Pose& );

  void assign_passport( ClaimingMoverOP, core::environment::DofPassportCOP );

  void cancel_passports();

  void remove_chainbreak_variants( core::pose::Pose&, core::Size up_res, core::Size down_res ) const;

  EnvClaimBrokerOP broker_;
  core::kinematics::FoldTreeCOP input_ft_;
  ClaimingMoverOP loop_closer_;
  std::set<ClaimingMoverOP> registered_movers_;

  mutable std::set< ProtectedConformationAP > pconfs_;

}; // end Environment base class

} // moves
} // protocols

#endif //INCLUDED_protocols_environment_Environment_HH
