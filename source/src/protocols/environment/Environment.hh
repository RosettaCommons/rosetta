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
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MoverContainer.fwd.hh>

#include <core/environment/DofPassport.fwd.hh>

#include <core/environment/SequenceAnnotation.hh>

#include <protocols/environment/ClaimingMover.fwd.hh>
#include <protocols/environment/ProtectedConformation.fwd.hh>
#include <protocols/environment/EnvClaimBroker.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
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

class Environment : public core::environment::EnvCore, public utility::pointer::enable_shared_from_this< Environment >
{
  typedef core::environment::EnvCore Parent;
  typedef core::environment::SequenceAnnotationCOP SequenceAnnotationCOP;
  typedef core::environment::SequenceAnnotationOP SequenceAnnotationOP;
  typedef core::environment::SequenceAnnotation SequenceAnnotation;

  typedef core::conformation::Conformation Conformation;
  typedef core::conformation::ConformationAP ConformationAP;
  typedef core::conformation::ConformationOP ConformationOP;
  typedef core::conformation::ConformationCOP ConformationCOP;

public:
  Environment( std::string name );
  virtual ~Environment();
  
  //@brief register the given mover and, recursively, and submovers yielded by yield_submovers.
  void register_mover( moves::MoverOP );

  //@brief register multiple movers (repetitively calls register_mover)
  //@note must be implemented in header for linker
  template< class Iterator >
  void register_movers( Iterator start, Iterator end ){
    for( Iterator it = start; it != end; ++it ){
      register_mover( *it );
    }
  }

  bool is_registered( ClaimingMoverOP ) const;

  core::pose::Pose start( core::pose::Pose const& );

  core::pose::Pose end( core::pose::Pose const& );

  EnvironmentCAP superenv() const;

  EnvClaimBrokerCOP broker() const { return broker_; }

  SequenceAnnotationCOP annotations() const { return ann_; }

  bool auto_cut() const { return bAutoCut_; }
  bool inherit_cuts() const { return bInheritCuts_; }
  bool allow_pure_movers() const { return bAllowPureMovers_; }

  void auto_cut( bool );
  void inherit_cuts( bool );
  void allow_pure_movers( bool );

  void pconf_destruction( Conformation * ptr ) const {
    pconfs_.erase( ptr );
  }

  void pconf_creation( Conformation * ptr ) const {
    pconfs_.insert( ptr );
  }

  /// self pointers
  inline EnvironmentCOP get_self_ptr() const { return shared_from_this(); }
  inline EnvironmentOP get_self_ptr() { return shared_from_this(); }
  inline EnvironmentCAP get_self_weak_ptr() const { return EnvironmentCAP( shared_from_this() ); }
  inline EnvironmentAP get_self_weak_ptr() { return EnvironmentAP( shared_from_this() ); }

private:
  
  core::conformation::ConformationOP end( ProtectedConformationCOP );

  core::pose::Pose broker( core::pose::Pose const& );

  void remove_nonpermenant_features( core::pose::Pose& );

  void assign_passport( ClaimingMoverOP, core::environment::DofPassportCOP );

  void cancel_passports();

  void remove_chainbreak_variants( core::pose::Pose&, core::Size up_res, core::Size down_res ) const;

  EnvClaimBrokerOP broker_;
  core::pose::Pose input_pose_;
  std::set<ClaimingMoverOP> registered_movers_;

  SequenceAnnotationOP ann_;

  bool bAutoCut_;
  bool bInheritCuts_;
  bool bAllowPureMovers_;

  mutable std::set< Conformation * > pconfs_;


}; // end Environment base class

} // moves
} // protocols

#endif //INCLUDED_protocols_environment_Environment_HH
