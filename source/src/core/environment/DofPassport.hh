// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/environment/Enironment.hh
/// @brief An environment that automatically distributes rights to shared degrees of freedom (e.g. fold tree)
///
/// @author Justin Porter

#ifndef INCLUDED_core_environment_DofPassport_hh
#define INCLUDED_core_environment_DofPassport_hh

// Unit Headers
#include <core/environment/DofPassport.fwd.hh>

// Package headers
#include <core/environment/EnvCore.fwd.hh>

// Project headers
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>

#include <core/conformation/Conformation.fwd.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/JumpID.hh>

#include <utility/vector1.fwd.hh>

// C++ Headers
#include <ostream>
#include <set>

// ObjexxFCL Headers

namespace core {
namespace environment {


class DofPassport : public utility::pointer::ReferenceCount {

  //make EnvCore a friend so that we can call the constructor of DofPassport
  friend class core::environment::EnvCore;

public:
  virtual ~DofPassport();

  void show( std::ostream& ) const;

  core::kinematics::MoveMapOP render_movemap() const;

  // typedef to make iteration through availiable dofs easy for client movers
  typedef std::set< core::id::DOF_ID >::const_iterator const_iterator;

  std::set< core::id::DOF_ID >::const_iterator begin() const;

  std::set< core::id::DOF_ID >::const_iterator end() const { return accessible_dofs_.end(); }

  ///@brief configure passport to allow access to a bond length or angle dof_id
  ///@param id the DOF_ID to allow access to
  void add_dof_access( id::DOF_ID const& id );

  void revoke_all_access();

  bool has_jump_access( int jump_num ) const;

  utility::vector1< int > active_jumps() const;

  //@brief Check access in this passport for an id
  bool dof_access( core::id::DOF_ID const& id ) const;

  //@brief Check access in this passport for both correct environment and correct id.
  bool dof_access( EnvCore const&, core::id::DOF_ID const& id ) const;

  std::string const& mover() const;

  core::Size const& env_id() const { return env_id_; }

  conformation::ConformationCOP reference_conformation() const;

  void reference_conformation( conformation::ConformationCOP conf );

private:

  DofPassport( std::string const& mover, Size env_id);

  DofPassport( DofPassport const& );

  //General code used by most (all?) access checks.
  bool access_check( EnvCore const& env, bool type_specific_check ) const;

  std::string mover_;
  core::Size env_id_;

  // Capable of storing all movable dofs.
  std::set< core::id::DOF_ID > accessible_dofs_;

  conformation::ConformationCOP conf_;

}; // end DofPassport base class

extern std::ostream& operator<<( std::ostream&, DofPassport const& );
} // environment
} // core

#endif //INCLUDED_core_environment_DofPassport_HH
