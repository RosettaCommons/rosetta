// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file BrokerElements
/// @brief
/// @author Justin Porter

#ifndef INCLUDED_protocols_environment_claims_BrokerElements_hh
#define INCLUDED_protocols_environment_claims_BrokerElements_hh

// Unit Headers
// #include <protocols/environment/claims/BrokerElements.fwd.hh>

// Package Headers
#include <protocols/environment/ClaimingMover.fwd.hh>

// Project Headers
#include <core/environment/LocalPosition.hh>

// ObjexxFCL Headers

// Utility headers

//// C++ headers
#include <string>

namespace protocols {
namespace environment {
namespace claims {

  // ControlStrengths indicate the level of control over the dof desired by the mover.
  // DOES_NOT_CONTROL is the lowest level, used for protocols that exert no control
  //   over the claimed dof (perhaps, for example, because they only require that it exists).
  // CAN_CONTROL is also a low, uncommon level indicating the mover would like to control
  //   (i.e. samples) the dof, but needn't. Fragment Insertion is often reliant on claims
  //   of this type, because other movers can take priority (i.e. if a region is fixed).
  // MUST_CONTROL asserts a need to sample the degree of freedom (e.g. jumps in docking).
  //   This level allows the mover, however, to share the dof with other movers.
  // EXCLUSIVE is as MUST_CONTROL, except that no other mover can be granted access to the dof.
  enum ControlStrength  {
    DOES_NOT_CONTROL = 0,
    CAN_CONTROL = 1,
    MUST_CONTROL,
    EXCLUSIVE
  };

  enum InitializationStrength {
    DOES_NOT_INITIALIZE = 0,
    CAN_INITIALIZE = 1,
    MUST_INITIALIZE
  };

  struct ResidueElement {
    ResidueElement() : label() {}
    std::string label;
    static std::string const type;
  };

  struct JumpElement {
    JumpElement() : label(), p1(), p2(), atom1(), atom2(), has_physical_cut(false) {}
    std::string label;
    core::environment::LocalPosition p1;
    core::environment::LocalPosition p2;
    std::string atom1;
    std::string atom2;
    bool has_physical_cut;
    static std::string const type;
  };

  struct CutElement {
    CutElement() : p(), physical( true ) {}
    core::environment::LocalPosition p;
    bool physical; //indicates if this cut will be closed later/should be scored as a chainbreak
    static std::string const type;
  };

  struct CutBiasElement {
    CutBiasElement() : p(), bias( 0.0 ) {}
    core::environment::LocalPosition p;
    core::Real bias;
    static std::string const type;
  };

  struct RTElement {
    // no initializer because then we'd have to include the whole ClaimingMover header.
    ClaimingMoverOP owner;
    std::string label;
    ControlStrength c_str;
    InitializationStrength i_str;
    static std::string const type;
  };

  struct TorsionElement {
    // no initializer because then we'd have to include the whole ClaimingMover header.
    ClaimingMoverOP owner;
    core::environment::LocalPosition p;
    ControlStrength c_str;
    InitializationStrength i_str;
    static std::string const type;
  };

  typedef utility::vector1< claims::ResidueElement > ResidueElements;
  typedef utility::vector1< claims::JumpElement > JumpElements;
  typedef utility::vector1< claims::CutElement > CutElements;
  typedef utility::vector1< claims::CutBiasElement > CutBiasElements;
  typedef utility::vector1< claims::RTElement > RTElements;
  typedef utility::vector1< claims::TorsionElement > TorsionElements;

  typedef utility::vector1< claims::ControlStrength > ControlStrengths;
  typedef utility::vector1< claims::ControlStrength > InitializationStrengths;

}
}
}

#endif
