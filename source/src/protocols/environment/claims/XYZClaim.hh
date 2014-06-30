// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Oliver Lange


#ifndef INCLUDED_protocols_environment_claims_XYZClaim_hh
#define INCLUDED_protocols_environment_claims_XYZClaim_hh


// Unit Headers
#include <protocols/environment/claims/XYZClaim.fwd.hh>

// Package Headers
#include <core/environment/LocalPosition.hh>
#include <core/environment/FoldTreeSketch.hh>

#include <protocols/environment/claims/EnvClaim.hh>
#include <protocols/environment/claims/BrokerElements.hh>

// Project Headers
#include <core/id/types.hh>
#include <core/id/DOF_ID.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

//#include <basic/options/option_macros.hh>

//// C++ headers
#include <string>
#include <sstream>
#include <utility>

namespace protocols {
namespace environment {
namespace claims {

class XYZClaim : public EnvClaim {
  typedef core::environment::FoldTreeSketch FoldTreeSketch;
  typedef EnvClaim Parent;

public:
  typedef core::environment::LocalPosition LocalPosition;
  typedef core::environment::LocalPositions LocalPositions;

  XYZClaim( ClaimingMoverOP owner,
            utility::tag::TagCOP tag );

  // Initializer for an empty XYZ claim
  XYZClaim( ClaimingMoverOP owner );

  // Initializer for a single backbone angle
  XYZClaim( ClaimingMoverOP owner,
            LocalPosition const& local_pos );

  // Initializer for a contiguous range of residues.
  XYZClaim( ClaimingMoverOP owner,
            std::string const& label,
            std::pair< core::Size, core::Size > const& range );

  virtual void yield_elements( ProtectedConformationCOP const&,
                               DOFElements& elements ) const;

  void add_position( LocalPosition const& p );

  LocalPositions const& positions() const;

  ControlStrength const& ctrl_strength() const;

  void strength( ControlStrength const&, ControlStrength const& );

  ControlStrength const& init_strength() const;

  virtual EnvClaimOP clone() const;

  virtual std::string str_type() const;

  virtual void show( std::ostream& os ) const;

protected:
  virtual
  DOFElement wrap_dof_id( core::id::DOF_ID const& id ) const;

  void build_bond_length_elements( core::Size seqpos,
                                   ProtectedConformationCOP const&,
                                   DOFElements& elements ) const;

  void build_bond_angle_elements( core::Size seqpos,
                                  ProtectedConformationCOP const&,
                                  DOFElements& elements ) const;

  void build_bond_torsion_elements( core::Size seqpos,
                                    ProtectedConformationCOP const&,
                                    DOFElements& elements ) const;

private:
  LocalPositions local_positions_;
  ControlStrength c_str_;
  ControlStrength i_str_;
  bool internal_only_;

}; //class XYZClaim


}
}
}

#endif
