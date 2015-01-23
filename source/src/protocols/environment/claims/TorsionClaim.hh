// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/environment/claims/TorsionClaim.hh
/// @brief An EnvClaim object for claiming arbitrary regions and types of torsional angles within a pose.
/// @author Justin R. Porter


#ifndef INCLUDED_protocols_environment_claims_TorsionClaim_hh
#define INCLUDED_protocols_environment_claims_TorsionClaim_hh


// Unit Headers
#include <protocols/environment/claims/TorsionClaim.fwd.hh>

// Package Headers
#include <protocols/environment/ProtectedConformation.fwd.hh>

#include <protocols/environment/claims/EnvClaim.hh>
#include <protocols/environment/claims/BrokerElements.hh>

#include <core/environment/LocalPosition.hh>
#include <core/environment/FoldTreeSketch.hh>

#include <core/conformation/Conformation.hh>

// Project Headers
#include <core/pack/task/residue_selector/ResidueSelector.hh>

#include <core/id/types.hh>

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

class TorsionClaim : public EnvClaim {
  typedef core::environment::FoldTreeSketch FoldTreeSketch;
  typedef EnvClaim Parent;

  typedef core::pack::task::residue_selector::ResidueSelectorCOP ResidueSelectorCOP;

public:
  typedef core::environment::LocalPosition LocalPosition;
  typedef core::environment::LocalPositions LocalPositions;

  TorsionClaim( ClientMoverOP owner,
                utility::tag::TagCOP tag,
                basic::datacache::DataMap& );

  TorsionClaim( ClientMoverOP owner,
                core::pack::task::residue_selector::ResidueSelectorCOP );

  // Initializer for a single backbone angle
  TorsionClaim( ClientMoverOP owner,
                LocalPosition const& local_pos );

  // Initializer for a contiguous range of backbone angles.
  TorsionClaim( ClientMoverOP owner,
                std::string const& label,
                std::pair< core::Size, core::Size > const& range );

  TorsionClaim( ClientMoverOP owner,
                LocalPositions const& positions );

  virtual void yield_elements( core::pose::Pose const&, DOFElements& elements ) const;

  ControlStrength const& ctrl_strength() const;

  void claim_sidechain( bool in ) { claim_sidechain_ = in; }

  bool claim_sidechain() const { return claim_sidechain_; }

  void claim_backbone( bool in ) { claim_backbone_ = in; }

  bool claim_backbone() const { return claim_backbone_; }

  /// @brief set the initialization and control strength of the TorsionClaim.
  void strength( ControlStrength const& control_strength,
                 ControlStrength const& initialization_strength );

  ControlStrength const& init_strength() const;

  virtual EnvClaimOP clone() const;

  virtual std::string type() const;

  virtual ResidueSelectorCOP selector() const { return selector_; }

  virtual void show( std::ostream& os ) const;

protected:
  virtual
  DOFElement wrap_dof_id( core::id::DOF_ID const& id ) const;

  void insert_dof_element( core::conformation::Conformation const& conf,
                           DOFElements& elements,
                           core::Size seqpos,
                           core::id::TorsionType type,
                           core::Size torsion_number) const;


private:
  ResidueSelectorCOP selector_;

  ControlStrength c_str_;
  ControlStrength i_str_;

  bool claim_sidechain_;
  bool claim_backbone_;

}; //class TorsionClaim


}
}
}

#endif
