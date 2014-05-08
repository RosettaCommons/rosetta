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


#ifndef INCLUDED_protocols_environment_claims_JumpClaim_hh
#define INCLUDED_protocols_environment_claims_JumpClaim_hh


// Unit Headers
#include <protocols/environment/claims/JumpClaim.fwd.hh>
#include <protocols/environment/claims/EnvClaim.hh>

// Package Headers
#include <core/environment/LocalPosition.hh>

// Project Headers

// ObjexxFCL Headers

// Utility headers

//// C++ headers

// option key includes


namespace protocols {
namespace environment {
namespace claims {

class JumpClaim : public EnvClaim {
  typedef EnvClaim Parent;
  typedef core::environment::LocalPosition LocalPosition;

public:

  JumpClaim( ClaimingMoverOP owner,
             std::string const& label,
             LocalPosition const& jpos1,
             LocalPosition const& jpos2,
             LocalPosition const& cutp = LocalPosition( "", 0 ) );

  virtual void yield_elements( core::environment::FoldTreeSketch const& fts, JumpElements& elements ) const;

  virtual void yield_elements( core::environment::FoldTreeSketch const& fts, CutElements& elements ) const;

  virtual void yield_elements( ProtectedConformationCOP const&, DOFElements& elements ) const;

  /// @brief set the two atom names to use as the start and beginning of the jump. If set to the empty string
  ///        the broker will use the FoldTree::put_stubs_intra_residue to pick one.
  void set_atoms( std::string const& a1, std::string const& a2 );

  void ctrl_strength( ControlStrength const& str );

  void init_strength( InitializationStrength const& str );

  bool physical() const { return physical_cut_; }

  void physical( bool setting ) { physical_cut_ = setting; }

  std::string const& label() const;

  LocalPosition const& pos1() const;

  LocalPosition const& pos2() const;

  std::string const& atom1() const;

  std::string const& atom2() const;

  virtual EnvClaimOP clone() const;

  virtual std::string str_type() const;

  virtual void show( std::ostream& os ) const;

protected:

  //protected and not public because it can return garbage if the value isn't set.
  LocalPosition const& cut_position() const;

private:

  std::string label_;

  LocalPosition pos1_;
  LocalPosition pos2_;

  LocalPosition cut_;

  std::string atom1_;
  std::string atom2_;

  bool physical_cut_;

  ControlStrength c_str_;
  InitializationStrength i_str_;

}; //JumpClaim

}
}
}

#endif
