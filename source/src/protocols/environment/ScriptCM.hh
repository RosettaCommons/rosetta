// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/ScriptCM.hh
/// @brief header file for CoMTracker
/// @author Justin R. Porter
/// @author Brian D. Weitzner
/// @author Oliver F. Lange

#ifndef INCLUDED_protocols_environment_ScriptCM_hh
#define INCLUDED_protocols_environment_ScriptCM_hh

// Unit Headers
#include <protocols/environment/ScriptCM.fwd.hh>
#include <protocols/environment/ClaimingMover.hh>

// Package headers#
#include <protocols/moves/MoveMapMover.hh>
#include <protocols/environment/claims/EnvClaim.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/WriteableCacheableMap.fwd.hh>

// C++ Headers

// ObjexxFCL Headers

namespace protocols {
namespace environment {

class ScriptCM : public environment::ClaimingMover {
  typedef environment::claims::EnvClaims EnvClaims;

public:
  ScriptCM();

  virtual
  ~ScriptCM() {};

  virtual
  EnvClaims yield_claims( core::pose::Pose const&,
                          basic::datacache::WriteableCacheableMapOP );

  virtual std::string get_name() const;

  virtual void initialize( core::pose::Pose& pose );

  virtual void apply( core::pose::Pose& );

  virtual void
  parse_my_tag( utility::tag::TagCOP tag,
                basic::datacache::DataMap & data,
                protocols::filters::Filters_map const & filters,
                protocols::moves::Movers_map const & movers,
                core::pose::Pose const & pose );

  std::string const& name() const { return name_; }

  void name( std::string const& name ) { name_ = name; }

  virtual
  moves::MoverOP fresh_instance() const;

  virtual
  moves::MoverOP clone() const;

protected:
  virtual void passport_updated();

  void set_client( moves::MoverOP );

  void add_claim( claims::EnvClaimOP claim );

  moves::MoveMapMoverOP client() const { return client_; }

private:
  std::string name_;
  EnvClaims claim_list_;
  moves::MoveMapMoverOP client_;

}; // end ScriptCM base class

} // environment
} // protocols

#endif //INCLUDED_protocols_environment_ScriptCM_hh
