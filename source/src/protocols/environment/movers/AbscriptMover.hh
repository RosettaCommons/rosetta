// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/AbscriptMover.hh
/// @author Justin Porter

#ifndef INCLUDED_protocols_environment_AbscriptMover_hh
#define INCLUDED_protocols_environment_AbscriptMover_hh

// Unit Headers
#include <protocols/environment/movers/AbscriptMover.fwd.hh>

// Package headers
#include <protocols/environment/ClaimingMover.hh>
#include <protocols/environment/movers/StageID.hh>
#include <protocols/environment/movers/AbscriptStageMover.hh>
#include <protocols/environment/movers/StagePreparer.fwd.hh>

#include <protocols/environment/claims/EnvClaim.fwd.hh>

#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/MoverContainer.fwd.hh>

// Project headers
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility Headers
#include <utility/vector0.fwd.hh>

// C++ Headers
#include <set>

// ObjexxFCL Headers

namespace protocols {
namespace environment {

class AbscriptMover : public protocols::environment::ClaimingMover {

public:
  AbscriptMover();

  AbscriptMover( AbscriptMover const& );

  virtual void apply( core::pose::Pose& pose );

  virtual std::string get_name() const;

  // the Abscript mover does not make any claims itself.
  virtual claims::EnvClaims yield_claims( core::pose::Pose& );

  virtual void yield_submovers( std::set< ClaimingMoverOP >& ) const;

  // the Abscript mover does not make any claims, and should never be given
  // initialization rights
  virtual void initialize( Pose& ){ runtime_assert( false ); }

  virtual void
  parse_my_tag( utility::tag::TagCOP const tag,
                basic::datacache::DataMap & data,
                protocols::filters::Filters_map const & filters,
                protocols::moves::Movers_map const & movers,
                core::pose::Pose const & pose );

  virtual
  moves::MoverOP fresh_instance() const;

  virtual
  moves::MoverOP clone() const;

private:

  class StageTracker;

  typedef std::set<ClaimingMoverOP> MoverSet;
  typedef std::map< StageID, MoverSet > IDMoverSetMap;

  void add_default_frags( std::string const& small_frags, std::string const& large_frags );

  void parse_subtags( utility::vector0< utility::tag::TagPtr > const&,
                      protocols::moves::Movers_map const& );

  void register_submover( protocols::moves::MoverOP, StageIDs const&, core::Real weight );

  void register_preparer( protocols::moves::MoverOP, StageIDs const& );

  std::map< Size, Size > calculate_iterations( core::pose::Pose const& );

  StageIDs parse_stage_id( std::string const& ) const;

  std::map< std::string, StageID > id_map_;
  std::map< StageID, AbscriptStageMoverOP > stage_movers_;
  moves::MonteCarloOP mc_;
  ClaimingMoverOP closer_;
}; // end AbscriptMover base class

} // environment
} // protocols

#endif //INCLUDED_protocols_environment_AbscriptMover_hh
