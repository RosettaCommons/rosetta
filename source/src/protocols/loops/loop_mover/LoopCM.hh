// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/LoopCM.hh
/// @brief A loop claiming mover that handles KIC and CCD perturb and refine movers.
/// @author Justin R. Porter

#ifndef INCLUDED_protocols_loops_loop_mover_LoopCM_hh
#define INCLUDED_protocols_loops_loop_mover_LoopCM_hh

#include <protocols/loops/loop_mover/LoopCM.fwd.hh>
#include <protocols/environment/ClientMover.hh>

#include <protocols/environment/claims/EnvClaim.hh>
#include <protocols/loops/loop_mover/IndependentLoopMover.hh>

#include <core/pack/task/residue_selector/ResidueSelector.hh>

#include <basic/datacache/WriteableCacheableMap.fwd.hh>

#ifdef WIN32
	#include <basic/datacache/WriteableCacheableMap.hh>
#endif

// C++ Headers

///////////////////////////////////////////////////////////////////////////////
namespace protocols {
namespace loops {
namespace loop_mover {

class LoopCM: public protocols::environment::ClientMover {
  typedef protocols::environment::ClientMover Parent;

  typedef environment::claims::EnvClaims EnvClaims;
public:
  LoopCM();

  virtual ~LoopCM() {}

  virtual
  void apply( core::pose::Pose& pose );

  virtual
  void parse_my_tag( TagCOP,
                     basic::datacache::DataMap&,
                     Filters_map const&,
                     moves::Movers_map const&,
                     Pose const& );

  virtual
  EnvClaims yield_claims( core::pose::Pose const&,
                          basic::datacache::WriteableCacheableMapOP );

  virtual void passport_updated();

  virtual void initialize( core::pose::Pose& ) {}

  virtual
  std::string get_name() const;

private:
  void build_mover( LoopsOP loops );

  std::string algorithm_;
  std::string style_;

  LoopMoverOP mover_;
  core::pack::task::residue_selector::ResidueSelectorCOP selector_;

}; // class LoopCM

} //namespace loop_mover
} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_loops_loop_mover_LoopCM_HH
