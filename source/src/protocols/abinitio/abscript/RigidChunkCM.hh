// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/abinitio/abscript/RigidChunkCM.hh
/// @author Justin Porter

#ifndef INCLUDED_protocols_abinitio_abscript_RigidChunkCM_hh
#define INCLUDED_protocols_abinitio_abscript_RigidChunkCM_hh

// Unit Headers
#include <protocols/abinitio/abscript/RigidChunkCM.fwd.hh>
#include <protocols/environment/ClaimingMover.hh>

// Package headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.tmpl.hh>


#ifdef WIN32
  #include <basic/datacache/WriteableCacheableMap.hh>
#endif

// Project headers

// C++ Headers

// ObjexxFCL Headers

namespace protocols {
namespace abinitio {
namespace abscript {

class RigidChunkCM : public protocols::environment::ClaimingMover {
  typedef ClaimingMover Parent;
  typedef environment::claims::EnvClaims EnvClaims;

public:
  RigidChunkCM();

  RigidChunkCM( std::string const& label,
                loops::Loops const& rigid_core,
                core::pose::Pose const& template_pose );

  virtual ~RigidChunkCM() {};

  virtual EnvClaims yield_claims( core::pose::Pose const&,
                                  basic::datacache::WriteableCacheableMapOP );

  virtual std::string get_name() const;

  virtual void initialize( Pose& pose );

  virtual void apply( core::pose::Pose& );

  std::string const& label() const { return label_; }

  void label( std::string const& label ) { label_ = label; }

  virtual void
  parse_my_tag( utility::tag::TagCOP tag,
               basic::datacache::DataMap & data,
               protocols::filters::Filters_map const & filters,
               protocols::moves::Movers_map const & movers,
               core::pose::Pose const & pose );

  virtual
  moves::MoverOP clone() const;

  loops::Loops select_parts( loops::Loops const& rigid_core,
                             core::Size random_grow_loops_by );

private:
  std::string label_;
  EnvClaims claims_;
  loops::Loops rigid_core_;
  core::pose::PoseCOP template_;

}; // end RigidChunkCM base class

} // abscript
} // abinitio
} // protocols

#endif //INCLUDED_protocols_abinitio_abscript_RigidChunkCM_hh
