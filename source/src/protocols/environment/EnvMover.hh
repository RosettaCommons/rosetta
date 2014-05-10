// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/EnvMover.hh
/// @author Justin Porter

#ifndef INCLUDED_protocols_environment_EnvMover_hh
#define INCLUDED_protocols_environment_EnvMover_hh

// Unit Headers
#include <protocols/environment/EnvMover.fwd.hh>

// Package headers
#include <protocols/moves/Mover.hh>
#include <protocols/environment/Environment.fwd.hh>

// Project headers
#include <protocols/moves/MoverContainer.hh>

// C++ Headers

// ObjexxFCL Headers

namespace protocols {
namespace environment {

class EnvMover : public moves::Mover {

public:
  EnvMover();

  EnvMover( EnvironmentOP env );

  virtual void
  parse_my_tag( utility::tag::TagCOP tag,
               basic::datacache::DataMap & data,
               protocols::filters::Filters_map const & filters,
               protocols::moves::Movers_map const & movers,
               core::pose::Pose const& pose );

  virtual ~EnvMover();

  virtual void apply( Pose& pose );

  virtual std::string get_name() const;

  virtual moves::MoverOP clone() const;

private:
  void parse_subtag( utility::tag::TagCOP tag,
                     basic::datacache::DataMap & data,
                     protocols::filters::Filters_map const & filters,
                     protocols::moves::Movers_map const & movers,
                     core::pose::Pose const& pose );

  EnvironmentOP env_;
  moves::SequenceMover movers_;

}; // end EnvMover base class

} // environment
} // protocols

#endif //INCLUDED_protocols_environment_EnvMover_hh
