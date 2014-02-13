// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/EnvEndMover.hh
/// @author Justin Porter

#ifndef INCLUDED_protocols_environment_EnvEndMover_hh
#define INCLUDED_protocols_environment_EnvEndMover_hh

// Unit Headers
#include <protocols/environment/EnvEndMover.fwd.hh>

// Package headers
#include <protocols/moves/Mover.hh>

#include <protocols/environment/Environment.fwd.hh>

// Project headers

// C++ Headers

// ObjexxFCL Headers

namespace protocols {
namespace environment {

class EnvEndMover : public moves::Mover {

public:
  EnvEndMover( EnvironmentOP env );

  virtual ~EnvEndMover();

  virtual void apply( Pose& pose );

  virtual std::string get_name() const;

  virtual moves::MoverOP clone() const;

private:
  EnvironmentOP env_;

}; // end EnvEndMover base class

} // environment
} // protocols

#endif //INCLUDED_protocols_environment_EnvEndMover_hh
