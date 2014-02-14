// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/abinitio/abscript/StagePreparer.hh
/// @author Justin Porter

#ifndef INCLUDED_protocols_abinitio_abscript_StagePreparer_hh
#define INCLUDED_protocols_abinitio_abscript_StagePreparer_hh

// Unit Headers
#include <protocols/abinitio/abscript/StagePreparer.fwd.hh>

// Package headers
#include <protocols/environment/ClaimingMover.hh>

// Project headers
#include <core/pose/Pose.hh>

// C++ Headers

// ObjexxFCL Headers

namespace protocols {
namespace abinitio {
namespace abscript{

class StagePreparer : public protocols::environment::ClaimingMover {

public:
  virtual void prepare( core::pose::Pose& pose, core::Real progress ) = 0;

}; // end StagePreparer base class

} // abscript
} // abinitio
} // protocols

#endif //INCLUDED_protocols_abinitio_abscript_StagePreparer_hh
