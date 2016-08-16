// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/environment/DofUnlock.hh
/// @brief A class that is used to express a mover-specific DoF-unlock on a ProtectedPose. Its destruction expresses a re-locking.
///
/// @author Justin Porter

#ifndef INCLUDED_protocols_environment_DofUnlock_hh
#define INCLUDED_protocols_environment_DofUnlock_hh

// Unit Headers
#include <protocols/environment/DofUnlock.fwd.hh>

// Package headers

// Project headers
#include <core/environment/DofPassport.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/conformation/Conformation.fwd.hh>

// C++ Headers

// ObjexxFCL Headers

namespace protocols {
namespace environment {

class DofUnlock : public utility::pointer::ReferenceCount {

public:
	DofUnlock( core::conformation::Conformation&, core::environment::DofPassportCOP );
	~DofUnlock();

private:
	core::conformation::Conformation& conformation_;
	core::environment::DofPassportCOP pass_;

}; // end DofUnlock base class

} // environment
} // protocols

#endif //INCLUDED_protocols_environment_DofUnlock_hh
