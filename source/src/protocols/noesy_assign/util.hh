// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file CrossPeakList.hh
/// @author Oliver Lange

#ifndef INCLUDED_protocols_noesy_assign_util_hh
#define INCLUDED_protocols_noesy_assign_util_hh


// Unit Headers
// Package Headers


// Project Headers
#include <core/id/NamedAtomID.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>

//Auto Headers
//// C++ headers

namespace protocols {
namespace noesy_assign {

/// @brief WARNING WARNING WARNING THREAD UNSAFE covalent_compliance RELIES ON THREAD-UNSAFE SINGLETON CovalentCompliance
bool covalent_compliance( core::id::NamedAtomID const& atom1, core::id::NamedAtomID const& atom2);
//, core::pose::Pose const& pose, Real dmax );
//core::Real compute_RPF_score( CrossPeakList const& cpl, core::pose::Pose const& pose, core::Real dcut = 5 );

}
}
#endif
