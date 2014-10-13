// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/polar_hydrogens/util.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_modeler_polar_hydrogens_util_HH
#define INCLUDED_protocols_stepwise_modeler_polar_hydrogens_util_HH

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

namespace protocols {
namespace stepwise {
namespace modeler {
namespace polar_hydrogens {

	// move to core::chemical::ResidueType?
	core::Size
	check_if_proton_chi_atom( core::pose::Pose const & pose, core::Size const rsd, core::Size const atomno );

	void
	pack_polar_hydrogens( core::pose::Pose & pose,
												bool allow_virtual_o2prime_hydrogens_ = false );

} //polar_hydrogens
} //modeler
} //stepwise
} //protocols

#endif
