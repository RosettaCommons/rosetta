// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file util.hh
/// @brief
/// @details
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_rna_util_HH
#define INCLUDED_protocols_rna_util_HH

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

namespace protocols {
namespace rna {

void
setup_base_pair_constraints(
	core::pose::Pose & pose,
	utility::vector1< std::pair< core::Size, core::Size > > const &  pairings,
	core::Real const scale_factor = 1.0,
	bool const use_flat_harmonic = false );

void
get_base_pairing_list(
	core::pose::Pose & pose,
	utility::vector1< std::pair< core::Size, core::Size> > & base_pairing_list );

void
assert_phosphate_nomenclature_matches_mini( core::pose::Pose const & pose);

} //rna
} //protocols

#endif