// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/util.hh
/// @brief Useful functions for rotamer sampler.
/// @detailed
/// @author Fang-Chieh Chou

#ifndef INCLUDED_protocols_rotamer_sampler_util_HH
#define INCLUDED_protocols_rotamer_sampler_util_HH

#include <core/types.hh>
#include <utility/vector1.fwd.hh>

namespace protocols {
namespace rotamer_sampler {

/// @brief Add torsions from center value, max range and bin size to a
//         existing torsion_list
void add_values_from_center(
	utility::vector1<core::Real> & torsions,
	core::Real const center,
	core::Real const max_range,
	core::Real const bin_size
);

}
}

#endif
