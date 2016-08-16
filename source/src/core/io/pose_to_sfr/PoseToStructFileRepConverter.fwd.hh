// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/pose_to_sfr/PoseToStructFileRepConverter.fwd.hh
/// @brief Forward declarations for class to convert a pose to a StructFileRep.
/// @details This conversion is a first step in PDB or mmCIF output.  It could be useful for other
/// input/output, too.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_core_io_pose_to_sfr_PoseToStructFileRepConverter_fwd_hh
#define INCLUDED_core_io_pose_to_sfr_PoseToStructFileRepConverter_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace io {
namespace pose_to_sfr {

class PoseToStructFileRepConverter;

typedef utility::pointer::shared_ptr< PoseToStructFileRepConverter > PoseToStructFileRepConverterOP;
typedef utility::pointer::shared_ptr< PoseToStructFileRepConverter const > PoseToStructFileRepConverterCOP;

} // namespace silent
} // namespace io
} // namespace core

#endif
