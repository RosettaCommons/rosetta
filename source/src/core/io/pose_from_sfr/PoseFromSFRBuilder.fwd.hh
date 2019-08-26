// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/PoseFromSFRBuilder.fwd.hh
/// @brief FWD header for PoseFromSFRBuilder
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_io_PoseFromSFRBuilder_fwd_hh
#define INCLUDED_core_io_PoseFromSFRBuilder_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace core {
namespace io {
namespace pose_from_sfr {


class PoseFromSFRBuilder;

typedef utility::pointer::shared_ptr< PoseFromSFRBuilder > PoseFromSFRBuilderOP;
typedef utility::pointer::shared_ptr< PoseFromSFRBuilder const > PoseFromSFRBuilderCOP;

} //pose_from_sfr
} //core
} //io

#endif //INCLUDED_core_io_PoseFromSFRBuilder_fwd_hh
