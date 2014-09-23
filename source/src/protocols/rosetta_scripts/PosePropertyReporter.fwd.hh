// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/rosetta_scripts/PosePropertyReporter.hh
/// @brief  Forward declarations for PosePropertyReporter
/// @author Luki Goldschmidt <lugo@uw.edu>

#ifndef INCLUDED_protocols_rosetta_scripts_PosePropertyReporter_fwd_hh
#define INCLUDED_protocols_rosetta_scripts_PosePropertyReporter_fwd_hh

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace rosetta_scripts {

class PosePropertyReporter;
typedef utility::pointer::shared_ptr< PosePropertyReporter > PosePropertyReporterOP;
typedef utility::pointer::shared_ptr< PosePropertyReporter const > PosePropertyReporterCOP;

} // rosetta_scripts
} // protocols

#endif //INCLUDED_protocols_rosetta_scripts_PosePropertyReporter_fwd_hh
