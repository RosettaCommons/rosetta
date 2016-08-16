// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SilentFilePoseInputStream.fwd.hh
/// @brief
/// @author James Thompson

#ifndef INCLUDED_core_import_pose_pose_stream_SilentFilePoseInputStream_FWD_HH
#define INCLUDED_core_import_pose_pose_stream_SilentFilePoseInputStream_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace import_pose {
namespace pose_stream {

class SilentFilePoseInputStream;
typedef utility::pointer::shared_ptr< SilentFilePoseInputStream > SilentFilePoseInputStreamOP;

} // pose_input_stream
} // import_pose
} // core

#endif
