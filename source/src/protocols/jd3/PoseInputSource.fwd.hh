// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/PoseInputSource.fwd.hh
/// @brief  Forward declaration of the %PoseInputSource class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_PoseInputSource_HH
#define INCLUDED_protocols_jd3_PoseInputSource_HH

//unit headers
#include <protocols/jd3/PoseInputSource.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>


namespace protocols {
namespace jd3 {

/// @brief Enum that allows the user to direct the kind of PoseInputter to use
/// to initialize a pose.
enum PoseInputKind {
	pik_unknown,
	pik_pdb_file,
	pik_silent_file
};

enum PoseInputSourceOrigin {
	piso_unknown,
	piso_command_line,
	piso_job_definition_file
};


class PoseInputSource;

typedef utility::pointer::shared_ptr< PoseInputSource > PoseInputSourceOP;
typedef utility::pointer::shared_ptr< PoseInputSource const > PoseInputSourceCOP;

typedef utility::vector1< PoseInputSourceOP > PoseInputSources;

} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_PoseInputSource_HH
