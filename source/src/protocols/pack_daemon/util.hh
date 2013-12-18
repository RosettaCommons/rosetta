// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/pack_daemon/util.hh
/// @brief  Declaration for utility functions for working with the multistate-design classes
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_pack_daemon_util_HH
#define INCLUDED_protocols_pack_daemon_util_HH

// Package headers
#include <protocols/pack_daemon/EntityCorrespondence.fwd.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/ResfileReader.fwd.hh>

// Protocols headers
#include <protocols/genetic_algorithm/EntityRandomizer.fwd.hh>

// C++ headers
#include <istream>
#include <string>

namespace protocols {
namespace pack_daemon {

void
create_entity_resfile_contents(
	std::istream & resfile,
	std::string const & resfile_name,
	core::pack::task::ResfileContentsOP & entity_resfile_contents,
	core::pack::task::PackerTaskOP & entity_task,
	core::Size & num_entities
);

void
initialize_task_from_entity_resfile_and_secondary_resfile(
	core::pose::Pose const & pose,
	EntityCorrespondenceCOP ec,
	core::pack::task::ResfileContents const & entity_resfile_contents,
	core::pack::task::ResfileContents const & secondary_resfile_contents,
	core::pack::task::PackerTaskOP task
);

void
initialize_ga_randomizer_from_entity_task(
	protocols::genetic_algorithm::PositionSpecificRandomizerOP rand,
	core::pack::task::PackerTaskOP entity_task
);

}
}

#endif
