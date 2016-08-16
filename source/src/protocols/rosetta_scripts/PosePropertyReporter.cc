// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/rosetta_scripts/PosePropertyReporter.hh
/// @brief  Base class for pose selectors used by MultiplePoseMover
/// @author Luki Goldschmidt <lugo@uw.edu>

#ifndef INCLUDED_protocols_rosetta_scripts_PosePropertyReporter_CC
#define INCLUDED_protocols_rosetta_scripts_PosePropertyReporter_CC

// Unit Headers
#include <protocols/rosetta_scripts/PosePropertyReporter.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/excn/Exceptions.hh>

// C++ Headers
#include <string>

static THREAD_LOCAL basic::Tracer TR( "protocols.rosetta_scripts.PosePropertyReporter" );

namespace protocols {
namespace rosetta_scripts {

PosePropertyReporter::PosePropertyReporter()
{
}

PosePropertyReporter::~PosePropertyReporter()
{
}

core::Real PosePropertyReporter::report_property( core::pose::Pose & ) const
{
	throw utility::excn::EXCN_RosettaScriptsOption("Pose property reporter " + get_name() + " cannot report property for a single pose -- not implemented" );
}

core::Real PosePropertyReporter::report_property( core::pose::Pose & , core::pose::Pose & ) const
{
	throw utility::excn::EXCN_RosettaScriptsOption("Pose property reporter " + get_name() + " cannot compare two poses -- not implemented" );
}


void PosePropertyReporter::parse_my_tag(
	utility::tag::TagCOP,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {
	TR << "***WARNING!!!! WARNING!!!*** parse_my_tag has been invoked for this PosePropertyReporter but it hasn't been defined. Are you sure this is appropriate?" << std::endl;
}

} // rosetta_scripts
} // protocols

#endif //INCLUDED_protocols_rosetta_scripts_PosePropertyReporter_CC
