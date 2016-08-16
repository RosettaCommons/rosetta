// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/architects/PoseArchitect.cc
/// @brief Design segments based on a pose
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/denovo_design/architects/PoseArchitect.hh>
#include <protocols/denovo_design/architects/PoseArchitectCreator.hh>

// Protocol headers
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/rosetta_scripts/util.hh>

// Core headers
#include <core/select/residue_selector/TrueResidueSelector.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.architects.PoseArchitect" );

namespace protocols {
namespace denovo_design {
namespace architects {

PoseArchitect::PoseArchitect( std::string const & id_value ):
	DeNovoArchitect( id_value )
{}

PoseArchitect::~PoseArchitect()
{}

DeNovoArchitectOP
PoseArchitect::clone() const
{
	return DeNovoArchitectOP( new PoseArchitect( *this ) );
}

std::string
PoseArchitect::type() const
{
	return architect_name();
}

void
PoseArchitect::parse_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data )
{
	core::select::residue_selector::ResidueSelectorCOP selector_ =
		protocols::rosetta_scripts::parse_residue_selector( tag, data );
}

components::StructureDataOP
PoseArchitect::design( core::pose::Pose const & pose, core::Real & ) const
{
	components::StructureDataOP sd = components::StructureDataFactory::get_instance()->create_from_pose( pose, id() );

	// add template segments
	for ( SegmentNameList::const_iterator s=sd->segments_begin(); s!=sd->segments_end(); ++s ) {
		core::Size const start = sd->segment( *s ).start();
		core::Size const stop = sd->segment( *s ).stop();
		sd->set_template_pose( *s, pose, start, stop );
	}

	return sd;
}

///////////////////////////////////////////////////////////////////////////////
/// Creator
///////////////////////////////////////////////////////////////////////////////
std::string
PoseArchitectCreator::keyname() const
{
	return PoseArchitect::architect_name();
}

DeNovoArchitectOP
PoseArchitectCreator::create_architect( std::string const & id ) const
{
	return DeNovoArchitectOP( new PoseArchitect( id ) );
}

} //architects
} //denovo_design
} //protocols

