// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Sarel Fleishman
#include <protocols/protein_interface_design/filters/RelativeSegmentFilter.hh>
#include <protocols/protein_interface_design/filters/RelativeSegmentFilterCreator.hh>


// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/import_pose/import_pose.hh>

namespace protocols {
namespace protein_interface_design {
namespace filters {

static THREAD_LOCAL basic::Tracer TR( "protocols.protein_interface_design.filters.RelativeSegmentFilter" );

using core::pose::Pose;

protocols::filters::FilterOP
RelativeSegmentFilterCreator::create_filter() const { return protocols::filters::FilterOP( new RelativeSegmentFilter ); }

std::string
RelativeSegmentFilterCreator::keyname() const { return "RelativeSegment"; }

bool
RelativeSegmentFilter::apply( Pose const & pose ) const {
	core::pose::Pose source_pose;
	core::import_pose::pose_from_file( source_pose, source_pose_ , core::import_pose::PDB_file);

	core::Size const nearest_to_from( protocols::rosetta_scripts::find_nearest_res( pose, source_pose, start_res() ) );
	core::Size const nearest_to_to( protocols::rosetta_scripts::find_nearest_res( pose, source_pose, stop_res() ) );

	if ( nearest_to_from == 0 || nearest_to_to == 0 ) {
		TR<<"nearest_to_from: "<<nearest_to_from<<" nearest_to_to: "<<nearest_to_to<<". Failing"<<std::endl;
		return false;
	}
	runtime_assert( nearest_to_to > nearest_to_from );
	TR<<"Residues on source pose: ";
	for ( core::Size i = nearest_to_from; i <= nearest_to_to; ++i ) {
		TR<<i<<",";
	}
	TR<<std::endl;
	return true;
}

void
RelativeSegmentFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & pose )
{
	source_pose( tag->getOption< std::string >( "source_pose" ) );
	start_res( tag->getOption< core::Size >( "start_res" ) );
	stop_res( tag->getOption< core::Size >( "stop_res" ));
	runtime_assert( stop_res() > start_res() );
	runtime_assert( start_res() > 0 );
	runtime_assert( stop_res() <= pose.size() );
	TR<<"source_pose: "<<source_pose()<<" start_res: "<<start_res()<<" stop res: "<<stop_res()<<std::endl;
}

void
RelativeSegmentFilter::report( std::ostream &, core::pose::Pose const & ) const {
}

core::Real
RelativeSegmentFilter::report_sm( core::pose::Pose const & ) const {
	return( 0.0 );
}

RelativeSegmentFilter::~RelativeSegmentFilter() {}

}
} // protein_interface_design
} // devel
