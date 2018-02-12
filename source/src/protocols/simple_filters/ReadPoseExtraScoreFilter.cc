// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/ReadPoseExtraScoreFilter.cc
/// @author Jack Maguire, jackmaguire1444@gmail.com


//Unit Headers
#include <protocols/simple_filters/ReadPoseExtraScoreFilter.hh>
#include <protocols/simple_filters/ReadPoseExtraScoreFilterCreator.hh>

//Project Headers
#include <basic/datacache/DataMap.hh>
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <protocols/filters/filter_schemas.hh>
#include <protocols/rosetta_scripts/util.hh>

namespace protocols {
namespace simple_filters {

ReadPoseExtraScoreFilter::ReadPoseExtraScoreFilter() :
	Filter( "ReadPoseExtraScoreFilter" ),
	term_name_( "" ),
	threshold_( 0 )
{}

ReadPoseExtraScoreFilter::ReadPoseExtraScoreFilter( ReadPoseExtraScoreFilter const & src ) :
	Filter( "ReadPoseExtraScoreFilter" ),
	term_name_( src.term_name_ ),
	threshold_( src.threshold_ )
{}

ReadPoseExtraScoreFilter::ReadPoseExtraScoreFilter(
	std::string term_name,
	core::Real threshold
) :
	Filter( "ReadPoseExtraScoreFilter" ),
	term_name_( std::move( term_name ) ),
	threshold_( threshold )
{}

ReadPoseExtraScoreFilter::~ReadPoseExtraScoreFilter(){}

void
ReadPoseExtraScoreFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	filters::Filters_map const &,
	moves::Movers_map const &,
	core::pose::Pose const &
) {
	term_name_ = tag->getOption< std::string >( "term_name" );
	threshold_ = tag->getOption< core::Real >( "threshold" );
}

bool
ReadPoseExtraScoreFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const score( compute( pose ) );
	return score <= threshold_;
}

void
ReadPoseExtraScoreFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	out << "Pose Extra Score named " << term_name_ << " is " << compute( pose ) << std::endl;
}

core::Real
ReadPoseExtraScoreFilter::report_sm( core::pose::Pose const & pose ) const {
	return( compute( pose ) );
}

core::Real
ReadPoseExtraScoreFilter::compute( core::pose::Pose const & pose ) const {
	runtime_assert( term_name_.size() );

	// First look to see if the term exists as a double (real).
	core::Real real_result;
	// If not, look to see if it exists as a string.
	std::string string_result;

	if ( core::pose::getPoseExtraScore( pose, term_name_, real_result ) ) {
		return real_result;
	} else if ( core::pose::getPoseExtraScore( pose, term_name_, string_result ) ) {
		core::Real result;
		try {
			result = std::stod( string_result );
		} catch ( ... ){
			utility_exit_with_message( string_result + " could not be converted to a double." );
		}
		return result;
	} else {
		utility_exit_with_message( term_name_ + " is not present in this pose." );
		return 9999.9;
	}
}

void ReadPoseExtraScoreFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "term_name" , xs_string ,
		"Name of the extra score term being searched for"  )
		+ XMLSchemaAttribute::required_attribute( "threshold" , xsct_real ,
		"If that energy is less than or equal to this threshold, returns true." );

	filters::xsd_type_definition_w_attributes(
		xsd,
		class_name(),
		"This filter does no scoring, just looks for extra score terms previously placed in the pose object.",
		attlist
	);
}

}
}
