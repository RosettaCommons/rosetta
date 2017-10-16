// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/SequenceDistanceFilter.cc
/// @brief
/// @author Christoffer Norn


//Unit Headers
#include <protocols/simple_filters/SequenceDistanceFilter.hh>
#include <protocols/simple_filters/SequenceDistanceFilterCreator.hh>
#include <utility/tag/Tag.hh>
//Project Headers
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
//#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace simple_filters {

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_filters.SequenceDistance" );

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP SequenceDistanceFilterCreator::create_filter() const { return protocols::filters::FilterOP( new SequenceDistance ); }

// XRW TEMP std::string
// XRW TEMP SequenceDistanceFilterCreator::keyname() const { return "SequenceDistance"; }

//default ctor
SequenceDistance::SequenceDistance() :
	protocols::filters::Filter( "SequenceDistance" ),
	sequence_comment_id_( "" ),
	target_seq_( "" ),
	pose_seq_( "" ),
	threshold_( 100000 )
{
}

SequenceDistance::~SequenceDistance() = default;

void
SequenceDistance::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & pose )
{
	threshold( tag->getOption< core::Size >( "threshold", 8000 ) );

	runtime_assert_string_msg( !(tag->hasOption("sequence_comment_id") && tag->hasOption("native_sequence")), "You can only compare to one sequence: Do you want to load a sequence from comments (save_sequence_name) OR specify the sequence directly in your XML script (native_sequence)" );

	if ( tag->hasOption("sequence_comment_id") ) {
		sequence_comment_id( tag->getOption< std::string >( "sequence_comment_id","" ) );
		target_seq(core::pose::get_all_comments(pose)[ sequence_comment_id() ]);
	} else if ( tag->hasOption("target_sequence") ) {
		target_seq( tag->getOption< std::string >( "target_sequence","" ) );
	}

	pose_seq( pose.sequence() );

	if ( target_seq().length() != pose_seq().length() ) {
		utility_exit_with_message("SequenceDistance filter cannot compare sequences of different length");
	}

}

bool
SequenceDistance::apply( core::pose::Pose const & pose ) const {
	core::Size val = compute( pose );
	TR << "Sequence distance reports " << val << " mutations." << std::endl;
	return( compute( pose ) >= threshold() );
}

void
SequenceDistance::report( std::ostream & o, core::pose::Pose const & pose ) const {
	bool const val = ( compute( pose ) >= threshold() );
	o << "SequenceDistance returns " << val << std::endl;
}

core::Real
SequenceDistance::report_sm( core::pose::Pose const & pose ) const {
	return( core::Real (compute( pose )) );
}

core::Size
SequenceDistance::compute(
	core::pose::Pose const & pose
) const {
	core::Size distance = 0;

	std::string pose_sequence = pose.sequence(); // we already this in pose_sequence, but I need to use pose for something, otherwise compiler warning... Silly.

	for ( core::Size i=1; i <= pose_sequence.length(); i++ ) {
		if ( pose_sequence[i] != target_seq()[i] ) distance++;
	}

	//TR << "The sequence distance is " << distance << std::endl;
	return distance;
}

std::string SequenceDistance::name() const {
	return class_name();
}

std::string SequenceDistance::class_name() {
	return "SequenceDistance";
}

void SequenceDistance::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default(
		"sequence_comment_id", xs_string,
		"Name of sequence to load from pose comments",
		"")
		+ XMLSchemaAttribute::attribute_w_default(
		"target_sequence", xs_string,
		"Specify a sequence to compare to",
		"")
		+ XMLSchemaAttribute::attribute_w_default(
		"threshold", xsct_non_negative_integer,
		"check only whether the comment exists or not, regardless of its content",
		"1000");

	protocols::filters::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"Measures the distance between a sequence stored in pose comments"
		"and the current sequence of pose. Useful together with SaveSequenceToCommentsMover.",
		attlist );
}

std::string SequenceDistanceFilterCreator::keyname() const {
	return SequenceDistance::class_name();
}

protocols::filters::FilterOP
SequenceDistanceFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new SequenceDistance );
}

void SequenceDistanceFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SequenceDistance::provide_xml_schema( xsd );
}


}
}
