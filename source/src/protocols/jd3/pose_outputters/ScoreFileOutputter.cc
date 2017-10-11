// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/pose_outputters/ScoreFileOutputter.cc
/// @brief  Definition of the %ScoreFileOutputter class's methods
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

//unit headers
#include <protocols/jd3/pose_outputters/ScoreFileOutputter.hh>
#include <protocols/jd3/pose_outputters/ScoreFileOutputterCreator.hh>

// package headers
#include <protocols/jd3/InnerLarvalJob.hh>
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/pose_outputters/PoseOutputterFactory.hh>

#include <core/types.hh>
#include <core/io/raw_data/ScoreFileData.hh>
#include <core/io/raw_data/ScoreMap.hh>

#include <basic/datacache/ConstDataMap.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
//#include <basic/options/keys/jd3.OptionKeys.gen.hh>

#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/options/keys/OptionKeyList.hh>

namespace protocols {
namespace jd3 {
namespace pose_outputters {

ScoreFileOutputter::ScoreFileOutputter() {}

ScoreFileOutputter::~ScoreFileOutputter() { flush(); }

std::string
ScoreFileOutputter::outputter_for_job(
	utility::tag::TagCOP output_tag,
	utility::options::OptionCollection const & job_options,
	InnerLarvalJob const & job
) const
{
	utility::file::FileName scorefile_name_for_job = filename_for_job( output_tag, job_options, job );
	return scorefile_name_for_job.name();
}

bool
ScoreFileOutputter::outputter_specified_by_command_line()
{
	return ! basic::options::option[ basic::options::OptionKeys::run::no_scorefile ];
}

void
ScoreFileOutputter::write_output_pose(
	LarvalJob const & job,
	JobOutputIndex const & output_index,
	utility::options::OptionCollection const & options,
	core::pose::Pose const & pose
)
{
	if ( scorefile_name_.name() == ""  ) {
		utility::tag::TagCOP job_output_tag;
		if ( job.inner_job()->jobdef_tag() ) {
			utility::tag::TagCOP job_tags = job.inner_job()->jobdef_tag();
			if ( job_tags->hasTag( "SecondaryOutput" ) ) {
				job_output_tag = job_tags->getTag( "SecondaryOutput" );
			}
		}
		scorefile_name_ = filename_for_job( job_output_tag, options, *job.inner_job() );
	}

	core::io::raw_data::ScoreFileData sfd( scorefile_name_ );

	std::map < std::string, core::Real > score_map;
	std::map < std::string, std::string > string_map;
	core::io::raw_data::ScoreMap::score_map_from_scored_pose( score_map, pose );
	core::io::raw_data::ScoreMap::add_arbitrary_score_data_from_pose( pose, score_map );
	core::io::raw_data::ScoreMap::add_arbitrary_string_data_from_pose( pose, string_map );

	// TO DO:
	// buffer the score file output and then write out the contents
	// in "flush," which is called only sporadically

	sfd.write_pose( pose, score_map, job.status_prefix() + job.job_tag_with_index_suffix( output_index ) + job.status_suffix(), string_map );

}


void
ScoreFileOutputter::flush()
{
	// TO DO: once silent-file data is being buffered in this class,
	// then the stored lines should be written out to disk.
	// currently, a noop.
}

std::string
ScoreFileOutputter::class_key() const
{
	return keyname();
}

std::string
ScoreFileOutputter::keyname()
{
	return "ScoreFile";
}

void ScoreFileOutputter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	//XMLSchemaRestriction string_w_one_dollarsign;
	//string_w_one_dollarsign.name( "string_w_one_dollarsign" );
	//string_w_one_dollarsign.base_type( xs_string );
	//// this reads "anything but a dollar sign or nothing followed by a dollar sign followed by anything but a dollar sign or nothing
	//string_w_one_dollarsign.add_restriction( xsr_pattern, "[^$]*$[^$]*" );
	//xsd.add_top_level_element( string_w_one_dollarsign );

	AttributeList output_scorefile_attributes;
	output_scorefile_attributes
		+ XMLSchemaAttribute::required_attribute( "filename", xs_string , "XRW TO DO" ) // this stops being required if we instead accept a filename-pattern
		+ XMLSchemaAttribute( "path", xs_string , "XRW TO DO" );
	XMLSchemaComplexTypeGenerator output_scorefile;
	output_scorefile.element_name( keyname() )
		.description( "XRW TO DO" )
		.complex_type_naming_func( & PoseOutputterFactory::complex_type_name_for_secondary_pose_outputter )
		.add_attributes( output_scorefile_attributes )
		.write_complex_type_to_schema( xsd );
}

void ScoreFileOutputter::list_options_read( utility::options::OptionKeyList & read_options )
{
	using namespace basic::options::OptionKeys;
	read_options
		+ run::no_scorefile
		+ out::file::scorefile
		+ out::path::score
		+ out::path::all
		+ out::prefix
		+ out::suffix;
}


utility::file::FileName
ScoreFileOutputter::filename_for_job(
	utility::tag::TagCOP output_tag,
	utility::options::OptionCollection const & job_options,
	InnerLarvalJob const & /*job*/
) const
{
	using namespace basic::options::OptionKeys;

	utility::file::FileName scorefile_name;

	if ( output_tag && output_tag->hasTag( keyname() ) ) {
		utility::tag::TagCOP score_tag = output_tag->getTag( keyname() );
		if ( score_tag->hasOption( "filename" ) ) {
			scorefile_name = score_tag->getOption< std::string >( "filename" );
			if ( score_tag->hasOption( "path" ) ) {
				utility::file::FileName scorefile_path( score_tag->getOption< std::string >( "path" ));
				scorefile_name.path( scorefile_path.path() );
			}
			return scorefile_name;
		}
	}

	if ( job_options[ run::no_scorefile ] ) {
		scorefile_name = "(none)";
		return scorefile_name;
	}

	if ( job_options[ out::file::scorefile ].user() ) {
		scorefile_name = job_options[ out::file::scorefile ]();
		if ( job_options[ out::path::score ].user() ) {
			scorefile_name.path( job_options[ out::path::score ]().path() );
		} else if ( ! scorefile_name.absolute() && job_options[ out::path::all ].user() ) {
			scorefile_name.path( job_options[ out::path::all ]().path() + "/" + scorefile_name.path() );
		}
	} else {
		scorefile_name = "score"; // default name "score.sc"
		std::ostringstream oss;
		//prefix, suffix
		oss << job_options[ out::prefix ]() << scorefile_name.base()
			<< job_options[ out::suffix ]();
		scorefile_name.base( oss.str() );
		//path
		if ( job_options[ out::path::score ].user() ) {
			scorefile_name.path( job_options[ out::path::score ]().path() );
			scorefile_name.vol( job_options[ out::path::score ]().vol() );
		} else if ( job_options[ out::path::all ].user() ) {
			scorefile_name.path( job_options[ out::path::all ]().path() + "/" + scorefile_name.path() );
		}

		// Deviation from JD2's score file naming scheme:
		// Do not put .fasc at the end of a score file if the out::file::fullatom flag
		// is on the command line; instead, if the user wants to
		// say that their score file should should be named with the .fasc extension
		// then they can give a name for the file explicitly.
		scorefile_name.ext( ".sc" );
	}

	return scorefile_name;
}


SecondaryPoseOutputterOP ScoreFileOutputterCreator::create_outputter() const {
	return SecondaryPoseOutputterOP( new ScoreFileOutputter );
}

std::string ScoreFileOutputterCreator::keyname() const
{
	return ScoreFileOutputter::keyname();
}

void ScoreFileOutputterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ScoreFileOutputter::provide_xml_schema( xsd );
}

void ScoreFileOutputterCreator::list_options_read( utility::options::OptionKeyList & read_options ) const
{
	ScoreFileOutputter::list_options_read( read_options );
}

bool ScoreFileOutputterCreator::outputter_specified_by_command_line() const {
	return ScoreFileOutputter::outputter_specified_by_command_line();
}

} // namespace pose_outputters
} // namespace jd3
} // namespace protocols
