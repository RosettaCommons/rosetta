// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/SilentFilePoseInputter.cc
/// @brief
/// @author Andy Watkins (amw579@stanford.edu)

///Unit headers
#include <protocols/jd3/pose_inputters/SilentFilePoseInputter.hh>
#include <protocols/jd3/pose_inputters/SilentFilePoseInputterCreator.hh>

// Package headers
#include <protocols/jd3/PoseInputSource.hh>
#include <protocols/jd3/pose_inputters/PoseInputterFactory.hh>
#include <protocols/jd3/pose_inputters/pose_inputter_schemas.hh>

///Project headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>

///Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/keys/OptionKeyList.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

///C++ headers
#include <string>

// External headers
#include <boost/algorithm/string/predicate.hpp>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/pose/symmetry/util.hh>

static THREAD_LOCAL basic::Tracer tr( "protocols.jd3.SilentFilePoseInputter" );

namespace protocols {
namespace jd3 {
namespace pose_inputters {

SilentFilePoseInputter::SilentFilePoseInputter()
{
	tr.Debug << "Instantiate SilentFilePoseInputter" << std::endl;
}

SilentFilePoseInputter::~SilentFilePoseInputter() {}

bool SilentFilePoseInputter::job_available_on_command_line() const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	return option[ in::file::silent ].user();
}

PoseInputSources SilentFilePoseInputter::pose_input_sources_from_command_line()
{
	using namespace core::io::silent;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	initialize_sfd_from_options_and_tag( basic::options::option, utility::tag::TagOP() );

	PoseInputSources input_sources;
	for ( auto iter : *sfd_ ) {
		std::string const & tag = iter->decoy_tag();
		if ( boost::starts_with( tag, "W_" ) && option[ in::file::skip_failed_simulations ] ) {
			continue;
		}

		PoseInputSourceOP input_source( new PoseInputSource );
		input_source->origin( keyname() );
		input_source->input_tag( tag );
		tr << "Tag: " << tag << std::endl;
		input_sources.push_back( input_source );
	}
	return input_sources;
}

PoseInputSources
SilentFilePoseInputter::pose_input_sources_from_tag(
	utility::options::OptionCollection const & opts,
	utility::tag::TagCOP tag
)
{
	using namespace core::io::silent;
	using namespace basic::options::OptionKeys;

	initialize_sfd_from_options_and_tag( opts, tag );

	PoseInputSources input_sources;
	for ( auto iter : *sfd_ ) {
		std::string const & decoy_tag = iter->decoy_tag();
		if ( boost::starts_with( decoy_tag, "W_" ) && ( opts[ in::file::skip_failed_simulations]
				|| tag->getOption< bool >( "skip_failed_simulations", false ) ) ) {
			continue;
		}

		PoseInputSourceOP input_source( new PoseInputSource );
		input_source->origin( keyname() );
		input_source->input_tag( decoy_tag );
		input_sources.push_back( input_source );
	}
	return input_sources;
}


/// @details This function will first see if the pose already exists in the Job.
/// If not, it will read it into the pose reference, and hand a COP cloned from
/// that pose to the Job. If the pose pre-exists it just copies the COP's pose
/// into it.
core::pose::PoseOP
SilentFilePoseInputter::pose_from_input_source(
	PoseInputSource const & input_source,
	utility::options::OptionCollection const & options,
	utility::tag::TagCOP tag // possibly null-pointing tag pointer
)
{
	if ( ! sfd_ ) { initialize_sfd_from_options_and_tag( options, tag ); }

	debug_assert( sfd_->has_tag( input_source.input_tag() ));

	using namespace core::io::silent;
	SilentStruct const & silent_struct( sfd_->get_structure( input_source.input_tag() ) );
	core::pose::PoseOP pose( new core::pose::Pose );
	silent_struct.fill_pose( *pose );
	return pose;
}

std::string SilentFilePoseInputter::keyname() { return "Silent"; }

/// @brief returns the schema for the PDB element used in a job-definition file
/// including all options that govern how a PDB is loaded.
void
SilentFilePoseInputter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	typedef XMLSchemaAttribute Attr;
	AttributeList attributes;
	attributes
		+ Attr::required_attribute( "silent_files", xs_string, "Comma-separated list of silent files to use" )
		+ Attr( "tags", xs_string, "Comma-separated list of tags specifying the subset of Poses that should be"
		" processed from the input silent file(s). If neither this attribute, nor the 'tagfile' attribute are used,"
		" then all Poses in the input silent file(s) are used." )
		+ Attr( "tagfile", xs_string, "File name whose contents lists a set of whitespace-separated tags"
		" specifying the subset of Poses that should be processed from the input silent file(s). If neither this"
		" attribute, nor the 'tags' attribute are used, then all Poses in the input silent file(s) are used." )
		+ Attr::attribute_w_default( "skip_failed_simulations", xsct_rosetta_bool, "Skip processing of input"
		" Poses if the tag starts with 'W_'", "false" );
	core::io::silent::SilentFileOptions::append_attributes_for_tag_parsing( xsd, attributes );

	pose_inputter_xsd_type_definition_w_attributes(
		xsd,
		keyname(),
		"Inputter for poses originating from silent files. By default, each Pose in the file will"
		" be input and will be used by the JobQueen; however, if you use the 'tags' or 'tagfile' attribute, you"
		" can specify a subset of Poses that will be used instead of the full set.",
		attributes );

}

void
SilentFilePoseInputter::list_options_read( utility::options::OptionKeyList & read_options )
{
	using namespace basic::options::OptionKeys;

	core::io::silent::SilentFileOptions::list_read_options( read_options );
	read_options
		+ in::file::silent
		+ in::file::tags
		+ in::file::tagfile
		+ in::file::skip_failed_simulations;

}

void
SilentFilePoseInputter::initialize_sfd_from_options_and_tag(
	utility::options::OptionCollection const & options,
	utility::tag::TagCOP tag
)
{
	using namespace core::io::silent;
	using namespace basic::options::OptionKeys;

	sf_opts_ = SilentFileOptionsOP( new SilentFileOptions( options ) );
	if ( tag ) sf_opts_->read_from_tag( tag );

	sfd_ = SilentFileDataOP( new SilentFileData( *sf_opts_ ));

	utility::vector1< utility::file::FileName > silent_files;
	utility::vector1< std::string > tags_vector;
	if ( tag ) {
		std::string files = tag->getOption< std::string >( "silent_files", "" );
		if ( files == "" ) {
			throw utility::excn::EXCN_Msg_Exception( "The 'silent_files' attribute must be provided to the SilentFilePoseInputer" );
		}
		utility::vector1< std::string > files_vector = utility::string_split( files, ',', std::string() );
		silent_files.reserve( files_vector.size() );
		for ( auto const & filename : files_vector ) {
			silent_files.push_back( utility::file::FileName( filename ) );
		}
	} else {
		silent_files = options[ in::file::silent ]();
	}

	// Tag subset: the utility::tag::TagCOP takes precedence over anything in the OptionCollection
	if ( tag && tag->hasOption( "tags" ) ) {
		std::string tags = tag->getOption< std::string >( "tags", "" );
		tags_vector = utility::string_split( tags, ',', std::string() );
	} else {
		if ( options[ in::file::tags ].user() ) {
			tags_vector = options[ in::file::tags ];
		}
	}

	std::string tagfile_name;
	if ( tag && tag->hasOption( "tagfile" ) ) {
		tagfile_name = tag->getOption< std::string >( "tagfile" );
	} else if ( options[ in::file::tagfile ].user() ) {
		tagfile_name = options[ in::file::tagfile ]();
	}
	if ( tagfile_name != "" ) {
		utility::io::izstream tag_file( tagfile_name );
		std::copy( std::istream_iterator< std::string >( tag_file ),
			std::istream_iterator< std::string >(),
			std::back_inserter( tags_vector ) );
	}

	initialize_sfd_from_files_and_tags( silent_files, tags_vector );

}

void
SilentFilePoseInputter::initialize_sfd_from_files_and_tags(
	utility::vector1< utility::file::FileName > const & silent_files,
	utility::vector1< std::string > const & tags
)
{
	for ( auto const & silent_file : silent_files ) {
		tr.Debug << "reading " << silent_file << std::endl;
		if ( tags.size() > 0 ) {
			sfd_->read_file( silent_file, tags );
		} else {
			sfd_->read_file( silent_file );
		}
	}
}

protocols::jd3::pose_inputters::PoseInputterOP
SilentFilePoseInputterCreator::create_inputter() const {
	return PoseInputterOP( new SilentFilePoseInputter );
}

std::string
SilentFilePoseInputterCreator::keyname() const
{
	return SilentFilePoseInputter::keyname();
}

void SilentFilePoseInputterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SilentFilePoseInputter::provide_xml_schema(xsd);
}

void SilentFilePoseInputterCreator::list_options_read(
	utility::options::OptionKeyList & read_options
) const
{
	SilentFilePoseInputter::list_options_read( read_options );
}


} // pose_inputters
} // jd3
} // protocols
