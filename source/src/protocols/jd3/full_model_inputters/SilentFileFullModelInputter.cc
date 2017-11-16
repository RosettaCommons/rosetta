// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/SilentFileFullModelInputter.cc
/// @brief
/// @author Andy Watkins (amw579@stanford.edu)

///Unit headers
#include <protocols/jd3/full_model_inputters/SilentFileFullModelInputter.hh>
#include <protocols/jd3/full_model_inputters/SilentFileFullModelInputterCreator.hh>

// Package headers
#include <protocols/jd3/full_model_inputters/FullModelInputSource.hh>
#include <protocols/jd3/full_model_inputters/FullModelInputterFactory.hh>
#include <protocols/jd3/full_model_inputters/full_model_inputter_schemas.hh>

///Project headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
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

static basic::Tracer tr( "protocols.jd3.SilentFileFullModelInputter" );

namespace protocols {
namespace jd3 {
namespace full_model_inputters {

SilentFileFullModelInputter::SilentFileFullModelInputter()
{
	tr.Debug << "Instantiate SilentFileFullModelInputter" << std::endl;
}

SilentFileFullModelInputter::~SilentFileFullModelInputter() {}

bool SilentFileFullModelInputter::job_available_on_command_line() const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	return option[ in::file::silent ].user();
}

FullModelInputSources SilentFileFullModelInputter::full_model_input_sources_from_command_line()
{
	using namespace core::io::silent;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	initialize_sfd_from_options_and_tag( basic::options::option, utility::tag::TagOP() );

	// From stepwise convention: when you send in a silent file for input, you ONLY READ THE FIRST.
	// This sounds counterintuitive. It's not the way it's done in rna_denovo! But this is necessary
	// for some silent file formats / integration tests... for now.
	// AMW TODO: investigate how this does for inputting a vector of silent files....

	FullModelInputSources input_sources;
	input_sources.emplace_back( new FullModelInputSource( keyname() ) );

	for ( auto iter : *sfd_ ) {
		std::string const & tag = iter->decoy_tag();
		if ( boost::starts_with( tag, "W_" ) && option[ in::file::skip_failed_simulations ] ) {
			continue;
		}

		input_sources[1]->input_tag( tag );
		tr << "Tag: " << tag << std::endl;
		break; // only one!
	}
	return input_sources;
}

FullModelInputSources
SilentFileFullModelInputter::full_model_input_sources_from_tag(
	utility::options::OptionCollection const & opts,
	utility::tag::TagCOP tag
)
{
	using namespace core::io::silent;
	using namespace basic::options::OptionKeys;

	initialize_sfd_from_options_and_tag( opts, tag );

	FullModelInputSources input_sources;
	input_sources.emplace_back( new FullModelInputSource( keyname() ) );
	for ( auto iter : *sfd_ ) {
		std::string const & decoy_tag = iter->decoy_tag();
		if ( boost::starts_with( decoy_tag, "W_" ) && ( opts[ in::file::skip_failed_simulations]
				|| tag->getOption< bool >( "skip_failed_simulations", false ) ) ) {
			continue;
		}

		input_sources[1]->input_tag( decoy_tag );
		break; // only one!
	}
	return input_sources;
}


/// @details This function will first see if the pose already exists in the Job.
/// If not, it will read it into the pose reference, and hand a COP cloned from
/// that pose to the Job. If the pose pre-exists it just copies the COP's pose
/// into it.
core::pose::PoseOP
SilentFileFullModelInputter::full_model_from_input_source(
	FullModelInputSource const & input_source,
	utility::options::OptionCollection const & options,
	utility::tag::TagCOP tag // possibly null-pointing tag pointer
)
{
	if ( ! sfd_ ) { initialize_sfd_from_options_and_tag( options, tag ); }

	debug_assert( sfd_->has_tag( input_source.input_tag() ));

	using namespace core::io::silent;
	//SilentStruct const & silent_struct( sfd_->get_structure( input_source.input_tag() ) );
	core::chemical::ResidueTypeSetCOP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	/*core::pose::PoseOP pose( new core::pose::Pose );
	silent_struct.fill_pose( *pose );
	return pose;
	*/


	// AMW EDIT HERE!!!
	//debug_assert( input_source.string_string_map().find( "filename" ) != input_source.string_string_map().end() );
	//core::import_pose::ImportPoseOptions import_opts( options );

	return core::import_pose::initialize_pose_and_other_poses_from_options( rsd_set, options );

}

std::string SilentFileFullModelInputter::keyname() { return "Silent"; }

/// @brief returns the schema for the PDB element used in a job-definition file
/// including all options that govern how a PDB is loaded.
void
SilentFileFullModelInputter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
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

	full_model_inputter_xsd_type_definition_w_attributes(
		xsd,
		keyname(),
		"Inputter for poses originating from silent files. By default, each Pose in the file will"
		" be input and will be used by the JobQueen; however, if you use the 'tags' or 'tagfile' attribute, you"
		" can specify a subset of Poses that will be used instead of the full set.",
		attributes );

}

void
SilentFileFullModelInputter::list_options_read( utility::options::OptionKeyList & read_options )
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
SilentFileFullModelInputter::initialize_sfd_from_options_and_tag(
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
			throw utility::excn::EXCN_Msg_Exception( "The 'silent_files' attribute must be provided to the SilentFileFullModelInputer" );
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
SilentFileFullModelInputter::initialize_sfd_from_files_and_tags(
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

protocols::jd3::full_model_inputters::FullModelInputterOP
SilentFileFullModelInputterCreator::create_inputter() const {
	return FullModelInputterOP( new SilentFileFullModelInputter );
}

std::string
SilentFileFullModelInputterCreator::keyname() const
{
	return SilentFileFullModelInputter::keyname();
}

void SilentFileFullModelInputterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SilentFileFullModelInputter::provide_xml_schema(xsd);
}

void SilentFileFullModelInputterCreator::list_options_read(
	utility::options::OptionKeyList & read_options
) const
{
	SilentFileFullModelInputter::list_options_read( read_options );
}


} // full_model_inputters
} // jd3
} // protocols
