// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/pose_inputters/PDBPoseInputter.cc
/// @brief  %PDBPoseInputter class definition
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Steven Lewis (smlewi@gmail.com)

// Unit headers
#include <protocols/jd3/pose_inputters/PDBPoseInputter.hh>
#include <protocols/jd3/pose_inputters/PDBPoseInputterCreator.hh>

// Package headers
#include <protocols/jd3/PoseInputSource.hh>
#include <protocols/jd3/pose_inputters/PoseInputterFactory.hh>
#include <protocols/jd3/pose_inputters/pose_inputter_schemas.hh>

// Project headers
#include <core/pose/Pose.hh>
//#include <core/io/StructFileReaderOptions.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/import_pose_options.hh>
#include <protocols/jd2/util.hh>

//utility headers
#include <utility/assert.hh>
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

namespace protocols {
namespace jd3 {
namespace pose_inputters {


PDBPoseInputter::PDBPoseInputter() {}
PDBPoseInputter::~PDBPoseInputter() {}

bool PDBPoseInputter::job_available_on_command_line() const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	return option[ in::file::s ].user() || option[ in::file::l ].user();
}

PoseInputSources
PDBPoseInputter::pose_input_sources_from_command_line()
{
	utility::vector1< utility::file::FileName > filenames_from_command_line = jd2::input_pdb_files_from_command_line();
	PoseInputSources input_sources;
	input_sources.reserve( filenames_from_command_line.size() );
	for ( core::Size ii = 1; ii <= filenames_from_command_line.size(); ++ii ) {
		utility::file::FileName ii_filename = filenames_from_command_line[ ii ];
		PoseInputSourceOP ii_source( new PoseInputSource( keyname() ));
		ii_source->input_tag( ii_filename.base() );
		ii_source->store_string_pair( "filename", ii_filename.name() );
		input_sources.push_back( ii_source );
	}
	return input_sources;
}

PoseInputSources
PDBPoseInputter::pose_input_sources_from_tag(
	utility::options::OptionCollection const & /*opts*/,
	utility::tag::TagCOP tag
)
{
	// Note: we claim that "path" is a useful tag, so, let's read it at some point!
	PoseInputSources input_sources;
	if ( tag->hasOption( "filename" ) ) {

		utility::file::FileName fname( tag->getOption< std::string >( "filename" ) );
		if ( tag->hasOption( "path" ) ) {
			fname.path( tag->getOption< std::string >( "path" ) );
		}
		PoseInputSourceOP input_source( new PoseInputSource( keyname() ));
		input_source->input_tag( fname.base() );
		input_source->store_string_pair( "filename", fname.name() );
		input_sources.push_back( input_source );

	} else if ( tag->hasOption( "listfile" ) ) {

		std::string list_fname( tag->getOption< std::string >( "listfile" ) );
		utility::io::izstream list_stream( list_fname.c_str() );
		if ( ! list_stream.good() ) {
			throw utility::excn::EXCN_Msg_Exception( "Unable to open list file \"" + list_fname + "\"" );
		}
		std::string line;
		while ( getline( list_stream, line ) ) {
			utility::file::FileName fname( line );
			if ( tag->hasOption( "path" ) ) {
				fname.path( tag->getOption< std::string >( "path" ) );
			}
			PoseInputSourceOP input_source( new PoseInputSource( keyname() ));
			input_source->input_tag( fname.base() );
			input_source->store_string_pair( "filename", fname.name() );
			input_sources.push_back( input_source );
		}

	} else {
		throw utility::excn::EXCN_Msg_Exception( "Did not find either a \"filename\" or a \"listfile\" option in the PDB input tag" );
	}
	return input_sources;
}

/// @details This is only a stub implementation -- the logic here ought to be robust enough
/// to handle non-fullatom poses, but it currently isn't.  This will improve.
/// One thought is that the options that control how poses are loaded from disk could
/// be included in the per-job option set, and then passed into this function via an
/// OptionCollection...
core::pose::PoseOP
PDBPoseInputter::pose_from_input_source(
	PoseInputSource const & input_source,
	utility::options::OptionCollection const & options,
	utility::tag::TagCOP /*tag*/ // possibly null-pointing tag pointer
)
{
	assert( input_source.string_string_map().find( "filename" ) != input_source.string_string_map().end() );
	core::import_pose::ImportPoseOptions import_opts( options );
	return core::import_pose::pose_from_file(
		input_source.string_string_map().find( "filename" )->second,
		import_opts,
		import_opts.read_fold_tree(),
		core::import_pose::PDB_file );
}

std::string
PDBPoseInputter::keyname() { return "PDB"; }

void
PDBPoseInputter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;

	AttributeList input_pdb_attributes;
	input_pdb_attributes
		+ XMLSchemaAttribute( "filename", xs_string , "XRW TO DO" )
		+ XMLSchemaAttribute( "listfile", xs_string , "XRW TO DO" )
		+ XMLSchemaAttribute( "path", xs_string , "XRW TO DO" );

	pose_inputter_xsd_type_definition_w_attributes(
		xsd,
		keyname(),
		"XRW TO DO",
		input_pdb_attributes );

}

void
PDBPoseInputter::list_options_read(
	utility::options::OptionKeyList & read_options
)
{
	core::import_pose::ImportPoseOptions::list_options_read( read_options );
}


PoseInputterOP PDBPoseInputterCreator::create_inputter() const
{
	return PoseInputterOP( new PDBPoseInputter );
}

std::string PDBPoseInputterCreator::keyname() const
{
	return PDBPoseInputter::keyname();
}

void PDBPoseInputterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PDBPoseInputter::provide_xml_schema( xsd );
}

void PDBPoseInputterCreator::list_options_read( utility::options::OptionKeyList & read_options ) const
{
	PDBPoseInputter::list_options_read( read_options );
}

} // namespace pose_inputters
} // namespace jd3
} // namespace protocols

