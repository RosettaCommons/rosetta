// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/PDBPoseOutputter.cc
/// @brief  Definition of the %PDBPoseOutputter class's methods
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

//unit headers
#include <protocols/jd3/pose_outputters/PDBPoseOutputter.hh>
#include <protocols/jd3/pose_outputters/PDBPoseOutputterCreator.hh>

//package headers
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/InnerLarvalJob.hh>
#include <protocols/jd3/pose_outputters/PoseOutputterFactory.hh>
#include <protocols/jd3/pose_outputters/pose_outputter_schemas.hh>

//project headers
#include <core/pose/Pose.hh>
#include <core/io/StructFileRepOptions.hh>
#include <core/io/pdb/pdb_writer.hh>

// ObjexxFCL
#include <ObjexxFCL/string.functions.hh>

// basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/FileName.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/options/keys/OptionKeyList.hh>

namespace protocols {
namespace jd3 {
namespace pose_outputters {

PDBPoseOutputter::PDBPoseOutputter() {}
PDBPoseOutputter::~PDBPoseOutputter() {}

bool
PDBPoseOutputter::outputter_specified_by_command_line()
{
	return basic::options::option[ basic::options::OptionKeys::out::pdb ].user() ||
		basic::options::option[ basic::options::OptionKeys::out::pdb_gz ].user();
}

void
PDBPoseOutputter::determine_job_tag(
	utility::tag::TagCOP pdb_tag,
	utility::options::OptionCollection const & /*job_options*/,
	InnerLarvalJob & job
) const
{
	if ( pdb_tag ) {
		using namespace utility::tag;
		if ( pdb_tag->hasOption( "filename" ) ) {
			utility::file::FileName fname( pdb_tag->getOption< std::string >( "filename" ));
			job.job_tag( fname.base() );
		} else if ( pdb_tag->hasOption( "filename_pattern" ) ) {
			// replace the single "$" in the filename_pattern with the input_tag already set for the inner larval job
			// we know there is only a single dollar sign because the "string_w_one_dollarsign" restriction is
			// formulated in such a way that forbids multiple dollar signs and requires there to be at least one
			std::string pattern = pdb_tag->getOption< std::string >( "filename_pattern" );
			std::string job_tag = utility::replace_in( pattern, "$", job.input_tag() );
			job.job_tag( job_tag );
		}
	} else {
		job.job_tag( job.input_tag() );
	}
}

/// @details In returning the empty string, we signal to the JobQueen (or any other
/// user) that all %PDBPoseOutputters are interchangable -- that they do not buffer data
/// destined for a single file, and so a single %PDBPoseOutputter may be used for all
/// output Poses.
std::string
PDBPoseOutputter::outputter_for_job(
	utility::tag::TagCOP,
	utility::options::OptionCollection const &,
	InnerLarvalJob const &
) const
{
	return "";
}

bool PDBPoseOutputter::job_has_already_completed( LarvalJob const & /*job*/ ) const
{
	// STUBBED OUT!
	return false;
}


void PDBPoseOutputter::mark_job_as_having_started( LarvalJob const & /*job*/ ) const
{
	// STUBBED OUT!
}

std::string
PDBPoseOutputter::class_key() const
{
	return keyname();
}

void
PDBPoseOutputter::write_output_pose(
	LarvalJob const & job,
	utility::options::OptionCollection const & job_options,
	utility::tag::TagCOP tag, // possibly null-pointing tag pointer
	core::pose::Pose const & pose )
{
	std::string out_fname = output_pdb_name( job, job_options, tag );
	core::io::StructFileRepOptionsOP sfr_opts( new core::io::StructFileRepOptions( job_options ) );
	core::io::pdb::dump_pdb( pose, out_fname, sfr_opts );
}

std::string
PDBPoseOutputter::output_pdb_name(
	LarvalJob const & job,
	utility::options::OptionCollection const & options,
	utility::tag::TagCOP tag
) const
{
	utility::file::FileName fn;
	if ( tag && tag->getOption< bool >( "pdb_gz", false )) {
		fn.ext( ".pdb.gz" );
	} else if ( options[ basic::options::OptionKeys::out::pdb_gz ] ) {
		fn.ext( ".pdb.gz" );
	} else {
		fn.ext( ".pdb" );
	}

	fn.base( ( job.status_prefix() == "" ? "" : ( job.status_prefix() + "_" ) )
		+ job.nstruct_suffixed_job_tag() );

	if ( tag && tag->hasOption( "path" ) ) {
		fn.path( tag->getOption< std::string >( "path" ) );
	}
	return fn();
}

void PDBPoseOutputter::flush() {}

std::string
PDBPoseOutputter::keyname() { return "PDB"; }

void
PDBPoseOutputter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;

	XMLSchemaRestriction string_w_one_dollarsign;
	string_w_one_dollarsign.name( "string_w_one_dollarsign" );
	string_w_one_dollarsign.base_type( xs_string );
	// this reads "anything but a dollar sign or nothing followed by a dollar sign followed by anything but a dollar sign or nothing
	string_w_one_dollarsign.add_restriction( xsr_pattern, "[^$]*$[^$]*" );
	xsd.add_top_level_element( string_w_one_dollarsign );

	AttributeList output_pdb_attributes;
	output_pdb_attributes
		+ XMLSchemaAttribute( "filename", xs_string , "The name to of the file write the output structure to"
		" -- only works correctly so long as there is only one input structure, otherwise the output"
		" structures would pile up on top of each other (the input tags distinguishing them are ignored.)."
		" Should not include any directory names; use the 'path' attribute instead for that."
		" Cannot be combined with the 'filename_pattern' attribute, which is typically preferrable to this one." )
		+ XMLSchemaAttribute( "filename_pattern", "string_w_one_dollarsign",
		"If you want to name the output pdb files for a job with some permutation on the input tag"
		" (i.e. the input pdb name) and then something that identifies something particular about"
		" the job (e.g. '1abc_steal_native_frags_0001.pdb') then use the filename_pattern string. It"
		" expects a string that has a single dolar sign; the original job tag will be substituted for"
		" the dolar sign. E.g. '$_steal_native_frags' would produce pdbs named"
		" '1abc_steal_native_frags_0001.pdb', '1abc_steal_native_frags_0002.pdb', ...,"
    " '2def_steal_native_frags_0001.pdb', etc. if it were used with input tags '1abc' and '2def'. Cannot be"
		" combined with the 'filename' attribute." )
		+ XMLSchemaAttribute( "path", xs_string , "XRW TO DO" )
	  + XMLSchemaAttribute::attribute_w_default( "pdb_gz", xsct_rosetta_bool, "Should the output PDB file be written as a .gz?", "false" );

	pose_outputter_xsd_type_definition_w_attributes( xsd, keyname(),
		"The (typically) default PoseOutputter that writes the structure out as a PDB-formatted file", output_pdb_attributes );
}

void
PDBPoseOutputter::list_options_read(
	utility::options::OptionKeyList & read_options
)
{
	core::io::StructFileRepOptions::list_options_read( read_options );
	read_options + basic::options::OptionKeys::out::pdb_gz;
}

PoseOutputterOP PDBPoseOutputterCreator::create_outputter() const
{
	return PoseOutputterOP( new PDBPoseOutputter );
}

std::string PDBPoseOutputterCreator::keyname() const
{
	return PDBPoseOutputter::keyname();
}

void PDBPoseOutputterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PDBPoseOutputter::provide_xml_schema( xsd );
}

void PDBPoseOutputterCreator::list_options_read( utility::options::OptionKeyList & read_options ) const
{
	PDBPoseOutputter::list_options_read( read_options );
}

bool
PDBPoseOutputterCreator::outputter_specified_by_command_line() const
{
	return PDBPoseOutputter::outputter_specified_by_command_line();
}


} // namespace pose_outputters
} // namespace jd3
} // namespace protocols
