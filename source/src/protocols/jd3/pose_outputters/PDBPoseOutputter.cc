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
#include <protocols/jd3/pose_outputters/PDBPoseOutputSpecification.hh>
#include <protocols/jd3/JobOutputIndex.hh>
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/InnerLarvalJob.hh>
#include <protocols/jd3/pose_outputters/PoseOutputterFactory.hh>
#include <protocols/jd3/pose_outputters/pose_outputter_schemas.hh>
#include <protocols/jd3/standard/MoverAndPoseJob.hh>

//project headers
#include <core/pose/Pose.hh>
#include <core/io/StructFileRepOptions.hh>
#include <core/io/pdb/pdb_writer.hh>

// ObjexxFCL
#include <ObjexxFCL/string.functions.hh>

// basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/options/keys/OptionKeyList.hh>

static basic::Tracer TR( "protocols.jd3.PDBPoseOutputter" );

namespace protocols {
namespace jd3 {
namespace pose_outputters {

PDBPoseOutputter::PDBPoseOutputter() = default;
PDBPoseOutputter::~PDBPoseOutputter() = default;

bool
PDBPoseOutputter::outputter_specified_by_command_line()
{
	return basic::options::option[ basic::options::OptionKeys::out::pdb ].user() ||
		basic::options::option[ basic::options::OptionKeys::out::pdb_gz ].user();
}

void
PDBPoseOutputter::determine_job_tag(
	utility::tag::TagCOP output_tag,
	utility::options::OptionCollection const & job_options,
	InnerLarvalJob & job
) const
{
	using namespace basic::options;
	utility::tag::TagCOP pdb_tag;
	if ( output_tag && output_tag->hasTag( keyname() ) ) {
		pdb_tag = output_tag->getTag( keyname() );
	}

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

	//Deal with prefix and suffix options.

	if ( pdb_tag and pdb_tag->hasOption( "prefix" ) ) {
		job.job_tag( pdb_tag->getOption< std::string >( "prefix" ) + job.job_tag() );
	} else if ( job_options[ OptionKeys::out::prefix].user() ) {
		job.job_tag( job_options[ OptionKeys::out::prefix]() + job.job_tag() );
	}

	if ( pdb_tag and pdb_tag->hasOption( "suffix ") ) {
		job.job_tag(  job.job_tag() + pdb_tag->getOption< std::string >( "suffix" ) );
	} else if ( job_options[ OptionKeys::out::suffix].user() ) {
		job.job_tag( job.job_tag() + job_options[ OptionKeys::out::suffix]());
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

std::string
PDBPoseOutputter::outputter_for_job(
	PoseOutputSpecification const &
) const
{
	return "";
}

bool PDBPoseOutputter::job_has_already_completed( LarvalJob const & job, utility::options::OptionCollection const & options ) const
{
	using namespace basic::options::OptionKeys;
	InnerLarvalJobCOP inner_job( job.inner_job() );

	utility::tag::TagCOP jobdef_tag( inner_job->jobdef_tag() );

	//TR << "jobdef_tag " << jobdef_tag << std::endl;
	utility::tag::TagCOP pdb_output_tag;
	if ( jobdef_tag && jobdef_tag->hasTag( "Output" ) && jobdef_tag->getTag( "Output" )->hasTag( keyname() ) ) {
		pdb_output_tag = jobdef_tag->getTag( "Output" )->getTag( keyname() );
	}

	// PDB tag "overwrite" attribute takes precedence over the options system
	// if the tag says "do not overwrite," then do not look at the options system.
	if ( pdb_output_tag && pdb_output_tag->hasOption( "overwrite" ) ) {
		if ( pdb_output_tag->getOption< bool >( "overwrite", false ) ) {
			return false;
		}
	} else if ( options[ basic::options::OptionKeys::out::overwrite ] ) { return false; }

	// Assume that each job is only going to produce a single PDB file -- if a job were
	// to form more than one PDB file output, then the output structures that it would write
	// would have different file names.
	std::string filename( output_pdb_name( job, options, pdb_output_tag ));

	bool exists = utility::file::file_exists( filename );
	if ( exists ) {
		TR << std::endl << "Skipping "<< filename << ". Please pass the overwrite option/tag" << std::endl;
	}
	return exists;
}


void PDBPoseOutputter::mark_job_as_having_started( LarvalJob const & /*job*/, utility::options::OptionCollection const & ) const
{
	// STUBBED OUT!
}

std::string
PDBPoseOutputter::class_key() const
{
	return keyname();
}

/// @brief Create the PoseOutputSpecification for a particular job
PoseOutputSpecificationOP
PDBPoseOutputter::create_output_specification(
	LarvalJob const & job,
	JobOutputIndex const & output_index,
	utility::options::OptionCollection const & job_options,
	utility::tag::TagCOP pdb_output_tag // possibly null-pointing tag pointer
)
{
	PDBPoseOutputSpecificationOP spec( new PDBPoseOutputSpecification );
	spec->out_fname( output_pdb_name( job, output_index, job_options, pdb_output_tag ) );
	spec->sfr_opts( core::io::StructFileRepOptions( job_options ) );
	return spec;
}

/// @brief Write a pose out to permanent storage (whatever that may be).
void PDBPoseOutputter::write_output(
	output::OutputSpecification const & spec,
	JobResult const & result
)
{
	using standard::PoseJobResult;
	debug_assert( dynamic_cast< PoseJobResult const * > ( &result ));
	auto const & pose_result( static_cast< PoseJobResult const & > ( result ));
	core::pose::Pose const & pose( *pose_result.pose() );

	using PDBPOS = PDBPoseOutputSpecification;
	debug_assert( dynamic_cast< PDBPOS const * > ( & spec ) );
	PDBPOS pdb_spec( static_cast< PDBPOS const & > ( spec ) );
	core::io::pdb::dump_pdb( pose, pdb_spec.out_fname(), pdb_spec.sfr_opts() );
}

std::string
PDBPoseOutputter::output_pdb_name(
	LarvalJob const & job,
	utility::options::OptionCollection const & options,
	utility::tag::TagCOP tag
) const
{
	JobOutputIndex faux_output_index;
	faux_output_index.primary_output_index = job.nstruct_index();
	faux_output_index.n_primary_outputs = job.nstruct_max();
	return output_pdb_name( job, faux_output_index, options, tag );
}

std::string
PDBPoseOutputter::output_pdb_name(
	LarvalJob const & job,
	JobOutputIndex const & output_index,
	utility::options::OptionCollection const & options,
	utility::tag::TagCOP tag
) const
{
	debug_assert( !tag || tag->getName() ==  keyname() );
	utility::file::FileName fn;
	if ( tag && tag->getOption< bool >( "pdb_gz", false ) ) {
		fn.ext( ".pdb.gz" );
	} else if ( options[ basic::options::OptionKeys::out::pdb_gz ] ) {
		fn.ext( ".pdb.gz" );
	} else {
		fn.ext( ".pdb" );
	}

	fn.base( ( job.status_prefix() == "" ? "" : ( job.status_prefix() + "_" ) )
		+ job.job_tag_with_index_suffix( output_index ) );

	// Priority: ask for the options in order:
	// 1. The path specified in the tag,
	// 2. The path specified in out::path::pdb <-- most specific
	// 3. The path specified in out::path::all <-- more general
	// 4. The path specified in out::path::path <-- dunno if this is more or less general than out::path::all?

	if ( tag && tag->hasOption( "path" ) ) {
		fn.path( tag->getOption< std::string >( "path" ) );
	} else if ( options[ basic::options::OptionKeys::out::path::pdb ].user() ) {
		fn.path( options[ basic::options::OptionKeys::out::path::pdb ]() );
	} else if ( options[ basic::options::OptionKeys::out::path::all ].user() ) {
		fn.path( options[ basic::options::OptionKeys::out::path::all ]() );
	} else if ( options[ basic::options::OptionKeys::out::path::path ].user() ) {
		fn.path( options[ basic::options::OptionKeys::out::path::path ]() );
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
		" combined with the 'filename' attribute. Cannot be combined with either the 'prefix' or 'suffix' attributes."
		"  Cannot be combined with either the out:prefix or out:suffix flags that might be provided on the command line" )
		+ XMLSchemaAttribute( "path", xs_string , "Give the directory to which the output .pdb file should be written."
		" Note that the output path does not become part of the job name, so if you have two jobs with the same job"
		" name written to different directories, then your log file and your score file (and any other secondary pose"
		" outputter) will not distinguish between which of the two jobs it is writing output for" )
		+ XMLSchemaAttribute( "overwrite", xsct_rosetta_bool , "If this is set to 'true', then the job(s) will run"
		" even if an output file with the name that this job would produce exists, and that previously-existing"
		" output file will be overwritten with the new output file." )
		+ XMLSchemaAttribute::attribute_w_default( "pdb_gz", xsct_rosetta_bool, "Should the output PDB file be written as a .gz?", "false" )
		+ XMLSchemaAttribute( "prefix", xs_string, "Set output PDB Prefix. Cannot be combined with the 'filename_pattern' attribute")
		+ XMLSchemaAttribute( "suffix", xs_string, "Set output PDB Suffix. Cannot be combined with the 'filename_pattern' attribute");

	pose_outputter_xsd_type_definition_w_attributes( xsd, keyname(),
		"The (typically) default PoseOutputter that writes the structure out as a PDB-formatted file", output_pdb_attributes );
}

void
PDBPoseOutputter::list_options_read(
	utility::options::OptionKeyList & read_options
)
{
	core::io::StructFileRepOptions::list_options_read( read_options );
	read_options
		+ basic::options::OptionKeys::out::pdb_gz
		+ basic::options::OptionKeys::out::overwrite
		+ basic::options::OptionKeys::out::prefix
		+ basic::options::OptionKeys::out::suffix
		+ basic::options::OptionKeys::out::path::path
		+ basic::options::OptionKeys::out::path::pdb
		+ basic::options::OptionKeys::out::path::all;
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
