// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/SilentFilePoseOutputter.cc
/// @brief  Definition of the %SilentFilePoseOutputter class's methods
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), Andy Watkins (amw579@stanford.edu)

//unit headers
#include <protocols/jd3/pose_outputters/SilentFilePoseOutputter.hh>
#include <protocols/jd3/pose_outputters/SilentFilePoseOutputterCreator.hh>

//package headers
#include <protocols/jd3/pose_outputters/SilentFilePoseOutputSpecification.hh>
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/InnerLarvalJob.hh>
#include <protocols/jd3/pose_outputters/PoseOutputterFactory.hh>
#include <protocols/jd3/pose_outputters/pose_outputter_schemas.hh>
#include <protocols/jd3/full_model_inputters/FullModelInputSource.hh>
#include <protocols/jd3/full_model/MoverAndFullModelJob.hh>
#include <protocols/jd3/standard/MoverAndPoseJob.hh>

//project headers
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/io/StructFileRepOptions.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentStructFactory.hh>

// ObjexxFCL
#include <ObjexxFCL/string.functions.hh>

// basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
//#include <basic/options/keys/run.OptionKeys.gen.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/keys/OptionKeyList.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/FileName.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>


namespace protocols {
namespace jd3 {
namespace pose_outputters {

using namespace core::io;
using namespace core::io::silent;

SilentFilePoseOutputter::SilentFilePoseOutputter() :
	buffer_limit_( 20 )
{}

SilentFilePoseOutputter::~SilentFilePoseOutputter() = default;

bool
SilentFilePoseOutputter::outputter_specified_by_command_line()
{
	return basic::options::option[ basic::options::OptionKeys::out::file::silent ].user();
}

void
SilentFilePoseOutputter::determine_job_tag(
	utility::tag::TagCOP /*output_tag*/,
	utility::options::OptionCollection const & /*job_options*/,
	InnerLarvalJob & job
) const {
	job.job_tag( job.input_tag() );
}

/// @details If this is just being given an InnerLarvalJob, that means that the
/// outputter is at most a function of the input source--as there's a bijection
/// b/n ILJ and input source. So, let's just use the job_tag. What else can we do?
std::string
SilentFilePoseOutputter::outputter_for_job(
	utility::tag::TagCOP sf_tag,
	utility::options::OptionCollection const & opts,
	InnerLarvalJob const &
) const
{
	if ( sf_tag ) {
		runtime_assert( sf_tag->hasOption( "filename" ));
		return sf_tag->getOption< std::string >( "filename" );
	} else {
		return opts[ basic::options::OptionKeys::out::file::silent ];
	}
}


std::string
SilentFilePoseOutputter::outputter_for_job(
	PoseOutputSpecification const & spec
) const
{
	using SFPOS = SilentFilePoseOutputSpecification;
	debug_assert( dynamic_cast< SFPOS const * > ( &spec ) );
	auto const & sf_spec( static_cast< SFPOS const & > ( spec ) );

	return sf_spec.out_fname();
}

bool SilentFilePoseOutputter::job_has_already_completed( LarvalJob const & /*job*/, utility::options::OptionCollection const & ) const
{
	// There is no "don't overwrite already generated structures" behavior for silent-file output.
	return false;
}


void SilentFilePoseOutputter::mark_job_as_having_started( LarvalJob const & /*job*/, utility::options::OptionCollection const & ) const
{
	// This is not a behavior supported by the SilentFilePoseOutputter
}

std::string
SilentFilePoseOutputter::class_key() const
{
	return keyname();
}


/// @brief Create the PoseOutputSpecification for a particular job
PoseOutputSpecificationOP
SilentFilePoseOutputter::create_output_specification(
	LarvalJob const & job,
	JobOutputIndex const & output_index,
	utility::options::OptionCollection const & job_options,
	utility::tag::TagCOP tag // possibly null-pointing tag pointer
)
{
	using namespace core::io::silent;
	debug_assert( !tag || tag->getName() == keyname() ); // I expect this Tag to point to my data

	core::io::silent::SilentFileOptions opts( job_options );
	std::string fname_out;
	core::Size buffer_limit(0);
	if ( tag ) {
		opts.read_from_tag( tag );
		fname_out = tag->getOption< std::string >( "filename" );
		if ( tag->hasOption( "buffer_limit" ) ) {
			buffer_limit = tag->getOption< core::Size >( "buffer_limit" );
		}
	} else {
		using namespace basic::options::OptionKeys;
		fname_out = job_options[ out::file::silent ]();
		// buffer limit? -- default of 0?
	}

	SilentFilePoseOutputSpecificationOP sf_pos( new SilentFilePoseOutputSpecification );
	sf_pos->sf_opts( opts );
	sf_pos->out_fname( fname_out );
	sf_pos->buffer_limit( buffer_limit );
	// There should ultimately be a better way of communicating. Meantime,
	// look for an option that's ubiquitous in stepwise runs but used little if at all
	// elsewhere.
	// Try to cast the input source?
	if ( dynamic_cast< full_model_inputters::FullModelInputSource const * >( &job.inner_job()->input_source() ) ) { //job_options[ ] ) {
		sf_pos->pose_tag( ( job.status_prefix() == "" ? "" : ( job.status_prefix() + "_" ) )
			+ job.job_tag_with_index_suffix( output_index, 6 ) );
	} else {
		sf_pos->pose_tag( ( job.status_prefix() == "" ? "" : ( job.status_prefix() + "_" ) )
			+ job.job_tag_with_index_suffix( output_index ) );
	}

	return sf_pos;
}

/// @brief Write a pose out to permanent storage (whatever that may be).
void
SilentFilePoseOutputter::write_output(
	output::OutputSpecification const & spec,
	JobResult const & result
)
{
	using namespace protocols::jd3::full_model;

	using standard::PoseJobResult;
	debug_assert( dynamic_cast< PoseJobResult const * > ( &result ) || dynamic_cast< FullModelJobResult const * > ( &result ) );
	auto const & pose_result( static_cast< PoseJobResult const & > ( result ));
	core::pose::Pose const & pose( *pose_result.pose() );

	using SFPOS = SilentFilePoseOutputSpecification;
	debug_assert( dynamic_cast< SFPOS const * > ( &spec ) );
	auto const & sf_spec( static_cast< SFPOS const & > ( spec ) );

	if ( ! opts_ ) {
		// I.e., we are a fresh outputter
		opts_.reset( new core::io::silent::SilentFileOptions( sf_spec.sf_opts() ) );
		fname_out_ = sf_spec.out_fname();
		buffer_limit_ = sf_spec.buffer_limit();
	}

	core::io::silent::SilentStructOP ss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out( pose, sf_spec.sf_opts() );
	ss->fill_struct( pose, sf_spec.pose_tag() );

	if ( dynamic_cast< FullModelJobResult const * >( &result ) && core::pose::full_model_info::full_model_info_defined( pose ) ) {
		ss->add_string_value( "missing", ObjexxFCL::string_of( core::pose::full_model_info::get_number_missing_residues_and_connections( pose ) ), -3 );
	}

	buffered_structs_.push_back( ss );

	if ( buffered_structs_.size() >= buffer_limit_ ) {
		flush();
	}

}

//void
//SilentFilePoseOutputter::write_output_pose(
// LarvalJob const & job,
// JobOutputIndex const & output_index,
// utility::options::OptionCollection const & job_options,
// utility::tag::TagCOP tag, // possibly null-pointing tag pointer
// core::pose::Pose const & pose
//)
//{
// if ( ! opts_ ) {
//  initialize_sf_options( job_options, tag );
// }
//
// std::string output_tag = ( job.status_prefix() == "" ? "" : ( job.status_prefix() + "_" ) )
//  + job.job_tag_with_index_suffix( output_index );
//
// ss->fill_struct( pose, output_tag );
//
//}

void SilentFilePoseOutputter::flush()
{
	if ( ! opts_ ) {
		debug_assert( buffered_structs_.size() == 0 );
		return;
	}

	// Here or in SilentFileData, we also need to EXPLICITLY output OTHER_STRUCTs
	core::io::silent::SilentFileData sfd( *opts_ );
	for ( auto const & iter : buffered_structs_ ) {
		sfd.add_structure( *iter );
		// One must explicitly add these because of the structure of
		// SilentFileData. Maybe it can or should change! Because then on
		// the other side, it has to change whether the struct being
		// added is an other_struct. This perturbs the code less though.
		for ( auto const & other_struct : iter->other_struct_list() ) {
			sfd.add_structure( *other_struct );
		}
	}
	sfd.write_all( fname_out_ );
	buffered_structs_.clear();
}

std::string
SilentFilePoseOutputter::keyname() { return "SilentFile"; }

void
SilentFilePoseOutputter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;

	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::required_attribute( "filename", xs_string , "The name of the output silent file that should be written to." )
		+ XMLSchemaAttribute( "buffer_limit", xsct_non_negative_integer, "The number of Poses that should be held in memory between each write to disk" );
	core::io::silent::SilentFileOptions::append_attributes_for_tag_parsing( xsd, attributes );


	pose_outputter_xsd_type_definition_w_attributes( xsd, keyname(),
		"A PoseOutputter that writes structures out in a rosetta-specific format; a single "
		" silent file can hold hundreds or thousdands of output structures, lessening the load"
		" on file systems, and making output management easier. Note that if two different Jobs"
		" defined in the Job-definition file intend to write their outputs to the same file, then"
		" the settings for the first Job will take precedence over the settings for the second Job."
		" This situation is complicated further if using MPI and multiple output/archive nodes:"
		" in this scenario, multiple silent files will be written (one per output node) and the first"
		" Job to be written to a silent file will determine which settings take precedence but that"
		" it is not knowable which of the two Jobs will be written first. To avoid confusion, it is"
		" recommended that either each Job write their outputs to a different file, or that the"
		" same options are used for all Jobs writing to the same file.",
		attributes );
}

void
SilentFilePoseOutputter::list_options_read(
	utility::options::OptionKeyList & read_options
)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::io::silent::SilentFileOptions::list_read_options( read_options );
	read_options
		+ out::silent_gz
		+ out::file::silent;
}


PoseOutputterOP SilentFilePoseOutputterCreator::create_outputter() const
{
	return PoseOutputterOP( new SilentFilePoseOutputter );
}

std::string SilentFilePoseOutputterCreator::keyname() const
{
	return SilentFilePoseOutputter::keyname();
}

void SilentFilePoseOutputterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SilentFilePoseOutputter::provide_xml_schema( xsd );
}

void SilentFilePoseOutputterCreator::list_options_read( utility::options::OptionKeyList & read_options ) const
{
	SilentFilePoseOutputter::list_options_read( read_options );
}

bool
SilentFilePoseOutputterCreator::outputter_specified_by_command_line() const
{
	return SilentFilePoseOutputter::outputter_specified_by_command_line();
}


} // namespace pose_outputters
} // namespace jd3
} // namespace protocols
