// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/mmTFPoseOutputter.cc
/// @brief  Definition of the %mmTFPoseOutputter class's methods
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com) - PDBPoseOutputter bases of this class

//unit headers
#include <protocols/jd3/pose_outputters/mmTFPoseOutputter.hh>
#include <protocols/jd3/pose_outputters/mmTFPoseOutputterCreator.hh>

//package headers
#include <protocols/jd3/pose_outputters/mmTFPoseOutputSpecification.hh>
#include <protocols/jd3/JobOutputIndex.hh>
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/InnerLarvalJob.hh>
#include <protocols/jd3/pose_outputters/PoseOutputterFactory.hh>
#include <protocols/jd3/pose_outputters/pose_outputter_schemas.hh>
#include <protocols/jd3/job_results/PoseJobResult.hh>

//project headers
#include <core/pose/Pose.hh>
#include <core/io/StructFileRepOptions.hh>
#include <core/io/mmtf/mmtf_writer.hh>

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

static basic::Tracer TR( "protocols.jd3.mmTFPoseOutputter" );

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {
namespace pose_outputters {

mmTFPoseOutputter::mmTFPoseOutputter():
	protocols::jd3::pose_outputters::PDBPoseOutputter()
{

}

mmTFPoseOutputter::~mmTFPoseOutputter() = default;

std::string
mmTFPoseOutputter::class_key() const
{
	return keyname();
}

bool
mmTFPoseOutputter::outputter_specified_by_command_line()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	return ( option[ out::mmtf ].user() && option[ out::mmtf ].value() );
}

/// @brief Create the PoseOutputSpecification for a particular job
PoseOutputSpecificationOP
mmTFPoseOutputter::create_output_specification(
	LarvalJob const & job,
	JobOutputIndex const & output_index,
	utility::options::OptionCollection const & job_options,
	utility::tag::TagCOP pdb_output_tag // possibly null-pointing tag pointer
)
{
	mmTFPoseOutputSpecificationOP spec( new mmTFPoseOutputSpecification );
	spec->out_fname( output_name( job, output_index, job_options, pdb_output_tag ) );
	spec->sfr_opts( core::io::StructFileRepOptions( job_options ) );
	return spec;
}

/// @brief Write a pose out to permanent storage (whatever that may be).
void mmTFPoseOutputter::write_output(
	output::OutputSpecification const & spec,
	JobResult const & result
)
{
	using job_results::PoseJobResult;
	debug_assert( dynamic_cast< PoseJobResult const * > ( &result ));
	auto const & pose_result( static_cast< PoseJobResult const & > ( result ));
	core::pose::Pose const & pose( *pose_result.pose() );

	using mmtfPOS = mmTFPoseOutputSpecification;
	debug_assert( dynamic_cast< mmtfPOS const * > ( & spec ) );
	mmtfPOS mmtf_spec( static_cast< mmtfPOS const & > ( spec ) );

	//Create directory if it does not exist
	utility::file::FileName full_fname( mmtf_spec.out_fname());
	std::string fpath = full_fname.path();
	if ( ! fpath.empty() && ! utility::file::is_directory(fpath) ) {
		bool success = utility::file::create_directory(fpath);
		if ( !success ) {
			utility_exit_with_message("Directory "+fpath+"for mmTF "+mmtf_spec.out_fname()+" is not present and could not be created" );
		}
	}

	core::io::mmtf::dump_mmtf( pose, mmtf_spec.out_fname(), mmtf_spec.sfr_opts() );
}

std::string
mmTFPoseOutputter::output_name(
	LarvalJob const & job,
	JobOutputIndex const & output_index,
	utility::options::OptionCollection const & options,
	utility::tag::TagCOP tag
) const
{
	debug_assert( !tag || tag->getName() ==  keyname() );
	utility::file::FileName fn;
	if ( tag && tag->getOption< bool >( "mmtf_gz", false ) ) {
		fn.ext( ".mmtf.gz" );
	} else if ( options[ basic::options::OptionKeys::out::mmtf_gz ]() ) {
		fn.ext( ".mmtf.gz" );
	} else {
		fn.ext( ".mmtf" );
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
	} else if ( options[ basic::options::OptionKeys::out::path::all ].user() ) {
		fn.path( options[ basic::options::OptionKeys::out::path::all ]() );
	} else if ( options[ basic::options::OptionKeys::out::path::path ].user() ) {
		fn.path( options[ basic::options::OptionKeys::out::path::path ]() );
	}

	return fn();
}

std::string
mmTFPoseOutputter::keyname() { return "mmTF"; }

void
mmTFPoseOutputter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
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
		" '1abc_steal_native_frags_0001.mmtf', '1abc_steal_native_frags_0002.pdb', ...,"
		" '2def_steal_native_frags_0001.mmtf', etc. if it were used with input tags '1abc' and '2def'."
		" prefix and suffix options either on the cmd line or here will be combined appropriately with this option - "
		" added to the beginning or end of the final name respectively.  If our prefix was `prefix_` and our suffix was "
		" `suffix_`, then our final PDB name would be: prefix_2def_steal_native_frags_suffix_0001`  ")

		+ XMLSchemaAttribute( "path", xs_string , "Give the directory to which the output .mmtf file should be written."
		" Note that the output path does not become part of the job name, so if you have two jobs with the same job"
		" name written to different directories, then your log file and your score file (and any other secondary pose"
		" outputter) will not distinguish between which of the two jobs it is writing output for" )

		+ XMLSchemaAttribute( "overwrite", xsct_rosetta_bool , "If this is set to 'true', then the job(s) will run"
		" even if an output file with the name that this job would produce exists, and that previously-existing"
		" output file will be overwritten with the new output file." )

		+ XMLSchemaAttribute::attribute_w_default( "mmtf_gz", xsct_rosetta_bool, "Should the output mmTF file be written as a .gz?", "false" )

		+ XMLSchemaAttribute( "prefix", xs_string, "Set output mmTF Prefix. Can be combined with the 'filename_pattern' attribute. "
		"Overrides any cmd-line prefix option set.")

		+ XMLSchemaAttribute( "suffix", xs_string, "Set output mmTF Suffix. Can be combined with the 'filename_pattern' attribute. "
		"Overrides any cmd-line prefix option set.");

	pose_outputter_xsd_type_definition_w_attributes( xsd, keyname(),
		"The (typically) default PoseOutputter that writes the structure out as a mmTF-formatted file", output_pdb_attributes );
}

void
mmTFPoseOutputter::list_options_read(
	utility::options::OptionKeyList & read_options
)
{
	core::io::StructFileRepOptions::list_options_read( read_options );
	read_options
		+ basic::options::OptionKeys::out::mmtf
		+ basic::options::OptionKeys::out::mmtf_gz
		+ basic::options::OptionKeys::out::overwrite
		+ basic::options::OptionKeys::out::prefix
		+ basic::options::OptionKeys::out::suffix
		+ basic::options::OptionKeys::out::path::path
		+ basic::options::OptionKeys::out::path::all;
}

PoseOutputterOP mmTFPoseOutputterCreator::create_outputter() const
{
	return utility::pointer::make_shared< mmTFPoseOutputter >();
}

std::string mmTFPoseOutputterCreator::keyname() const
{
	return mmTFPoseOutputter::keyname();
}

void mmTFPoseOutputterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	mmTFPoseOutputter::provide_xml_schema( xsd );
}

void mmTFPoseOutputterCreator::list_options_read( utility::options::OptionKeyList & read_options ) const
{
	mmTFPoseOutputter::list_options_read( read_options );
}

bool
mmTFPoseOutputterCreator::outputter_specified_by_command_line() const
{
	return mmTFPoseOutputter::outputter_specified_by_command_line();
}


} // namespace pose_outputters
} // namespace jd3
} // namespace protocols

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::jd3::pose_outputters::mmTFPoseOutputter::save( Archive & arc ) const {
	arc( cereal::base_class< protocols::jd3::pose_outputters::PDBPoseOutputter >( this ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::jd3::pose_outputters::mmTFPoseOutputter::load( Archive & arc ) {
	arc( cereal::base_class< protocols::jd3::pose_outputters::PDBPoseOutputter >( this ) );
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::pose_outputters::mmTFPoseOutputter );
CEREAL_REGISTER_TYPE( protocols::jd3::pose_outputters::mmTFPoseOutputter )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_jd3_pose_outputters_mmTFPoseOutputter )
#endif // SERIALIZATION
