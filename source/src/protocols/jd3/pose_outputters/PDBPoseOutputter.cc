// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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


namespace protocols {
namespace jd3 {
namespace pose_outputters {

PDBPoseOutputter::PDBPoseOutputter() {}
PDBPoseOutputter::~PDBPoseOutputter() {}

bool
PDBPoseOutputter::outputter_specified_by_command_line() const
{
	return basic::options::option[ basic::options::OptionKeys::out::pdb ].user() ||
		basic::options::option[ basic::options::OptionKeys::out::pdb_gz ].user();
}

void
PDBPoseOutputter::determine_job_tag(
	utility::tag::TagCOP output_tag,
	utility::options::OptionCollection const & /*job_options*/,
	InnerLarvalJob & job
) const
{
	if ( output_tag ) {
		using namespace utility::tag;
		runtime_assert( output_tag->hasTag( keyname() ) );
		TagCOP pdb_tag = output_tag->getTag( keyname() );
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
	core::pose::Pose const & pose )
{
	std::string out_fname = output_pdb_name( job );
	utility::io::ozstream ostream( out_fname );

	core::io::StructFileRepOptionsOP sfr_opts( new core::io::StructFileRepOptions( job_options ) );
	core::io::pdb::dump_pdb( pose, "", true, true, ostream, out_fname, sfr_opts );
}

std::string
PDBPoseOutputter::output_pdb_name( LarvalJob const & job ) const
{
	return ( job.status_prefix() == "" ? "" : ( job.status_prefix() + "_" ) ) + job.job_tag() + "_" +
		ObjexxFCL::lead_zero_string_of( job.nstruct_index(), std::max( 4, int( std::log10( job.nstruct_max() ))) ) +
		".pdb";
}

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
		+ XMLSchemaAttribute( "filename", xs_string )
		+ XMLSchemaAttribute( "filename_pattern", "string_w_one_dollarsign" )
		+ XMLSchemaAttribute( "path", xs_string );
	XMLSchemaComplexTypeGenerator output_pdb;
	output_pdb.element_name( keyname() )
		.complex_type_naming_func( & PoseOutputterFactory::complex_type_name_for_pose_outputter )
		.add_attributes( output_pdb_attributes )
		.write_complex_type_to_schema( xsd );
}

void
PDBPoseOutputter::list_options_read(
	utility::options::OptionKeyList & read_options
)
{
	core::io::StructFileRepOptions::list_options_read( read_options );
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


} // namespace pose_outputters
} // namespace jd3
} // namespace protocols
