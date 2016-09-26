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
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/InnerLarvalJob.hh>
#include <protocols/jd3/pose_outputters/PoseOutputterFactory.hh>

//project headers
#include <core/pose/Pose.hh>
#include <core/io/StructFileRepOptions.hh>
#include <core/io/silent/BinarySilentStruct.hh>

// ObjexxFCL
#include <ObjexxFCL/string.functions.hh>

// basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

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

SilentFilePoseOutputter::SilentFilePoseOutputter() {}
SilentFilePoseOutputter::~SilentFilePoseOutputter() {}

bool
SilentFilePoseOutputter::outputter_specified_by_command_line()
{
	return basic::options::option[ basic::options::OptionKeys::out::file::silent ].user();
}

void
SilentFilePoseOutputter::determine_job_tag(
	utility::tag::TagCOP output_tag,
	utility::options::OptionCollection const & /*job_options*/,
	InnerLarvalJob & job
) const {
	if ( output_tag ) {
		using namespace utility::tag;
		runtime_assert( output_tag->hasTag( keyname() ) );
		TagCOP silent_tag = output_tag->getTag( keyname() );
		if ( silent_tag->hasOption( "filename" ) ) {
			utility::file::FileName fname( silent_tag->getOption< std::string >( "filename" ));
			job.job_tag( fname.base() );
		} else if ( silent_tag->hasOption( "filename_pattern" )) {
			// replace the single "$" in the filename_pattern with the input_tag already set for the inner larval job
			// we know there is only a single dollar sign because the "string_w_one_dollarsign" restriction is
			// formulated in such a way that forbids multiple dollar signs and requires there to be at least one
			std::string pattern = silent_tag->getOption< std::string >( "filename_pattern" );
			std::string job_tag = utility::replace_in( pattern, "$", job.input_tag() );
			job.job_tag( job_tag );
		}
	} else {
		job.job_tag( job.input_tag() );
	}
}

/// @details If this is just being given an InnerLarvalJob, that means that the
/// outputter is at most a function of the input source--as there's a bijection
/// b/n ILJ and input source. So, let's just use the job_tag. What else can we do?
std::string
SilentFilePoseOutputter::outputter_for_job(
	utility::tag::TagCOP,
	utility::options::OptionCollection const &,
	InnerLarvalJob const & ilj
) const
{
	// Hope determine_job_tag has already been called?
	return ilj.job_tag();
}

bool SilentFilePoseOutputter::job_has_already_completed( LarvalJob const & /*job*/ ) const
{
	// STUBBED OUT!
	return false;
}


void SilentFilePoseOutputter::mark_job_as_having_started( LarvalJob const & /*job*/ ) const
{
	// STUBBED OUT!
}

std::string
SilentFilePoseOutputter::class_key() const
{
	return keyname();
}

void
SilentFilePoseOutputter::write_output_pose(
	LarvalJob const & job,
	utility::options::OptionCollection const & /*job_options*/,
	core::pose::Pose const & pose )
{
	std::string out_fname = output_silent_name( job );

	//StructFileRepOptionsOP sfr_opts( new StructFileRepOptions( job_options ) );

	// What kind of silent struct? Figure out here.
	// For now, assume binary--it'll never fail.

	// Just add pose to buffer!
	buffered_structs_.push_back( std::make_pair( BinarySilentStructOP( new BinarySilentStruct ), job.job_tag() ) );
	buffered_structs_[ buffered_structs_.size() ].first->fill_struct( pose, out_fname );
	//core::io::pdb::dump_pdb( pose, "", true, true, ostream, out_fname, sfr_opts );
}

std::string
SilentFilePoseOutputter::output_silent_name( LarvalJob const & job ) const
{
	return ( job.status_prefix() == "" ? "" : ( job.status_prefix() + "_" ) ) + job.job_tag() + "_" +
		ObjexxFCL::lead_zero_string_of( job.nstruct_index(), std::max( 4, int( std::log10( job.nstruct_max() ))) ) +
		".pdb";
}

void SilentFilePoseOutputter::flush() {

	using utility::vector1;
	using utility::file::FileName;

	typedef std::map< std::string, SilentFileData > SFD_MAP;
	typedef vector1< std::pair< SilentStructOP, FileName > >::iterator iter;
	SFD_MAP sfds;

	// There'll be an option read that involves disabling writing any structures;
	// look for that -- in jd2 it set bWriteNoStructures_

	//tr.Debug << "writing " << buffered_structs_.size() << " structs." << std::endl;
	for ( iter it = buffered_structs_.begin(), end = buffered_structs_.end(); it != end; ++it ) {
		//tr.Debug << "writing struct " << ss->decoy_tag() << std::endl;
		//tr.Debug << "writing struct " << (*it->first)->decoy_tag() << std::endl;
		//SilentStructOP ss = it->first;
		sfds[ it->second ].add_structure( (*it->first) );
	}
	for ( SFD_MAP::iterator it = sfds.begin(); it!=sfds.end(); ++it ) {
		it->second.write_all( it->first );
	}
	// very important to clear after writing!
	buffered_structs_.clear();

	//tr.Debug << "currently have " << saved_structs_.size() << " structs." << std::endl;
}

std::string
SilentFilePoseOutputter::keyname() { return "SilentFile"; }

void
SilentFilePoseOutputter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;

	XMLSchemaRestriction string_w_one_dollarsign;
	string_w_one_dollarsign.name( "string_w_one_dollarsign" );
	string_w_one_dollarsign.base_type( xs_string );
	// this reads "anything but a dollar sign or nothing followed by a dollar sign followed by anything but a dollar sign or nothing
	string_w_one_dollarsign.add_restriction( xsr_pattern, "[^$]*$[^$]*" );
	xsd.add_top_level_element( string_w_one_dollarsign );

	AttributeList output_silent_attributes;
	output_silent_attributes
		+ XMLSchemaAttribute( "filename", xs_string )
		+ XMLSchemaAttribute( "filename_pattern", "string_w_one_dollarsign" )
		+ XMLSchemaAttribute( "path", xs_string );
	XMLSchemaComplexTypeGenerator output_silent;
	output_silent.element_name( keyname() )
		.complex_type_naming_func( & PoseOutputterFactory::complex_type_name_for_pose_outputter )
		.add_attributes( output_silent_attributes )
		.write_complex_type_to_schema( xsd );
}

void
SilentFilePoseOutputter::list_options_read(
	utility::options::OptionKeyList & read_options
)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	StructFileRepOptions::list_options_read( read_options );
	read_options +
		+ out::silent_gz
		+ out::file::atom_tree_diff
		+ out::file::atom_tree_diff_bb
		+ out::file::atom_tree_diff_sc
		+ out::file::atom_tree_diff_bl
		+ out::file::silent
		+ out::file::silent_struct_type
		+ out::file::silent_print_all_score_headers
		+ out::file::raw /*(?)*/
		+ out::file::weight_silent_scores
		+ out::file::silent_preserve_H /*(secretly pdb?)*/
		+ run::no_prof_info_in_silentout;
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
