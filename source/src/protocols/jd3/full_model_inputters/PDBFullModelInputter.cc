// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/full_model_inputters/PDBFullModelInputter.cc
/// @brief  %PDBFullModelInputter class definition
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Steven Lewis (smlewi@gmail.com)

// Unit headers
#include <protocols/jd3/full_model_inputters/PDBFullModelInputter.hh>
#include <protocols/jd3/full_model_inputters/PDBFullModelInputterCreator.hh>

// Package headers
#include <protocols/jd3/full_model_inputters/FullModelInputSource.hh>
#include <protocols/jd3/full_model_inputters/FullModelInputterFactory.hh>
#include <protocols/jd3/full_model_inputters/full_model_inputter_schemas.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/chemical/ChemicalManager.hh>
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
#include <utility/options/keys/OptionKeyList.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

namespace protocols {
namespace jd3 {
namespace full_model_inputters {


PDBFullModelInputter::PDBFullModelInputter() = default;
PDBFullModelInputter::~PDBFullModelInputter() = default;

bool PDBFullModelInputter::job_available_on_command_line() const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// AMW: this only pertains to stepwise, where -fasta might be the only
	// input. It isn't properly a use of this inputter, but it's kinda all
	// lumped in.
	return ( option[ in::file::s ].user() || option[ in::file::l ].user() || option[ in::file::fasta ].user() ) && !option[ in::file::silent ].user();
}

// There is always one -- EVEN if there is no associated filename!
FullModelInputSources
PDBFullModelInputter::full_model_input_sources_from_command_line()
{
	utility::vector1< utility::file::FileName > filenames_from_command_line = jd2::input_pdb_files_from_command_line();
	FullModelInputSources input_sources;
	input_sources.emplace_back( new FullModelInputSource( keyname() ) );
	input_sources[ 1 ]->input_tag( "S" );

	// 'just fasta' from_scratch case: we must provide a dummy filename
	// that signals we are truly reading in nothing, but the command line
	// writ large is still 'one job'
	if ( filenames_from_command_line.empty() ) {
		input_sources[ 1 ]->store_string_pair( "filename", "fasta" );
	}

	// otherwise, mimic the PoseInputter, but we store all our filenames
	// into a single input source
	for ( utility::file::FileName const & ii_filename : filenames_from_command_line ) {
		input_sources[ 1 ]->store_string_pair( "filename", ii_filename.name() );
	}
	return input_sources;
}

FullModelInputSources
PDBFullModelInputter::full_model_input_sources_from_tag(
	utility::options::OptionCollection const & /*opts*/,
	utility::tag::TagCOP tag
)
{
	// At some point, we should probably make it so you can propose multiple
	// full models... maybe with another set of tags?
	FullModelInputSources input_sources;
	input_sources.emplace_back( new FullModelInputSource( keyname() ) );
	input_sources[ 1 ]->input_tag( "S" );
	if ( tag->hasOption( "filename" ) ) {

		utility::file::FileName fname( tag->getOption< std::string >( "filename" ) );
		if ( tag->hasOption( "path" ) ) {
			fname.path( tag->getOption< std::string >( "path" ) );
		}
		input_sources[ 1 ]->store_string_pair( "filename", fname.name() );

	} else if ( tag->hasOption( "listfile" ) ) {

		std::string list_fname( tag->getOption< std::string >( "listfile" ) );
		utility::io::izstream list_stream( list_fname.c_str() );
		if ( ! list_stream.good() ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Unable to open list file \"" + list_fname + "\"" );
		}
		std::string line;
		while ( getline( list_stream, line ) ) {
			utility::file::FileName fname( line );
			if ( tag->hasOption( "path" ) ) {
				fname.path( tag->getOption< std::string >( "path" ) );
			}
			input_sources[ 1 ]->store_string_pair( "filename", fname.name() );
		}

	} else {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Did not find either a \"filename\" or a \"listfile\" option in the PDB input tag" );
	}
	return input_sources;
}

/// @details This is only a stub implementation -- the logic here ought to be robust enough
/// to handle non-fullatom poses, but it currently isn't.  This will improve.
/// One thought is that the options that control how poses are loaded from disk could
/// be included in the per-job option set, and then passed into this function via an
/// OptionCollection...
core::pose::PoseOP
PDBFullModelInputter::full_model_from_input_source(
	FullModelInputSource const & input_source,
	utility::options::OptionCollection const & options,
	utility::tag::TagCOP /*tag*/ // possibly null-pointing tag pointer
)
{
	// AMW EDIT HERE!!!
	debug_assert( input_source.string_string_map().find( "filename" ) != input_source.string_string_map().end() );
	//core::import_pose::ImportPoseOptions import_opts( options );
	core::chemical::ResidueTypeSetCOP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	utility::vector1< core::pose::PoseOP > input_poses;
	for ( auto const & filename : input_source.string_string_map().at( "filename" ) ) {
		if ( filename == "fasta" ) continue;
		input_poses.push_back( core::import_pose::get_pdb_and_cleanup( filename, rsd_set ) );
	}

	return core::import_pose::initialize_pose_and_other_poses_from_options_and_input_poses( rsd_set, options, input_poses );
}

std::string
PDBFullModelInputter::keyname() { return "PDB"; }

void
PDBFullModelInputter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;

	AttributeList input_pdb_attributes;
	input_pdb_attributes
		+ XMLSchemaAttribute( "filename", xs_string , "XRW TO DO" )
		+ XMLSchemaAttribute( "listfile", xs_string , "XRW TO DO" )
		+ XMLSchemaAttribute( "path", xs_string , "XRW TO DO" );

	full_model_inputter_xsd_type_definition_w_attributes(
		xsd,
		keyname(),
		"XRW TO DO",
		input_pdb_attributes );

}

void
PDBFullModelInputter::list_options_read(
	utility::options::OptionKeyList & read_options
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	core::import_pose::ImportPoseOptions::list_options_read( read_options );
	// The following are a bunch of options called by the full model importing
	// functions in import_pose but not the other ones, i.e., they have nothing
	// to do with ImportPoseOptions
	// Theoretically, we could create a special ImportFullModelOptions class
	// to contain these. When we have finished science and are in search of beauty,
	// do this.
	read_options + in::file::s
		+ in::file::silent
		+ in::file::input_res
		+ full_model::other_poses
		+ full_model::cutpoint_open
		+ in::file::fasta
		+ edensity::mapfile
		+ full_model::extra_min_res
		+ full_model::sample_res
		+ full_model::working_res
		+ full_model::rna::terminal_res
		+ full_model::rna::block_stack_above_res
		+ full_model::rna::block_stack_below_res
		+ full_model::root_res
		+ full_model::jump_res
		+ full_model::cutpoint_closed
		+ full_model::cyclize
		+ full_model::fiveprime_cap
		+ full_model::rna::bulge_res
		+ full_model::extra_min_jump_res
		+ full_model::virtual_sugar_res
		+ full_model::global_seq_file
		+ OptionKeys::stepwise::alignment_anchor_res
		+ full_model::motif_mode
		+ full_model::calc_rms_res
		+ full_model::rna::force_syn_chi_res_list
		+ full_model::rna::force_anti_chi_res_list
		+ full_model::rna::force_north_sugar_list
		+ full_model::rna::force_south_sugar_list
		+ full_model::rna::sample_sugar_res
		+ OptionKeys::stepwise::protein::disulfide_file
		+ constraints::cst_file
		+ stepwise::monte_carlo::vary_loop_length_frequency;
}


FullModelInputterOP PDBFullModelInputterCreator::create_inputter() const
{
	return FullModelInputterOP( new PDBFullModelInputter );
}

std::string PDBFullModelInputterCreator::keyname() const
{
	return PDBFullModelInputter::keyname();
}

void PDBFullModelInputterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PDBFullModelInputter::provide_xml_schema( xsd );
}

void PDBFullModelInputterCreator::list_options_read( utility::options::OptionKeyList & read_options ) const
{
	PDBFullModelInputter::list_options_read( read_options );
}

} // namespace full_model_inputters
} // namespace jd3
} // namespace protocols

